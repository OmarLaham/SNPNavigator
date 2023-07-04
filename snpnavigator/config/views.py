from django.conf import settings
from os import path
from django.contrib.auth.decorators import login_required
from django.shortcuts import render, get_object_or_404, redirect
from django.template import loader
from django.http import HttpResponse, JsonResponse, Http404, HttpResponseRedirect
from django import template
from django.views import View

from django.conf import settings

from django.utils.decorators import method_decorator
from django.views.decorators.csrf import csrf_exempt

import pandas as pd
import numpy as np
import math
import random

import re #regex
import json
import gzip
from glob import glob, escape

#from natsort import natsorted #bette sorting for strings starting with numbers

import shutil

import os
from os import path

from . import helpers
from .helpers import log, LogStatus

gwas_pval_thresh = 5 * (10 ** -8)
peak_open_thresh = 10 #TODO: replace by q-val after peak calling


def filter_snps_in_ocrs(run_id, df_snps, df_peaks, peak_cell_types, peaks_count_matrix_column_names, open_peak_cell_types, close_to_another_ocr): #filter snps using OCRs (Open Chromatin Regions)

    # filter df_peaks to include only peaks that are open for selected cell types in (open_peak_cell_types). This will make iteration faster
    log("filter_snps_in_ocrs", "filter df_peaks", LogStatus.Start)
    df_peaks_filter_query = []
    for col_name in peaks_count_matrix_column_names:
        for open_peak_cell_type in open_peak_cell_types:
            if open_peak_cell_type in col_name:
                df_peaks_filter_query.append("{0} > {1}".format(col_name, peak_open_thresh))
    df_peaks_filter_query = " & ".join(df_peaks_filter_query)
    df_peaks_filtered = df_peaks.query(df_peaks_filter_query)
    log("filter_snps_in_ocrs", "filter df_peaks", LogStatus.End)
    log("len(filtered_peaks)", len(df_peaks_filtered))

    # open mapping_snps_to_peaks file to get IDs of SNPs that lay in the selected open peaks
    # TODO: create auto generation process for the "mapping_snps_to_peaks" file while preparing run.
    log("filter_snps_in_ocrs", "read mapping_snps_to_peaks_file", LogStatus.Start)
    snp_peaks_mapping_file_path = path.join(settings.MEDIA_ROOT, "runs", run_id, "auto_generated_files", "sz_mapping_snps_to_peaks.tsv")
    df_mapping_snps_to_peaks = pd.read_csv(snp_peaks_mapping_file_path,
                                           sep="\t",
                                           keep_default_na=False) #keep_default_na=False is important so empty string or "NA" strings are not converted to pd.NaN)

    #keep only SNPs laying in open OCRs of selected cell types
    df_mapping_snps_to_peaks = df_mapping_snps_to_peaks[df_mapping_snps_to_peaks["peak_name"].isin(df_peaks_filtered["name"].values.tolist())]

    # filter to SNPs that lay in OCR and close to another OCR
    if close_to_another_ocr == 1:
        log("filter_snps_in_ocrs", "filter close to another OCR", LogStatus.Start)
        df_mapping_snps_to_peaks = df_mapping_snps_to_peaks.query("adj_prev_peak_name != '' | adj_next_peak_name != ''")
        log("filter_snps_in_ocrs", "filter close to another OCR", LogStatus.End)

    log("filter_snps_in_ocrs", "read mapping_snps_to_peaks_file", LogStatus.End)

    # filter df_snps using IDs intersection
    log("filter_snps_in_ocrs", "filter df_snps", LogStatus.Start)
    df_filtered_snps = df_snps[df_snps["id"].isin(df_mapping_snps_to_peaks["snp_id"].values.tolist())]
    log("filter_snps_in_ocrs", "filter df_snps", LogStatus.End)

    return df_filtered_snps

def json_snp_query(request, run_id, spec_chr, open_peak_cell_types, cpg_island, close_to_another_ocr, diseases_peaks_match, diseases_peaks_mismatch):

    dict_run_config = helpers.get_run_config(run_id)

    gwas_pval_col = dict_run_config["condition_1_pval_col"]

    log("query", "load GWAS", LogStatus.Start)
    # TODO: commented till finding a solution to boost speed of loading df GWAS
    #df_gwas = pd.read_csv(path.join(settings.MEDIA_ROOT, "data", "gwas", dict_run_config["condition_1_name"], dict_run_config["condition_1_gwas_file"]),
    #                      sep="\t", skiprows=int(dict_run_config["condition_1_gwas_file_skiprows"]))
    # European SZ GWAS will be loaded on server startup as temp condition1 GWAS file
    df_gwas = settings.DF_GWAS
    log("query", "load GWAS", LogStatus.End)

    #rename some cols of gwas df for standardization
    df_gwas = df_gwas.rename(columns={
        dict_run_config["condition_1_pval_col"]: "pval",
        dict_run_config["condition_1_snp_id_col"]: "id",
        dict_run_config["condition_1_chrom_col"]: "chr",
        dict_run_config["conditoin_1_pos_col"]: "pos",
    })

    #filter for specific chr if passed
    if spec_chr != "NA":
        df_gwas = df_gwas.query("chr == {0}".format(spec_chr))

    #filter using GWAS pval thresh
    df_gwas = df_gwas.query("pval <= {0}".format(gwas_pval_thresh))

    #filter to only IDs that start with "rs". Some IDs in the database are not standard (e.g. 10:104427825_C_T).
    df_gwas = df_gwas[df_gwas["id"].str.startswith('rs', na=False)]

    #filter SNPs using selected criteria in the UI
    df_selected_snps = df_gwas[["id", "chr", "pos", "pval"]]

    #integrate ATAC-seq peaks if included in queries
    df_peaks = pd.read_csv(path.join(settings.MEDIA_ROOT, "data", "peaks", dict_run_config["condition_1_name"], dict_run_config["condition_1_peaks_file"]),
                          sep="\t", index_col=False)

    # rename some cols of peaks df for standardization
    df_peaks = df_peaks.rename(columns={
        dict_run_config["condition_1_peaks_chrom_col"]: "chr",
        dict_run_config["condition_1_peaks_start_col"]: "start",
        dict_run_config["condition_1_peaks_end_col"]: "end",
        dict_run_config["condition_1_peaks_gene_id_col"]: "geneId",
        dict_run_config["condition_1_peaks_associated_gene_name_col"]: "assoc_gene_name"
    })

    #get peaks cell types and count matrix column names from config
    peak_cell_types = dict_run_config["condition_1_peaks_cell_types"].split(",")# e.g. "GLUT,GABA,OLIG"
    peaks_count_matrix_column_names = dict_run_config["condition_1_count_matrix_column_names"].split(",")# e.g. "AVG_GLUT,AVG_GABA,AVG_OLIG"

    #filter using OCRs of selected cell types if provided
    if open_peak_cell_types != "NA":
        df_selected_snps = filter_snps_in_ocrs(run_id, df_selected_snps, df_peaks, peak_cell_types, peaks_count_matrix_column_names, open_peak_cell_types, close_to_another_ocr)

    #add "Operations" col to use in the UI
    df_selected_snps["operations"] = "<a href='javascript:;'>some link</a>"

    return JsonResponse({"snps": df_selected_snps.values.tolist()})

# def start(request):
#
#     if request.method == 'POST':
#         print("email:", request.POST["email"])
#         print("description:", request.POST["description"])
#         context = {}
#         html_template = loader.get_template('home/start.html')
#         return HttpResponse(html_template.render(context, request))
#     else:
#         context = {}
#         html_template = loader.get_template('home/start.html')
#         return HttpResponse(html_template.render(context, request))
#
# def upload(request, run_id):
#
#     context = {}
#     html_template = loader.get_template('home/upload.html')
#     return HttpResponse(html_template.render(context, request))
