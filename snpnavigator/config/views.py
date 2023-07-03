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

gwas_pval_thresh = 5 * (10 ** -8)


def filter_snps_in_ocrs(df_snps, df_peaks, peak_cell_types, peaks_count_matrix_column_names, open_peak_cell_types): #filter snps using OCRs (Open Chromatin Regions)

    #create set to store SNPs ids
    #for each open_peak_cell_type SNPs that lay in OCR
    #find intersection of SNPs ids
    #filter df_snps using IDs intersection
    #return

    return None
def json_snp_query(request, run_id, spec_chr, open_peak_cell_types, cpg_island, close_to_another_open_peak, diseases_peaks_match, diseases_peaks_mismatch):

    dict_run_config = helpers.get_run_config(run_id)

    gwas_pval_col = dict_run_config["condition_1_pval_col"]

    df_gwas = pd.read_csv(path.join(settings.MEDIA_ROOT, "data", "gwas", dict_run_config["condition_1_name"], dict_run_config["condition_1_gwas_file"]),
                          sep="\t", skiprows=int(dict_run_config["condition_1_gwas_file_skiprows"]))

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
                          sep="\t")

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
        df_selected_snps = filter_snps_in_ocrs(df_selected_snps, df_peaks, peak_cell_types, peaks_count_matrix_column_names, open_peak_cell_types)

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
