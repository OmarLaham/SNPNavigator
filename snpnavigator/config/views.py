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


import requests
import xmltodict

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

def filter_snps_for_cpg_islands(run_id, df_selected_snps, gwas_genome_version):

    # Algorithm:
    # step 1- check for mutation that introduces C (=> COULD introduce CpG if G follows) or removes C (=> Could introduce CpG if G follows)
    # step 2- for each selected mutation, query sequence online to check if G follows
    # step 3- return results

    # step 1
    df_selected_snps = df_selected_snps.query("origin_allele == 'C' | mutation_allele=='C'")
    if len(df_selected_snps) == 0:
        return df_selected_snps # We've got 0 matches ;)

    # step 2
    snps_cpg_islands = {}
    for index, row in df_selected_snps.iterrows():
        id = row["id"]
        chrom = row["chr"]
        pos = row["pos"]
        origin_allele = row["origin_allele"]
        mutation_allele = row["mutation_allele"]

        # query the UCSC Genome and get xml to check the next nucleotide
        fetch_seq_url = "http://genome.ucsc.edu/cgi-bin/das/{0}/dna?segment={1}:{2},{3}".format(gwas_genome_version, chrom, pos+1, pos+1)
        response = requests.get(fetch_seq_url)
        seq_data = xmltodict.parse(response.content)
        next_nucleotide = seq_data["DASDNA"]["SEQUENCE"]["DNA"]["#text"]

        if next_nucleotide.upper() == 'G':
            # decide if introducing or removing a CpG
            if mutation_allele == 'C':
                snps_cpg_islands[id] = "introducing_CpG"
            elif origin_allele == 'C':
                snps_cpg_islands[id] = "removing_CpG"

    cpg_islands_snps_ids = set(snps_cpg_islands.keys())
    df_selected_snps = df_selected_snps[df_selected_snps["id"].isin(cpg_islands_snps_ids)]

    # step 3
    return (df_selected_snps, snps_cpg_islands)

def json_snp_query(request, run_id, spec_chr, open_peak_cell_types, cpg_island, close_to_another_ocr, diseases_peaks_match, diseases_peaks_mismatch):

    dict_run_config = helpers.get_run_config(run_id)

    gwas_genome_version = dict_run_config["gwas_genome_version"]
    gwas_pval_col = dict_run_config["condition_1_pval_col"]

    log("query", "load GWAS", LogStatus.Start)
    # TODO: commented till finding a solution to boost speed of loading df GWAS
    #df_gwas = pd.read_csv(path.join(settings.MEDIA_ROOT, "data", "gwas", dict_run_config["condition_1_name"], dict_run_config["condition_1_gwas_file"]),
    #                      sep="\t", skiprows=int(dict_run_config["condition_1_gwas_file_skiprows"]))
    # European SZ GWAS will be loaded on server startup as temp condition1 GWAS file
    df_gwas = settings.DF_GWAS
    log("query", "load GWAS", LogStatus.End)

    # rename some cols of gwas df for standardization
    df_gwas = df_gwas.rename(columns={
        dict_run_config["condition_1_pval_col"]: "pval",
        dict_run_config["condition_1_snp_id_col"]: "id",
        dict_run_config["condition_1_chrom_col"]: "chr",
        dict_run_config["condition_1_pos_col"]: "pos",
        dict_run_config["condition_1_allele_origin_col"]: "origin_allele",
        dict_run_config["condition_1_allele_mutation_col"]: "mutation_allele"
    })

    # filter for specific chr if passed
    if spec_chr != "NA":
        df_gwas = df_gwas.query("chr == {0}".format(spec_chr))

    # filter using GWAS pval thresh
    df_gwas = df_gwas.query("pval <= {0}".format(gwas_pval_thresh))

    # filter to only IDs that start with "rs". Some IDs in the database are not standard (e.g. 10:104427825_C_T).
    df_gwas = df_gwas[df_gwas["id"].str.startswith('rs', na=False)]

    #filter SNPs using selected criteria in the UI
    df_snps = df_gwas[["id", "chr", "pos", "pval", "origin_allele", "mutation_allele"]]

    #calc -log10(pval)
    log("query", "calc snps -log10(pval)", LogStatus.Start)
    df_snps["-log10(pval)"] = ""
    df_snps["-log10(pval)"] = df_snps.apply(
        lambda row: -1 * math.log10(row["pval"])
        , axis=1
    )
    log("query", "calc snps -log10(pval)", LogStatus.End)

    # integrate ATAC-seq peaks if included in queries
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

    # get peaks cell types and count matrix column names from config
    peak_cell_types = dict_run_config["condition_1_peaks_cell_types"].split(",")# e.g. "GLUT,GABA,OLIG"
    peaks_count_matrix_column_names = dict_run_config["condition_1_count_matrix_column_names"].split(",")# e.g. "AVG_GLUT,AVG_GABA,AVG_OLIG"

    # filter using OCRs of selected cell types if provided
    df_selected_snps = df_snps.copy()
    snps_cpg_islands = None
    if open_peak_cell_types != "NA":
        log("query", "filter_snps_in_ocrs", LogStatus.Start)
        df_selected_snps = filter_snps_in_ocrs(run_id, df_snps, df_peaks, peak_cell_types, peaks_count_matrix_column_names, open_peak_cell_types, close_to_another_ocr)
        log("query", "filter_snps_in_ocrs", LogStatus.End)

        if cpg_island != 0:
            log("query", "filter_snps_for_cpg_islands", LogStatus.Start)
            df_selected_snps, snps_cpg_islands = filter_snps_for_cpg_islands(run_id, df_selected_snps, gwas_genome_version)
            log("query", "filter_snps_for_cpg_islands", LogStatus.End)


    # add "selected" col to df_snps, the df containing all gwas sig. snps
    log("query", "adding 'selected' col to df_snps", LogStatus.Start)
    selected_snps_ids = set(df_selected_snps.id.values.tolist())
    df_snps["selected"] = ""
    df_snps["selected"] = df_snps.apply(
        lambda row: True if row["id"] in selected_snps_ids else False
        , axis=1
    )
    log("query", "adding 'selected' col to df_snps", LogStatus.End)

    # group snps into chromosomes for Manhattan plot
    # there will be 2 series: selected and not-selected
    log("query", "grouping snps into series for Manhatan plot", LogStatus.Start)

    # sort by chr then by pos
    log("query", "sorting df snps by chr then pos", LogStatus.Start)
    df_snps = df_snps.sort_values(["chr", "pos"], ascending=True)
    df_snps = df_snps.reset_index()
    log("query", "sorting df snps by chr then pos", LogStatus.End)

    # create Manhattan plot data
    df_snps_manhattan = df_snps.copy()
    del df_snps_manhattan["index"]
    df_snps_manhattan = df_snps_manhattan.reset_index() # reset index to use the original index as x axis value
    df_snps_manhattan.rename({"index": "x", "-log10(pval)": "y"}, inplace=True)

    # reorder cols so firs two cols are x and y for Highcharts
    df_snps_manhattan = df_snps_manhattan[["index", "-log10(pval)", "id", "pval", "chr", "pos", "selected"]]

    # split into 2 manhattan series and keep only x and y cols. Now: Highcharts accepts only numbers and takes only the first 2 values (performance)
    manhattan_series_unselected = df_snps_manhattan[df_snps_manhattan["selected"] == False][["index", "-log10(pval)"]]
    manhattan_series_selected = df_snps_manhattan[df_snps_manhattan["selected"] == True] # col selection for selected SNPs will be done later after splitting by chr

    # create manhattan_series object for Highcharts and init it with unselected snps
    manhattan_series = [
        {
            "name": 'Unselected SNPs',
            "id": "Unselected SNPs",
            "marker": {
                "symbol": 'circle'  # can be e.g. triangle or square
            },
            "color": "lightgray",
            "data": manhattan_series_unselected.values.tolist()
        }
    ]

    # split selected into multiple series so they can get different colors
    manhattan_selected_colors = ["green","blue","yellow","red"]
    manhattan_selected_chroms = manhattan_series_selected["chr"].unique()
    for i in range(len(manhattan_selected_chroms)) :
        chrom = manhattan_selected_chroms[i]
        manhattan_series_selected_chr = manhattan_series_selected.query("chr=={0}".format(chrom))
        manhattan_series.append(
            {
                "name": 'Selected SNPs - Chr{0}'.format(chrom),
                "id": 'Selected SNPs - Chr{0}'.format(chrom),
                "marker": {
                    "symbol": 'circle'  # can be e.g. triangle or square
                },
                "color": manhattan_selected_colors[i % len(manhattan_selected_colors)],
                "data": manhattan_series_selected_chr[["index", "-log10(pval)"]].values.tolist()
            }
        )
    log("query", "grouping snps into series for Manhatan plot", LogStatus.End)

    # add "Operations" col to use in the UI
    df_selected_snps["operations"] = "<a href='javascript:;'>some link</a>"

    return JsonResponse({
        "selected_snps": df_selected_snps.values.tolist(),
        "snps_cpg_islands": snps_cpg_islands if snps_cpg_islands else [],
        "manhattan": {
            "peaks": [],
            "snps_details": df_snps_manhattan.to_dict('index'),# convert df into dict to query SNPs details while using Manhattan,
            "series": manhattan_series
        }
    })

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
