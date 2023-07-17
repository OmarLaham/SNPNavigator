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
# TODO: replace by q-val after peak calling or allow user to set count matrix min value to consider OCR (open) or cell-type-specific OCR
# note that according to the paper https://www.nature.com/articles/s41467-020-19319-2 they define a cell-type-specific OCR when the peak is differentially open \
# in a pairwise comparison with all other cell types
peak_open_thresh = 8 # I averaged cols from same cell type -> calculated mins of averages -> average all result mins together!


def filter_snps_in_ocrs(run_id, df_snps, df_peaks, peak_cell_types, peaks_count_matrix_column_names, open_peak_cell_types, cell_specific_ocrs, close_to_another_ocr): #filter snps using OCRs (Open Chromatin Regions)

    # filter df_peaks to include only peaks that are open for selected cell types in (open_peak_cell_types). This will make iteration faster
    log("filter_snps_in_ocrs", "filter df_peaks", LogStatus.Start)
    query_parts = []

    #split str into lst.
    open_peak_cell_types = open_peak_cell_types.split(",")

    if cell_specific_ocrs == "yes": # filter for open cell-type specific OCRs
        for col_name in peaks_count_matrix_column_names:
            for open_peak_cell_type in open_peak_cell_types:
                if open_peak_cell_type in col_name:
                    query_parts.append(open_peak_cell_type)
        df_peaks_filtered = df_peaks.query("specific_for_cell_type == @query_parts")
    elif cell_specific_ocrs == "no": # filter for OCRs regardless of cell-type specificity

        # cell_specific_ocr_thresh = { #TODO: clean if not used
        #     "GLU": 0,
        #     "GABA": 0,
        #     "OLIG": 0,
        #     "MGAS": 0
        # }

        #for col_name in peaks_count_matrix_column_names:
        for open_peak_cell_type in open_peak_cell_types:
            #if open_peak_cell_type in col_name:
            query_parts.append("ocr_{0} == 1".format(open_peak_cell_type))
        df_peaks_filtered = df_peaks.query(" | ".join(query_parts))


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
    # step 1- check for mutation that introduces C (=> COULD introduce CpG if G follows) or removes C (=> Could introduce CpG if G follows) then
    # for each selected mutation, query sequence online to check if G follows
    # step 2- check for mutation that introduces G (=> COULD introduce CpG if C preceeds) or removes G (=> Could introduce CpG if C preceeds) then
    # for each selected mutation, query sequence online to check if C preceeds
    # step 3- return results

    snps_cpg_islands = {}

    # step 1: 'C' mutations
    df_allele_specific_snps = df_selected_snps.query("origin_allele == 'C' | mutation_allele=='C'")

    # for each selected mutation, query sequence online to check if G follows
    for index, row in df_allele_specific_snps.iterrows():
        id = row["id"]
        chrom = row["chr"]
        pos = row["pos"]
        origin_allele = row["origin_allele"]
        mutation_allele = row["mutation_allele"]

        # check if 'G' follows
        # query the UCSC Genome and get xml to check the next nucleotide
        fetch_next_n_seq_url = "http://genome.ucsc.edu/cgi-bin/das/{0}/dna?segment={1}:{2},{3}".format(gwas_genome_version, chrom, pos+1, pos+1)
        response = requests.get(fetch_next_n_seq_url)
        seq_data = xmltodict.parse(response.content)
        next_nucleotide = seq_data["DASDNA"]["SEQUENCE"]["DNA"]["#text"]

        if next_nucleotide.upper() == 'G':
            # decide if introducing or removing a CpG
            if mutation_allele == 'C':
                snps_cpg_islands[id] = "introducing_CpG"
            elif origin_allele == 'C':
                snps_cpg_islands[id] = "removing_CpG"

    # step 2: 'G' mutations
    df_allele_specific_snps = df_selected_snps.query("origin_allele == 'G' | mutation_allele=='G'")

    # for each selected mutation, query sequence online to check if C preceeds
    for index, row in df_allele_specific_snps.iterrows():
        id = row["id"]
        chrom = row["chr"]
        pos = row["pos"]
        origin_allele = row["origin_allele"]
        mutation_allele = row["mutation_allele"]

        # check if 'C' preceeds
        # query the UCSC Genome and get xml to check the next nucleotide
        fetch_prev_n_seq_url = "http://genome.ucsc.edu/cgi-bin/das/{0}/dna?segment={1}:{2},{3}".format(
            gwas_genome_version, chrom, pos - 1, pos - 1)
        response = requests.get(fetch_prev_n_seq_url)
        seq_data = xmltodict.parse(response.content)
        prev_nucleotide = seq_data["DASDNA"]["SEQUENCE"]["DNA"]["#text"]

        if prev_nucleotide.upper() == 'C':
            # decide if introducing or removing a CpG
            if mutation_allele == 'G':
                snps_cpg_islands[id] = "introducing_CpG"
            elif origin_allele == 'G':
                snps_cpg_islands[id] = "removing_CpG"

    cpg_islands_snps_ids = set(snps_cpg_islands.keys())
    df_selected_snps = df_selected_snps[df_selected_snps["id"].isin(cpg_islands_snps_ids)]

    # step 3
    return (df_selected_snps, snps_cpg_islands)

def filter_snps_for_genomic_regions(run_id, df_selected_snps, gwas_genome_version, spec_gen_region):

    # load genes coords for selected genome version
    # TODO: create auto generation processes for the mapping file
    log("filter_snps_for_genomic_regions", "load snps to locus group mapping", LogStatus.Start)
    df_map_snp_locus_group = pd.read_csv(path.join(settings.MEDIA_ROOT, "runs", run_id, "auto_generated_files", "sz_mapping_snps_to_locus_group.tsv"), sep="\t")

    if spec_gen_region == "non-coding":
        df_map_snp_locus_group = df_map_snp_locus_group.query("locus_group!='protein-coding gene'")
    elif spec_gen_region == "protein-coding":
        df_map_snp_locus_group = df_map_snp_locus_group.query("locus_group=='protein-coding gene'")

    log("filter_snps_for_genomic_regions", "load snps to locus group mapping", LogStatus.End)

    selected_snps_ids = set(df_map_snp_locus_group.snp_id.values.tolist())

    df_selected_snps = df_selected_snps[df_selected_snps["id"].isin(selected_snps_ids)]
    return df_selected_snps

def filter_to_match_mismatch_other_condition(run_id, dict_run_config, df_selected_snps, condition_number, condition_match):

    # load condition_2 GWAS
    condition_gwas_sep = dict_run_config["condition_{0}_gwas_file_sep".format(condition_number)]
    log("query", "load condition {0} GWAS".format(condition_number), LogStatus.Start)
    df_gwas_condition = pd.read_csv(
        path.join(settings.MEDIA_ROOT, "runs", run_id, "auto_generated_files", "gwas", dict_run_config["condition_{0}_name".format(condition_number)],
                  dict_run_config["condition_{0}_gwas_file".format(condition_number)]),
        sep=condition_gwas_sep)

    # rename cols to standardize
    condition_pval_col = dict_run_config["condition_{0}_pval_col".format(condition_number)]
    condition_snp_id_col = dict_run_config["condition_{0}_snp_id_col".format(condition_number)]
    df_gwas_condition = df_gwas_condition.rename(columns={
        condition_pval_col: "pval",
        condition_snp_id_col: "id"
    })

    log("query", "load condition {0} GWAS".format(condition_number), LogStatus.End)
    log("query", "filter to match/mismatch condition {0}.".format(condition_number), LogStatus.Start)
    df_gwas_condition = df_gwas_condition.query("pval <= {0}".format(gwas_pval_thresh))
    if condition_match == "match":
        df_selected_snps = df_selected_snps[df_selected_snps["id"].isin(set(df_gwas_condition["id"].values.tolist()))]
    elif condition_match == "mismatch":
        df_selected_snps = df_selected_snps[~df_selected_snps["id"].isin(set(df_gwas_condition["id"].values.tolist()))]
    log("query", "filter to match/mismatch condition {0}.".format(condition_number), LogStatus.End)

    return df_selected_snps

def json_snp_query(request, run_id, spec_chr, spec_gen_region, overlap_eqtl, open_peak_cell_types, cell_specific_ocrs, cpg_island, close_to_another_ocr, condition_2_match, condition_3_match):

    dict_run_config = helpers.get_run_config(run_id)

    gwas_genome_version = dict_run_config["gwas_genome_version"]

    log("query", "load GWAS", LogStatus.Start)
    # TODO: implement creation of  GWAS-pval-thresholded files in auto_generated_files of the run
    df_snps = pd.read_csv(path.join(settings.MEDIA_ROOT, "runs", run_id, "auto_generated_files", "gwas", dict_run_config["condition_1_name"], dict_run_config["condition_1_gwas_file"]),
                         sep=dict_run_config["condition_1_gwas_file_sep"])
    log("query", "load GWAS", LogStatus.End)

    # rename some cols of gwas df for standardization
    df_snps = df_snps.rename(columns={
        dict_run_config["condition_1_pval_col"]: "pval",
        dict_run_config["condition_1_snp_id_col"]: "id",
        dict_run_config["condition_1_chrom_col"]: "chr",
        dict_run_config["condition_1_pos_col"]: "pos",
        dict_run_config["condition_1_allele_origin_col"]: "origin_allele",
        dict_run_config["condition_1_allele_mutation_col"]: "mutation_allele"
    })

    # filter to only IDs that start with "rs". Some IDs in the database are not standard (e.g. 10:104427825_C_T).
    df_snps = df_snps[df_snps["id"].str.startswith('rs', na=False)]

    # filter using GWAS pval thresh
    df_snps = df_snps.query("pval <= {0}".format(gwas_pval_thresh))

    # calc -log10(pval)
    log("query", "calc snps -log10(pval)", LogStatus.Start)
    df_snps["-log10(pval)"] = ""
    df_snps["-log10(pval)"] = df_snps.apply(
        lambda row: -1 * math.log10(row["pval"])
        , axis=1
    )
    log("query", "calc snps -log10(pval)", LogStatus.End)

    # sort by chr then by pos
    log("query", "sorting df snps by chr then pos", LogStatus.Start)
    df_snps = df_snps.sort_values(["chr", "pos"], ascending=True)
    log("query", "sorting df snps by chr then pos", LogStatus.End)

    # store all SNPs ids and df to be able to differentiate selected from unselected SNPs latter
    all_snps_ids = set(df_snps.id.values.tolist())
    df_all_snps = df_snps.copy()

    # filter for specific chr if passed
    if spec_chr != "NA":
        df_snps = df_snps.query("chr == {0}".format(spec_chr))

    # filter for match/mismatch with condition_2
    if condition_2_match != "NA":
        df_snps = filter_to_match_mismatch_other_condition(run_id, dict_run_config, df_snps, 2, condition_2_match)

    if condition_3_match != "NA":
        df_snps = filter_to_match_mismatch_other_condition(run_id, dict_run_config, df_snps, 3, condition_3_match)

    # filter for genomic region if passed
    if spec_gen_region != "NA":
        df_snps = filter_snps_for_genomic_regions(run_id, df_snps, gwas_genome_version, spec_gen_region)


    #keep cols that we really need
    df_snps = df_snps[["id", "chr", "pos", "pval", "eqtl_gene_id", "eqtl_fdr", "origin_allele", "mutation_allele"]]

    # select SNPs that overlap with eQTL (FDR < 0.05) if overlap_eqtl is selected
    if overlap_eqtl != "NA" and overlap_eqtl == "overlap-eqtl":
        log("query", "filter_overlap_eQTL", LogStatus.Start)
        df_snps = df_snps.query("eqtl_fdr != '-'") # NOTE: when I mapped SNPs to eQTLs I mapped only to eQTLs with FDR < 0.05 and if there is no overlapping eQTL I used FDR='-'
        log("query", "filter_overlap_eQTL", LogStatus.End)


    # filter using OCRs of selected cell types if provided
    #df_selected_snps = df_snps.copy()
    snps_cpg_islands = None
    if open_peak_cell_types != "NA":

        # integrate ATAC-seq peaks if included in queries
        df_peaks = pd.read_csv(path.join(settings.MEDIA_ROOT, "data", "peaks", dict_run_config["condition_1_name"],
                                         dict_run_config["condition_1_peaks_file"]),
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
        peak_cell_types = dict_run_config["condition_1_peaks_cell_types"].split(",")  # e.g. "GLUT,GABA,OLIG"
        peaks_count_matrix_column_names = dict_run_config["condition_1_count_matrix_column_names"].split(
            ",")  # e.g. "AVG_GLUT,AVG_GABA,AVG_OLIG"

        log("query", "filter_snps_in_ocrs", LogStatus.Start)
        df_snps = filter_snps_in_ocrs(run_id, df_snps, df_peaks, peak_cell_types, peaks_count_matrix_column_names, open_peak_cell_types, cell_specific_ocrs, close_to_another_ocr)
        log("query", "filter_snps_in_ocrs", LogStatus.End)

        if cpg_island != 0:
            log("query", "filter_snps_for_cpg_islands", LogStatus.Start)
            df_snps, snps_cpg_islands = filter_snps_for_cpg_islands(run_id, df_snps, gwas_genome_version)
            log("query", "filter_snps_for_cpg_islands", LogStatus.End)


    # group snps into chromosomes for Manhattan plot
    # there will be 2 series: selected and not-selected
    log("query", "grouping snps into series for Manhatan plot", LogStatus.Start)

    df_snps = df_snps.reset_index()

    # create Manhattan plot data
    df_snps_manhattan = df_all_snps.copy().reset_index()
    del df_snps_manhattan["index"]
    df_snps_manhattan = df_snps_manhattan.reset_index() # reset index to use the original index as x axis value
    df_snps_manhattan.rename({"index": "x", "-log10(pval)": "y"}, inplace=True)

    # reorder cols so firs two cols are x and y for Highcharts
    df_snps_manhattan = df_snps_manhattan[["index", "-log10(pval)", "id", "pval", "chr", "pos"]]

    # split into 2 manhattan series and keep only x and y cols. Now: Highcharts accepts only numbers and takes only the first 2 values (performance)
    selected_snps_ids = set(df_snps.id.values.tolist())
    unselected_snps_ids = all_snps_ids - selected_snps_ids
    manhattan_series_unselected = df_snps_manhattan[df_snps_manhattan["id"].isin(unselected_snps_ids)][["index", "-log10(pval)"]]
    manhattan_series_selected = df_snps_manhattan[df_snps_manhattan["id"].isin(selected_snps_ids)] # col selection for selected SNPs will be done later after splitting by chr
    # create manhattan_series object for Highcharts and init it with unselected snps - if any -
    manhattan_series = []
    if len(manhattan_series_unselected):
        manhattan_series.append(
            {
                "name": 'Unselected SNPs',
                "id": "Unselected SNPs",
                "marker": {
                    "symbol": 'circle'  # can be e.g. triangle or square
                },
                "color": "lightgray",
                "data": manhattan_series_unselected.values.tolist()
            }
        )

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

    # add "Operations" col to use in the UI.
    df_snps["operations"] = ""
    df_snps["operations"] = df_snps.apply(
        lambda row: "<a class='tbl-snps-lnk' target='_blank' href='https://www.ncbi.nlm.nih.gov/snp/?term={0}'>dbSNP</a>".format(row["id"]) + ", " +
                    "<a class='tbl-snps-lnk' target='_blank' href='https://www.ncbi.nlm.nih.gov/clinvar/?term={0}'>ClinVar</a>".format(row["id"]) + ", " +
                    "<a class='tbl-snps-lnk' target='_blank' href='https://www.ebi.ac.uk/gwas/variants/{0}'>GWAS Catalog</a>".format(row["id"])
        , axis=1
    )
    "<a href='javascript:;'>some link</a>"
    # modify "eQTL Gene ID" col as an HTML link
    df_snps["eqtl_gene_id"] = df_snps.apply(
        lambda row: "<a class='tbl-snps-lnk' target='_blank' href='http://www.ensembl.org/Homo_sapiens/Gene/Summary?g={0}'>{0}</a>".format(row["eqtl_gene_id"])
        , axis=1
    )



    return JsonResponse({
        "selected_snps": df_snps[["id", "chr", "pos", "pval", "eqtl_gene_id", "eqtl_fdr", "operations"]].values.tolist(),
        "snps_cpg_islands": snps_cpg_islands if snps_cpg_islands else [],
        "manhattan": {
            "peaks": [],
            "snps_details": df_snps_manhattan.to_dict('index'),# convert df into dict to query SNPs details while using Manhattan,
            "series": manhattan_series
        }
    })

