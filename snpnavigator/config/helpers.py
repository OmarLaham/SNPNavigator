import numpy as np
import pandas as pd

from django.conf import settings

import os
import os.path as path

#support for enumeration
from enum import Enum

def get_run_config(run_id):

    config_file_path = path.join(settings.MEDIA_ROOT, "runs", run_id, "{0}_config.csv".format(run_id))
    df_config = pd.read_csv(config_file_path, sep=",", keep_default_na=False) #keep_default_na=False is important so empty string or "NA" strings are not converted to pd.NaN creating invalid JSON response
    #zip config keys and values into dict
    dict_run_config = dict(zip(df_config.key, df_config.value))
    return dict_run_config


class LogStatus(Enum):
    Start = "Start"
    End = "End"
def log(main_message, sub_message, status = None):
    print("Debug: > {0}: [{1}] {2}".format(main_message, status if status else "Status NA", sub_message))


def format_eqtl_gene_ids_html(eqtl_gene_id_value):
    if eqtl_gene_id_value == "-":
        return ""

    eqtl_gene_ids = eqtl_gene_id_value.split(",")
    html_generated_links = []

    for eqtl_gene_id in eqtl_gene_ids:
        eqtl_gene_id = eqtl_gene_id.strip()
        html_generated_links.append(
            "<a class='tbl-snps-lnk' target='_blank' href='http://www.ensembl.org/Homo_sapiens/Gene/Summary?g={0}'>{0}</a>".format(
                eqtl_gene_id))

    return ", ".join(html_generated_links)
