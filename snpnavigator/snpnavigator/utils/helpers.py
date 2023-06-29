import numpy as np
import pandas as pd

from django.conf import settings

import os
import os.path as path

def get_run_config(run_id):

    df_config = pd.read_csv(path.join(settings.RUNS_DIR, run_id), sep=",")
    #zip config keys and values into dict
    dict_run_config = dict(zip(df_config.key, df_config.value))
    return dict_run_config
