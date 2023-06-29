import numpy as np
import pandas as pd

from django.conf import settings

import os
import os.path as path

def get_run_config(run_id):

    config_file_path = path.join(settings.MEDIA_ROOT, "runs", "{0}.csv".format(run_id))
    print("Config file path:", config_file_path)
    df_config = pd.read_csv(config_file_path, sep=",")
    #zip config keys and values into dict
    dict_run_config = dict(zip(df_config.key, df_config.value))
    return dict_run_config
