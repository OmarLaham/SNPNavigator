from django.conf import settings
from os import path
from django.contrib.auth.decorators import login_required
from django.shortcuts import render, get_object_or_404, redirect
from django.template import loader
from django.http import HttpResponse, JsonResponse, Http404, HttpResponseRedirect
from django import template
from django.views import View
from planet.models import Run

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

from natsort import natsorted #bette sorting for strings starting with numbers

import shutil

import os
from os import path

import utils.helpers as helpers

def json_snp_query(request, run_id, open_peak_cell_types, cpg_island, close_to_another_open_peak, diseases_peaks_match, diseases_peaks_mismatch):

    dict_run_config = helpers.get_run_config(run_id)

    return JsonResponse({"run_config": dict_run_config})

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
