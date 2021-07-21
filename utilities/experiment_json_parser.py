#!/usr/bin/env python

import json
import pandas as pd

def json_parser( j , utilitiesPath): 
    f = open(j , "r")
    experiment_json = json.load(f)
    exp_id = experiment_json["experiment_id"]
    category = experiment_json["data_category"]
    sample , library_prep , platform = get_libraries_metadata(experiment_json["libraries"], utilitiesPath)
    software = experiment_json["software"][0]["name"]
    challenge = experiment_json["challenge_id"]
    return exp_id, category, sample, library_prep, platform, software, challenge

def get_libraries_metadata(libraries_list, utilities):
    matrix_file = utilities + "/rnaseq-data-matrix.tsv"
    matrix_read = pd.read_csv(matrix_file , sep="\t")
    sample_list = set()
    library_prep_list = set()
    platform_list = set()
    for l in libraries_list:
        s = matrix_read.loc[matrix_read["file_acc"]==l]["sample"].to_string().split()[1]
        p = matrix_read.loc[matrix_read["file_acc"]==l]["library_prep"].to_string().split()[1]
        t = matrix_read.loc[matrix_read["file_acc"]==l]["platform"].to_string().split()[1]
        sample_list.add(s)
        library_prep_list.add(p)
        platform_list.add(t)
    if len(sample_list)>1:
        sample_name = "sample_kitchen_sink"
    elif len(sample_list)==1:
        sample_name = list(sample_list)[0]
    else:
        print("utilities/rnaseq-data-matrix.tsv not found. Please check that you have in your utilities folder the LRGASP RNA-Seq Data Matrix. You can also find it here: https://lrgasp.github.io/lrgasp-submissions/docs/rnaseq-data-matrix.html .")
    if len(library_prep_list)>1:
        library_prep_name = "library_prep_kitchen_sink"
    elif len(library_prep_list)==1:
        library_prep_name = list(library_prep_list)[0]
    else:
        print("utilities/rnaseq-data-matrix.tsv not found. Please check that you have in your utilities folder the LRGASP RNA-Seq Data Matrix. You can also find it here: https://lrgasp.github.io/lrgasp-submissions/docs/rnaseq-data-matrix.html .")
    if len(platform_list)>1:
        platform_name = "platform_kitchen_sink"
    elif len(platform_list)==1:
        platform_name = list(platform_list)[0]
    else:
        print("utilities/rnaseq-data-matrix.tsv not found. Please check that you have in your utilities folder the LRGASP RNA-Seq Data Matrix. You can also find it here: https://lrgasp.github.io/lrgasp-submissions/docs/rnaseq-data-matrix.html .")
    return sample_name, library_prep_name, platform_name


