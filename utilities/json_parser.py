#!/usr/bin/env python

import json

def json_parser( experiment, entry ):
    e = open(experiment , "r")
    f = open(entry, "r")
    experiment_json = json.load(e)
    entry_json = json.load(f)
    exp_id = experiment_json["experiment_id"]
    p = experiment_json["platforms"]
    if len(p)>1:
        platforms = '+'.join(p)
    else:
        platforms = p[0]
    ent_id = entry_json["entry_id"]
    return exp_id, ent_id, platforms
