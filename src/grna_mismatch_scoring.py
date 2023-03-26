# -*- coding: utf-8 -*-
#!/usr/bin/env python
# Author: Yixin Guo, Ziwei Xue
#
# ------------------------- #
# Python Modules
# ------------------------- #

import numpy as np
import sys
from functools import partial
from collections import Counter

# ------------------------- #
# Package Python Modules
# ------------------------- #

from .dbutils import *
from .grna_query import *

off_target_penalty = {"diffTE":.06, "sameTE":.03, 'exon':.09, 'intron':.09, 'promoter-TSS':.09, "intergenic": .0,}
mm_penalty = {'0': .004, '3': .001, '2':.002, '1':.003}

def map_annotation(te1, te2, x):
    anno = x[0].split('(')[0]
    if te2:
        if anno == "TE":
            if x[2] == te2:
                return (0,"sameTEdup",x[0].split('(')[1].split(';')[0])
            anno = x[1]
            if anno == te1:
                return (off_target_penalty["sameTE"], "sameTE", x[0].split('(')[1].split(';')[0])
            else:
                return (off_target_penalty["diffTE"], "diffTE", x[0].split('(')[1].split(';')[0])
        else:
            return (off_target_penalty[anno], anno, anno)
    else:
        if anno == "TE":
            if anno == te1:
                return (off_target_penalty["sameTE"], "sameTE", x[0].split('(')[1].split(';')[0])
            else:
                return (off_target_penalty["diffTE"], "diffTE", x[0].split('(')[1].split(';')[0])
        else:
            return (off_target_penalty[anno], anno, anno)


def calculate_offtarget_score(gids_mms, self_te_class, self_te_dup=None):
    result = {}
    par_func = partial(map_annotation, self_te_class, self_te_dup)
    for gid, mm in gids_mms.items():
        score = 0
        result[gid] = {}
        result[gid]["mm0"], result[gid]["mm1"], result[gid]["mm2"], result[gid]["mm3"] = {},{},{},{}
        result[gid]["mm0"]["brief"], result[gid]["mm1"]["brief"], result[gid]["mm2"]["brief"], result[gid]["mm3"]["brief"] = [],[],[],[]
        result[gid]["mm0"]["detail"], result[gid]["mm1"]["detail"], result[gid]["mm2"]["detail"], result[gid]["mm3"]["detail"] = [],[],[],[]

        for i in ["mm0", "mm1", "mm2", "mm3"]:
            annotation = list(map(lambda x:(x[1],x[2],x[3]), mm[i]))
            annotation = list(map(par_func, annotation))
            score += sum(list(map(lambda x:x[0], annotation))) * mm_penalty['0']
            result[gid][i]["brief"] =  list(map(lambda x:x[1], annotation))
            result[gid][i]["detail"] = list(map(lambda x:x[2], annotation))

            result[gid][i]["brief"] = dict(Counter(result[gid][i]["brief"]))
            result[gid][i]["detail"] = dict(Counter(result[gid][i]["detail"]))


            gids_mms[gid][i] = {}

            gids_mms[gid][i]["class"] = result[gid][i]["brief"]
            gids_mms[gid][i]["raw"] = result[gid][i]["detail"]

            gids_mms[gid]["score"] = gids_mms[gid]["gscore"] -score
        
    return gids_mms
