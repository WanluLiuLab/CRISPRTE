# -*- coding: utf-8 -*-
#!/usr/bin/env python
# Author: Ziwei Xue
#
# ------------------------- #
# Python Modules
# ------------------------- #


import numpy as np
import random

from .dbutils import *

random.seed(990211)
np.random.seed(990211)

def FLATTEN(x): return [i for s in x for i in s]

def GROUPBY2(iterable):
    result = {}
    for i in iterable:
        if i[-1] in result.keys():
            result[i[-1]] = result[i[-1]] + [i[:-1]]
        else:
            result[i[-1]] = [i[:-1]]
    return result

with open("./database/hg38_TEs.json","r") as f:
    te_info = json.loads(f.readlines()[0])


def query_mismatch_annotation(cursor, **args):
    mm_max = args["mm_max"] if "mm_max" in args.keys() else 3
    if 'te_dup' in args:
        targetTE_gids = SQLite.select_gid_seqinfo_table1_tedup(cursor, te_dup=args['te_dup'])
    elif 'te_class' in args:
        # print('te_class', te['te_class'])
        targetTE_gids = SQLite.select_gid_seqinfo_table1_teclass(cursor, te_class=args['te_class'])
        # Only shows first 2 grnas per duplicates
        targetTE_gids = GROUPBY2(targetTE_gids)
        for i in targetTE_gids.keys():
            targetTE_gids[i] = targetTE_gids[i][:2]
        targetTE_gids = FLATTEN(targetTE_gids.values())

        targetTE_gids = sorted(targetTE_gids, key=lambda x:x[2])[:100]
    result = {}
    for i in targetTE_gids:
        result[i[0]] = {}
        result[i[0]]['pos'] = i[1]
        result[i[0]]['gscore'] = i[2]
        result[i[0]]['pam'] = i[3]
        result[i[0]]['upstream'] = i[4]
        result[i[0]]['downstream'] = i[5]

    targetTE_gids = [i[0] for i in targetTE_gids]
    for targetTE_gid in targetTE_gids:
        mm0_annotate, mm1_annotate, mm2_annotate, mm3_annotate = [],[],[],[]
        mm_result = SQLite.select_table2_gid(cursor, targetTE_gid)
        mm0_result = SQLite.select_gid_anno_table1_gid(cursor, targetTE_gid)[0]
        mm0_annotate.append(mm0_result)
        if mm_result[2] != None and mm_max > 0:
            mm1_annotate = SQLite.select_gid_anno_table1_gids(cursor, set(mm_result[2]))
        if mm_result[3] != None and mm_max > 1:
            mm2_annotate = SQLite.select_gid_anno_table1_gids(cursor, set(mm_result[3]))
        if mm_result[4] != None and mm_max > 2:
            mm3_annotate = SQLite.select_gid_anno_table1_gids(cursor, set(mm_result[4]))
        
        result[targetTE_gid]['gseq'] = mm_result[1]
        result[targetTE_gid]['mm0'] = mm0_annotate
        result[targetTE_gid]['mm1'] = mm1_annotate
        result[targetTE_gid]['mm2'] = mm2_annotate
        result[targetTE_gid]['mm3'] = mm3_annotate
    # print("result", result)
    return result

# mm0_gid = SQLite.select_gid_table1_tedup

def calculate_mismatch(cursor, gids, mm_max = 1):
    if type(list(gids.values())[0]) == dict:
        for i in gids.items():
            targetTE_gid = i[1]['gid']
            mm0_annotate, mm1_annotate, mm2_annotate, mm3_annotate = [],[],[],[]
            mm_result = SQLite.select_table2_gid(cursor, targetTE_gid)
            if mm_result[2] != None and mm_max > 0:
                for mm1 in set(mm_result[2]):
                    mm1_result = SQLite.select_gid_anno_table1_gid(cursor, mm1)
                    mm1_annotate.append(mm1_result)
            if mm_result[3] != None and mm_max > 1:
                for mm2 in set(mm_result[3]):
                    mm2_result = SQLite.select_gid_anno_table1_gid(cursor, mm2)
                    mm2_annotate.append(mm2_result)
            if mm_result[4] != None and mm_max > 2:
                for mm3 in set(mm_result[4]):
                    mm3_result = SQLite.select_gid_anno_table1_gid(cursor, mm3)
                    mm3_annotate.append(mm3_result)
            mm0_result = SQLite.select_gid_anno_table1_gid(cursor, targetTE_gid)
            mm0_annotate.append(mm0_result)
            gids[i[0]]["mm0"] = mm0_annotate
            gids[i[0]]["mm1"] = mm1_annotate
            gids[i[0]]["mm2"] = mm2_annotate
            gids[i[0]]["mm3"] = mm3_annotate
            gids[i[0]]['gseq'] = mm_result[1]

    elif type(list(gids.values())[0]) == int:
        targetTE_gid = gids['gid']
        mm0_annotate, mm1_annotate, mm2_annotate, mm3_annotate = [],[],[],[]
        mm_result = SQLite.select_table2_gid(cursor, targetTE_gid)
        if mm_result[2] != None:
            for mm1 in set(mm_result[2]):
                mm1_result = SQLite.select_gid_anno_table1_gid(cursor, mm1)
                mm1_annotate.append(mm1_result)
        if mm_result[3] != None:
            for mm2 in set(mm_result[3]):
                mm2_result = SQLite.select_gid_anno_table1_gid(cursor, mm2)
                mm2_annotate.append(mm2_result)
        if mm_result[4] != None:
            for mm3 in set(mm_result[4]):
                mm3_result = SQLite.select_gid_anno_table1_gid(cursor, mm3)
                mm3_annotate.append(mm3_result)
        mm0_result = SQLite.select_gid_anno_table1_gid(cursor, targetTE_gid)
        mm0_annotate.append(mm0_result)
        gids["mm0"] = mm0_annotate
        gids["mm1"] = mm1_annotate
        gids["mm2"] = mm2_annotate            
        gids["mm3"] = mm3_annotate
        gids["gseq"] = SQLite.select_gseq_table2_gid(cursor,targetTE_gid)
    else:
        raise TypeError("Unexpected input type")


def get_combination(db_path, gids_set, mm_max = 1):
    conn = SQLite.connect(db_path)
    cursor = conn.cursor()
    gids = np.unique(gids_set)
    
    gseq_mms = list(map(lambda x:SQLite.select_table2_gid(cursor, x), gids))
    targetTE_gids = {}
    for i in gseq_mms:
        targetTE_gids[i[0]] = {}
        mm1_annotate, mm2_annotate, mm3_annotate = [],[],[]
        targetTE_gids[i[0]]['gseq'] = i[1]
        targetTE_gids[i[0]]['mm0'] = [SQLite.select_gid_anno_table1_gid(cursor, i[0])]
        if i[2] and mm_max > 0:
            for j in i[2]:
                mm1_annotate.append(SQLite.select_gid_anno_table1_gid(cursor, j))
        if i[3] and mm_max > 1:
            for j in i[3]:
                mm2_annotate.append(SQLite.select_gid_anno_table1_gid(cursor, j))
        if i[4] and mm_max > 2:
            for j in i[4]:
                mm3_annotate.append(SQLite.select_gid_anno_table1_gid(cursor, j))

        targetTE_gids[i[0]]['mm1'] = mm1_annotate
        targetTE_gids[i[0]]['mm2'] = mm2_annotate
        targetTE_gids[i[0]]['mm3'] = mm3_annotate
        targetTE_gids[i[0]]['gscore'] = round(np.mean(list(map(lambda x:x[2],SQLite.select_allinfo_table1_gid(cursor, i[0])))),1)

    return targetTE_gids