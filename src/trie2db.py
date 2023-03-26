# -*- coding: utf-8 -*-
#!/usr/bin/env python
# Author: Ziwei Xue
#
# ------------------------- #
# Python Modules
# ------------------------- #


import numpy as np
import sys
import pandas as pd
import pybedtools
import tqdm

# ------------------------- #
# Package Python Modules
# ------------------------- #

from dbutils import *
from grna_scoring import *
from grna_annotate import *


fa_path = ""
trie_data_path = ""

GROUPBY1 = lambda l: [(k, list([y for (x,y) in l if x == k])) for k in dict(l).keys()]
FLATTEN = lambda x: [i for s in x for i in s]

file_path = ""
sql = SQLite()
def trie2db_1(t,db_path):

    #######################################
    #    This step is Memory Consuming    #
    #######################################

    sql.create_table1(sql.connect(db_path))
    sql.create_table2(sql.connect(db_path))   
    conn =  sql.connect(db_path)
    cursor = conn.cursor()
    print("Building Trie Tree")
    t = trie.build(fa_path)
    print("Building Table2")
    for k,v in zip(t.values(),t.keys()):
        sql.insert_table2_nomm(cursor, (k,v))
    conn.commit()

    ########################
    #    Serialize Trie    #
    # Temporarily del trie #
    ########################
    
    # whether to use gzip?
    with open(trie_data_path, "wb+") as f:
        trie.save(f,t)



    
def trie2db_2(t,db_path):

    conn =  sql.connect(db_path)
    cursor = conn.cursor()
    
    df = pd.read_csv(file_path, sep="\t", header=None)
    df.iloc[:,1:3] = df.iloc[:,1:3].astype(np.int32)
    df.iloc[:,4] = df.iloc[:,4].astype(np.int32)

    df.columns=list(map(str,range(1,9)))
    
    #######################################
    #    This step is Memory Consuming    #
    #######################################

    # START ANNOTATION HERE #
    gRNAbedAnnotated = pd.DataFrame()
    for i in tqdm.trange(0,len(df),10000000):
        gRNAbed = pybedtools.BedTool.from_dataframe(df.iloc[i:i+10000000,:])
        gRNAbedAnnotatedSplit = grna_annotate(gRNAbed)
        gRNAbedAnnotatedSplit.to_csv("Homo_sapiens.GRCh38.97.dna.primary_assembly.gRNA.annotated." + str(i) + ".csv")
        gRNAbedAnnotated = pd.concat([gRNAbedAnnotated, gRNAbedAnnotatedSplit])
    # END ANNOTATION HERE #

    ###################################
    #    This step is time Consuming  #
    ###################################

    print("Calculating initial gRNA score")
    gRNAbedAnnotated.insert(5, "gscore_moreno", np.nan)
    gRNAbedAnnotated.insert(6, "gscore_azimuth", np.nan)

    import warnings
    warnings.filterwarnings("ignore")
    for j in tqdm.trange(0, len(gRNAbedAnnotated),100000):
        df =  gRNAbedAnnotated.iloc[j:j+100000, :]
        rows = []
        seqs = []
        for i in range(len(df)):
            row = list(df.iloc[i,:])
            seq = sql.select_gseq_table2_gid(cursor, row[1])
            seqs.append(seq)
            row[5] = calculate_moreno_mateos(seq, row[2], row[3], row[4])
            rows.append(row)
        az = predict(np.array(list(map(lambda x:rows[x[1]][3][-4:] + x[0] + rows[x[1]][2] + rows[x[1]][4][:3], zip(seqs,range(len(row)))))))
        for i,score in enumerate(az):
            rows[i][6] = score
        for row in rows:
            try:
                row[1] = int(row[1])
                sql.insert_table1(cursor, tuple(row))
            except:
                pass


    gRNAbedAnnotated = np.array(gRNAbedAnnotated)
    inserted_uid = set()
    for i in tqdm.trange(len(gRNAbedAnnotated)):
        i = gRNAbedAnnotated[i]
        k = sql.select_gseq_table2_gid(cursor, i[1])
        uid = i[1]
        if not uid in inserted_uid:
            if i[5].split('(')[0] == 'TE':
                offtargets = dict(GROUPBY1(list(map(lambda x: (x[2], str(t.get_uid(x[0]))), t.get_approximate_hamming(k, 3)))))
                mm1 = json.dumps(offtargets[1]) if 1 in offtargets.keys() else None
                mm2 = json.dumps(offtargets[2]) if 2 in offtargets.keys() else None
                mm3 = json.dumps(offtargets[3]) if 3 in offtargets.keys() else None
                sql.insert_table2_gid(cursor, uid, (mm1, mm2, mm3))
        inserted_uid.add(uid)

    conn.commit()
    conn.close()

def map_function(gids):
    for i in tqdm.trange(len(gids)):
        gid = gids[i]
        k = sql.select_gseq_table2_gid(cursor, gid)
        offtargets = dict(GROUPBY1(list(map(lambda x: (x[2], str(t.get_uid(x[0]))), t.get_approximate_hamming(k, 3)))))
        mm1 = json.dumps(offtargets[1]) if 1 in offtargets.keys() else None
        mm2 = json.dumps(offtargets[2]) if 2 in offtargets.keys() else None
        mm3 = json.dumps(offtargets[3]) if 3 in offtargets.keys() else None
        sql.insert_table2(gid, (mm1,mm2,mm3))
    conn.close()
    return gids

def map_function(fp, gids): 
    conn = sql.connect("Homo_sapiens.GRCh38.97.dna.primary_assembly.gRNA.db") 
    cursor = conn.cursor() 
    for i in tqdm.trange(len(gids)): 
        gid = gids[i] 
        k = sql.select_gseq_table2_gid(cursor, gid) 
        offtargets = dict(GROUPBY1(list(map(lambda x: (x[2], str(t.get_uid(x[0]))), t.get_approximate_hamming(k, 3))))) 
        mm1 = json.dumps(offtargets[1]) if 1 in offtargets.keys() else ''
        mm2 = json.dumps(offtargets[2]) if 2 in offtargets.keys() else ''
        mm3 = json.dumps(offtargets[3]) if 3 in offtargets.keys() else ''
        fp.write( ",".join([gid, mm1,mm2,mm3]) + '\n' ) 
    conn.close() 
    fp.close()