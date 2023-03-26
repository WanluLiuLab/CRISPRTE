from .dbutils import *
import pickle

with open("/root/Data/CRISPRTE/refData/te2int-hg38.pkl","rb") as f:
    te2int = pickle.load(f)
    int2te = {i:j for j, i in te2int.items()}

with open("/root/Data/CRISPRTE/refData/TE-copy-hg38.pkl","rb") as f:
    te_info = pickle.load(f)

with open("/root/Data/CRISPRTE/refData/te2int-mm10.pkl", "rb") as f: 
    te2int = pickle.load(f)
    int2te = {i:j for j, i in te2int.items()}

with open("/root/Data/CRISPRTE/refData/TE-copy-mm10.pkl","rb") as f:
    te_info = pickle.load(f)


def insert_table1(cursor, result):
    result = list(result)
    t = result[7].split("(")[0]
    if t == "TE":
        if result[8] not in te2int.keys():
            return None
        if len(result[9].split("_dup")) < 2:
            result[8] = te2int[result[9].split("_dup")[0]]
            result[9] = 0
        else:
            result[8] = te2int[result[9].split("_dup")[0]]
            result[9] = int(result[9].split("_dup")[1])
    elif t == "intron":
        result[8], result[9] = -1,0
    elif t == 'exon':
        result[8], result[9] = -2,0
    elif t == 'intergenic':
        result[8], result[9] = -3,0
    elif t == 'promoter-TSS':
        result[8], result[9] = -4,0
    else:
        raise ValueError("Unknown datatype")
    result = tuple(result)
    PGSQL.insert_table1(cursor, result)

def insert_table1_remedy(cursor, result):
    result = list(result)
    if result[7].split("(")[0] == "TE":
        if result[8] not in te2int.keys():
            return
        if len(result[9].split("_dup")) < 2:
            return (te2int[result[9]], 0)
        result[8] = te2int[result[9].split("_dup")[0]]
        result[9] = int(result[9].split("_dup")[1])
    else:
        return
    result = tuple(result)
    PGSQL.insert_table1(cursor, result)

def parseIfNone(r):
    if not r or r == "None":
        return '{' + '}'
    else:
        return  '{' + ",".join(json.loads(r)) + '}'

def insert_table2(cursor, result):
    if result[2] == result[3] == result[4] == None or result[2] == result[3] == result[4] == "None":
        result = (result[0], result[1], '{' + '}', '{' + '}', '{' + '}')
    else:
        result = (result[0], result[1], parseIfNone(result[2]), parseIfNone(result[3]), parseIfNone(result[4]))
    PGSQL.insert_table2(cursor, result)


