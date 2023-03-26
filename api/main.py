from fastapi import FastAPI, Request, Response, HTTPException
from fastapi import responses
from fastapi.responses import PlainTextResponse
from apiutils import PGSQL, te_info, combination_gids, parseGidAnno, parseGidAnno2, query_mismatch_annotation, calculate_offtarget_score
from typing import Counter, Optional, Union, List
import json
from pydantic import BaseModel
from pako import Pako, Js
import numpy as np
from scipy.signal import savgol_filter
import pandas as pd
import matplotlib.pyplot as plt
import os

app = FastAPI(
    title = "CRISPRTE openAPI",
    docs_url="/api/docs", 
    redoc_url="/api/redoc"
)



class TEDupItem(BaseModel):
    te_dup:str 

class TEClassItem(BaseModel):
    te_class:str 

class GidItem(BaseModel):
    gid: Union[int, List[int]]

class GidTEItem(BaseModel):
    gid: int
    te: str
    response: int

class GseqItem(BaseModel):
    gseq: Optional[str]
    gseqs: Optional[List[str]]

class GTFItem(BaseModel):
    type:str
    chrom:str
    start:int
    end:int

class TEItem(BaseModel):
    te_class: str
    te_copy: int

    

@app.post("/api/v3/validCombinationTE")
def getGidByTedup(ga: str):
    return Response(content = Js.btoa(Pako.deflate(json.dumps(list(combination_gids[ga].keys())))))

@app.post("/api/v1/getGidByTedup")
def getGidByTedup(ga: str, item: TEDupItem):
    te_dup = item.te_dup
    if '_dup' in te_dup:
        te_class, te_dup = te_dup.split("_dup")
        te_class = te_info[ga]["te2int"][te_class]
        te_int = int(te_dup)
    else:
        te_class = te_info[ga]["te2int"][te_dup]
        te_dup = 0
    pconn = PGSQL.connect(dbname = "crisprte{}".format(ga))
    pcursor = pconn.cursor()
    result = PGSQL.select_gid_table1_tedup(pcursor, te_class, te_dup)
    pconn.close()
    return json.dumps(result)

@app.post("/api/v1/getGidByTeclass")
def getGidByTeclass(ga:str, item: TEClassItem):
    te_class = item.te_class
    te_class = te_info[ga]["te2int"][te_class]
    pconn = PGSQL.connect(dbname = "crisprte{}".format(ga))
    pcursor = pconn.cursor()
    result = PGSQL.select_gid_table1_teclass(pcursor, te_class)
    pconn.close()
    return json.dumps(result)

@app.post("/api/v1/getGidSeqinfoByTedup")
def getGidSeqinfoByTedup(ga:str, item: TEDupItem):
    te_dup = item.te_dup
    if '_dup' in te_dup:
        te_class, te_dup = te_dup.split("_dup")
        te_class = te_info[ga]["te2int"][te_class]
        te_int = int(te_dup)
    else:
        te_class = te_info[ga]["te2int"][te_dup]
        te_dup = 0
    pconn = PGSQL.connect(dbname = "crisprte{}".format(ga))
    pcursor = pconn.cursor()
    result = PGSQL.select_gid_seqinfo_table1_tedup(pcursor, te_class, te_dup)
    result = list(map(lambda x:dict(zip(["Gid", "Location", "Moreno", "Azimuth", "Pam", "Upstream", "Downstream"], x)), result))
    pconn.close()
    return json.dumps(result)

@app.post("/api/v1/getGidSeqinfoByTeclass")
def getGidSeqinfoByTeclass(ga:str, item:TEClassItem):
    te_class = item.te_class
    te_class = te_info[ga]["te2int"][te_class]
    pconn = PGSQL.connect(dbname = "crisprte{}".format(ga))
    pcursor = pconn.cursor()
    result = PGSQL.select_gid_seqinfo_table1_teclass(pcursor, te_class)
    result = list(map(lambda x:dict(zip(["Gid", "Location", "Moreno", "Azimuth", "Pam", "Upstream", "Downstream"], x)), result))
    pconn.close()
    return json.dumps(result)

@app.post("/api/v3/getGRNACombinationByTeclass")
def getGRNACombinationByTeclasss(ga:str, item:TEClassItem):
    pconn = PGSQL.connect(dbname = "crisprte{}".format(ga))
    pcursor = pconn.cursor()
    result = combination_gids[ga][item.te_class]
    gseqs = dict(PGSQL.select_gseq_table2_gids(pcursor, result["gids"]))
    resp = {}
    for k,v in result["comb3"].items():
        gids = json.loads(k)
        resp[json.dumps(list(map(lambda x:gseqs[x], gids)))] = v
    return Response(content = Js.btoa(Pako.deflate(json.dumps(resp))))



@app.post("/api/v2/getGRNAByTedup")
def getGRNAByTedup(ga:str, item: TEDupItem):
    te_dup = item.te_dup
    if '_dup' in te_dup:
        te_class, te_dup = te_dup.split("_dup")
        te_class = te_info[ga]["te2int"][te_class]
        te_int = int(te_dup)
    else:
        te_class = te_info[ga]["te2int"][te_dup]
        te_dup = 0
    pconn = PGSQL.connect(dbname = "crisprte{}".format(ga))
    pcursor = pconn.cursor()
    result = query_mismatch_annotation(pcursor, te_class = te_class, te_dup = te_dup, bydup = True, ga = ga)
    pconn.close()
    return Response(content = Js.btoa(Pako.deflate(json.dumps(result))))

@app.post("/api/v3/getGRNAByTedup")
def getGRNAByTedup(ga:str, item: TEDupItem):
    te_dup = item.te_dup
    if '_dup' in te_dup:
        te_class_, te_dup = te_dup.split("_dup")
        if te_class_ not in te_info[ga]["te2int"].keys():
            response = {
                "err": 0x1,
                "message": "invalid TE class"
            }
            return Response(content = Js.btoa(Pako.deflate(json.dumps(response))))
        te_class = te_info[ga]["te2int"][te_class_]
        te_int = int(te_dup)
    else:
        te_class_ = te_dup
        if te_class_ not in te_info[ga]["te2int"].keys():
            response = {
                "err": 0x1,
                "message": "invalid TE class"
            }
            return Response(content = Js.btoa(Pako.deflate(json.dumps(response))))
        te_class = te_info[ga]["te2int"][te_class_]
        te_dup = 0
    pconn = PGSQL.connect(dbname = "crisprte{}".format(ga))
    pcursor = pconn.cursor()
    result = query_mismatch_annotation(pcursor, te_class = te_class, te_dup = te_dup, bydup = True, ga = ga)
    result = calculate_offtarget_score(result, te_class_, te_class_ + "_dup" + str(te_dup))
    pconn.close()
    return Response(content = Js.btoa(Pako.deflate(json.dumps(result))))

@app.post("/api/v1/getGidAnno")
def getGidAnno(ga:str, item:GidItem):
    gid = item.gid
    pconn = PGSQL.connect(dbname = "crisprte{}".format(ga))
    pcursor = pconn.cursor()
    if type(gid) == int:
        result = PGSQL.select_gid_anno_table1_gid(pcursor, gid)
    else:
        result = PGSQL.select_gid_anno_table1_gids(pcursor, gid)
    result = list(map(lambda x:dict(zip(["Gid", "AnnoClass", "AnnoInfo", "Dup"], parseGidAnno(ga, x))), result))
    pconn.close()
    return json.dumps(result)

@app.post("/api/v1/getGidInfo")
def getGidInfo(ga:str, item:GidItem) -> str:
    gid = item.gid
    pconn = PGSQL.connect(dbname = "crisprte{}".format(ga))
    pcursor = pconn.cursor()
    result = PGSQL.select_allinfo_table1_gid(pcursor, gid)
    result = list(map(lambda x:dict(zip(["Location", "Gid", "Moreno", "Azimuth", "Pam", "Upstream", "Downstream", "AnnoClass", "AnnoInfo", "Dup"], parseGidAnno(ga, x))), result))
    pconn.close()
    return json.dumps(result)

@app.post("/api/v1/getGseqMismatch")
def getGseqMismatch(ga:str, item:GidItem) -> str:
    gid = item.gid
    pconn = PGSQL.connect(dbname = "crisprte{}".format(ga))
    pcursor = pconn.cursor()
    result = PGSQL.select_table2_gid(pcursor, gid)
    result = dict(zip(["Gid", "Gseq", "Mismatch1", "Mismatch2", "Mismatch3"], result))
    pconn.close()
    return json.dumps(result)

@app.post("/api/v1/getGseq")
def getGseqMismatch(ga:str, item: GidItem, with_mismatch: bool = False) -> str:
    gid = item.gid
    pconn = PGSQL.connect(dbname = "crisprte{}".format(ga))
    pcursor = pconn.cursor()
    if type(gid) == int:
        if with_mismatch:
            result = PGSQL.select_gseq_table2_where_mm(pcursor, gid)
        else:
            result = PGSQL.select_gseq_table2_gid(pcursor, gid)
        pconn.close()
        return result

    else:
        if with_mismatch:
            result = PGSQL.select_gseq_table2_gids(pcursor, gid)
        else:
            result = PGSQL.select_gseqs_table2_where_mm(pcursor, gid)
        result = list(map(lambda x:dict(zip(["Gid", "Gseq"], parseGidAnno(ga, x))), result))
        pconn.close()
        return result

@app.post("/api/v2/getGid")
def getGid(ga:str, item: GseqItem):
    pconn = PGSQL.connect(dbname = "crisprte{}".format(ga))
    pcursor = pconn.cursor()
    gid = PGSQL.select_gid_gseq(pcursor, item.gseq)
    pconn.close()
    return gid

@app.post("/api/v3/getMismatchBedGseq")
def getMismatchBedGseq(ga:str, item: GseqItem):    

    result = {}
    # Single gseq
    if item.gseq and not item.gseqs:
        pconn = PGSQL.connect(dbname = "crisprte{}".format(ga))
        pcursor = pconn.cursor()
        gid = PGSQL.select_gid_gseq(pcursor, item.gseq)
        if not gid:
            response = {
                "err": 0x1,
                "message": "invalid TE class"
            }
            return Response(content = Js.btoa(Pako.deflate(json.dumps(response))))
        gid,gseq,mm1,mm2,mm3=PGSQL.select_table2_gid(pcursor, gid)
        result['mm0'] = pd.DataFrame(list(map(lambda x:x[0].split("_")[:-1] + ['.','.'] + x[0].split("_")[-1:] + [x[7].tobytes().decode()],PGSQL.select_allinfo_table1_gid(pcursor, gid)))).to_csv(sep="\t",header=None, index=False)
        result['mm1'] = pd.DataFrame(list(map(lambda x:x[0].split("_")[:-1] + ['.','.'] + x[0].split("_")[-1:] + [x[7].tobytes().decode()],PGSQL.select_allinfo_table1_gids(pcursor, mm1)))).to_csv(sep="\t",header=None, index=False)
        result['mm2'] = pd.DataFrame(list(map(lambda x:x[0].split("_")[:-1] + ['.','.'] + x[0].split("_")[-1:] + [x[7].tobytes().decode()],PGSQL.select_allinfo_table1_gids(pcursor, mm2)))).to_csv(sep="\t",header=None, index=False)
        result['mm3'] = pd.DataFrame(list(map(lambda x:x[0].split("_")[:-1] + ['.','.'] + x[0].split("_")[-1:] + [x[7].tobytes().decode()],PGSQL.select_allinfo_table1_gids(pcursor, mm3)))).to_csv(sep="\t",header=None, index=False)
        pconn.close()
        return Response(content = Js.btoa(Pako.deflate(json.dumps(result))))
    # Multiple gseqs
    elif item.gseqs and not item.gseq:
        pconn = PGSQL.connect(dbname = "crisprte{}".format(ga))
        pcursor = pconn.cursor()
        for gseq in item.gseqs:
            print(gseq, len(gseq));
            gid = PGSQL.select_gid_gseq(pcursor, gseq)
            print(gid)
            if not gid:
                continue
            result[gseq] = {}
            df = pd.DataFrame(list(map(lambda x:x[0].split("_")[:-1] + ['.','.'] + x[0].split("_")[-1:] + [x[7].tobytes().decode()],PGSQL.select_allinfo_table1_gid(pcursor, gid))))
            result[gseq]["mm0"] = {
                "target": df.to_csv(sep="\t",header=None, index=False),
                "raw": dict(Counter(list(map(parseGidAnno2, df.iloc[:,6]))))
            }
            _, _, mm1, mm2, mm3 = PGSQL.select_table2_gid(pcursor, gid)
            if len(mm1) > 0:
                df_mm1 = pd.DataFrame(list(map(lambda x:x[0].split("_")[:-1] + ['.','.'] + x[0].split("_")[-1:] + [x[7].tobytes().decode()],PGSQL.select_allinfo_table1_gids(pcursor, mm1))))
                result[gseq]["mm1"] = {
                    "target": df_mm1.to_csv(sep="\t",header=None, index=False),
                    "raw": dict(Counter(list(map(parseGidAnno2, df_mm1.iloc[:,6]))))
                }
            else:
                result[gseq]["mm1"] = {
                    "target": '',
                    "raw": ''
                }
            if len(mm2) > 0:
                df_mm2 = pd.DataFrame(list(map(lambda x:x[0].split("_")[:-1] + ['.','.'] + x[0].split("_")[-1:] + [x[7].tobytes().decode()],PGSQL.select_allinfo_table1_gids(pcursor, mm2))))
                result[gseq]["mm2"] = {
                    "target": df_mm2.to_csv(sep="\t",header=None, index=False),
                    "raw": dict(Counter(list(map(parseGidAnno2, df_mm1.iloc[:,6]))))
                }
            else:
                result[gseq]["mm2"] = {
                    "target": '',
                    "raw": ''
                }

            if len(mm3) > 0:
                df_mm3 = pd.DataFrame(list(map(lambda x:x[0].split("_")[:-1] + ['.','.'] + x[0].split("_")[-1:] + [x[7].tobytes().decode()],PGSQL.select_allinfo_table1_gids(pcursor, mm3))))
                result[gseq]["mm3"] = {
                    "target": df_mm3.to_csv(sep="\t",header=None, index=False),
                    "raw": dict(Counter(list(map(parseGidAnno2, df_mm1.iloc[:,6]))))
                }
            else:
                result[gseq]["mm3"] = {
                    "target": '',
                    "raw": ''
                }
        pconn.close()
        return Response(content = Js.btoa(Pako.deflate(json.dumps(result))))

@app.post("/api/v2/getMismatchBedGid")
def getMismatchBedGid(ga:str, item: GidItem):    
    pconn = PGSQL.connect(dbname = "crisprte{}".format(ga))
    pcursor = pconn.cursor()
    result = []
    gid,gseq,mm1,mm2,mm3=PGSQL.select_table2_gid(pcursor, item.gid)
    result.append(pd.DataFrame(list(map(lambda x:x[0].split("_")[:-1] + ['.','.'] + x[0].split("_")[-1:] + [x[7].tobytes().decode()],PGSQL.select_allinfo_table1_gid(pcursor, gid)))).to_csv(sep="\t",header=None, index=False))
    result.append(pd.DataFrame(list(map(lambda x:x[0].split("_")[:-1] + ['.','.'] + x[0].split("_")[-1:] + [x[7].tobytes().decode()],PGSQL.select_allinfo_table1_gids(pcursor, mm1)))).to_csv(sep="\t",header=None, index=False))
    result.append(pd.DataFrame(list(map(lambda x:x[0].split("_")[:-1] + ['.','.'] + x[0].split("_")[-1:] + [x[7].tobytes().decode()],PGSQL.select_allinfo_table1_gids(pcursor, mm2)))).to_csv(sep="\t",header=None, index=False))
    result.append(pd.DataFrame(list(map(lambda x:x[0].split("_")[:-1] + ['.','.'] + x[0].split("_")[-1:] + [x[7].tobytes().decode()],PGSQL.select_allinfo_table1_gids(pcursor, mm3)))).to_csv(sep="\t",header=None, index=False))
    pconn.close()
    return json.dumps(result)

@app.post("/api/v3/getOnTargetTE")
def getOnTargetTE(ga:str, item: GidTEItem):
    if "{}-{}.svg".format(item.gid, te_info[ga]["te2int"][item.te]) in os.listdir("/root/CRISPRTE/website-dev/static/img/svg/curves") and item.response == 2:
        return "{}-{}.svg".format(item.gid, te_info[ga]["te2int"][item.te])
    pconn = PGSQL.connect(dbname = "crisprte{}".format(ga))
    pcursor = pconn.cursor()
    target_consensus_length = PGSQL.select_te_ttf_by_te(pcursor, te_info[ga]["te2int"][item.te])['length']
    bed1 = PGSQL.select_bed_table1_gid(pcursor, item.gid)
    bed2 = PGSQL.select_gtf_by_tename(pcursor, item.te, True)
    bed3 = bed1.intersect(bed2, wa=True,wb=True).to_dataframe(names=['chrom', 'start', 'end', 'name', 'score', 'strand', 'thickChrom', 'thickStart', 'thickEnd', 'itemRgb', 'blockCount', 'blockStrand', 'blockName','TE'])
    bed3["Relative"] = list(map(lambda z: z["query_begin"] - z["target_begin"] if z and z["query_begin"] - z["target_begin"] > 0 else None, map(lambda x:PGSQL.select_te_dtf_by_te(pcursor, te_info[ga]["te2int"][x.split("_dup")[0]], int(x.split("_dup")[1])) if 'dup' in x else 0, bed3.blockName)))
    relative_distance = list(map(lambda x:x[0] - x[1] + x[4] if x[-1] == '+' else x[3] - x[2] + x[4], zip(bed3["start"], bed3["thickStart"], bed3["end"], bed3["thickEnd"],bed3["Relative"], bed3["blockStrand"])))
    stack = np.zeros(target_consensus_length)
    for i in relative_distance:
        if np.isnan(i):continue
        i = int(i)
        stack[i:i+20] += 1
    w = savgol_filter(stack, 9, 5)
    fig,ax = plt.subplots()
    ax.spines['right'].set_color('none')     
    ax.spines['top'].set_color('none')
    ax.spines['bottom'].set_color('none')     
    ax.spines['left'].set_color('none')
    ax.set_xticks([])
    ax.set_yticks([])
    plt.plot(range(target_consensus_length), w)
    fig.set_size_inches(3,1)
    plt.savefig("/root/CRISPRTE/website-dev/static/img/svg/curves/{}-{}.svg".format(item.gid, te_info[ga]["te2int"][item.te]))
    if item.response == 1:
        return Response(content = Js.btoa(Pako.deflate(json.dumps(bed3.to_csv()))))
    elif item.response == 2:
        return "{}-{}.svg".format(item.gid, te_info[ga]["te2int"][item.te])
    else:
        return None

@app.post("/api/v3/getGtfByRegion")
def getGtfByRegion(ga:str, item:GTFItem):
    pconn = PGSQL.connect(dbname = "crisprte{}".format(ga))
    pcursor = pconn.cursor()
    result = ''
    if item.type == "te":
        result = PGSQL.select_te_gtf_by_region(pcursor, "te", item.chrom, item.start, item.end)
        for i in range(len(result)):
            for key in result[i].keys():
                if type(result[i][key]) == memoryview:
                    result[i][key] = result[i][key].tobytes().decode()
            result[i]['feature'] = 'te'
    elif item.type == "pcg":
        result = PGSQL.select_gtf_by_region(pcursor, "pcg", item.chrom, item.start, item.end)
        for i in range(len(result)):
            for key in result[i].keys():
                if type(result[i][key]) == memoryview:
                    result[i][key] = result[i][key].tobytes().decode()
    elif item.type == "all":
        result_gene = PGSQL.select_gtf_by_region(pcursor, "pcg", item.chrom, item.start, item.end)
        result_te = PGSQL.select_te_gtf_by_region(pcursor, "te", item.chrom, item.start, item.end) 
        for i in range(len(result_gene)):
            for key in result_gene[i].keys():
                if type(result_gene[i][key]) == memoryview:
                    result_gene[i][key] = result_gene[i][key].tobytes().decode()

        for i in range(len(result_te)):
            for key in result_te[i].keys():
                if type(result_te[i][key]) == memoryview:
                    result_te[i][key] = result_te[i][key].tobytes().decode()
            result_te[i]['feature'] = 'te'
        result = result_gene + result_te
    pconn.close()
    return Response(content = Js.btoa(Pako.deflate(json.dumps(result))))

@app.post("/api/v3/getTtfByTE")
def getTtfByTE(ga:str, item: TEItem):
    pconn = PGSQL.connect(dbname = "crisprte{}".format(ga))
    pcursor = pconn.cursor()
    result = PGSQL.select_te_ttf_by_te(pcursor, te_info[ga]["te2int"][TEItem.te_class])
    return Response(content = Js.btoa(Pako.deflate(json.dumps(result))))

@app.post("/api/v3/getDtfByTE")
def getTtfByTE(ga:str, item: TEItem):
    pconn = PGSQL.connect(dbname = "crisprte{}".format(ga))
    pcursor = pconn.cursor()
    result = PGSQL.select_te_dtf_by_te(pcursor, te_info[ga]["te2int"][TEItem.te_class], item.te_copy)
    return Response(content = Js.btoa(Pako.deflate(json.dumps(result))))

