from itertools import combinations
import psycopg2
import json
import os
import pickle
from functools import partial
from collections import Counter
import pandas as pd
import pybedtools
from typing import List
from django.conf import settings
# BEGIN In-memory data
import tqdm

from cache_utils.decorators import cached

from .grna_scoring import calculate_azimuth

te_info = {"hg38": {}, "mm10": {}}

filemap = lambda x: os.path.join(settings.BASE_DIR, "api/" "refData", x)
with open(filemap("te2int-hg38.pkl"), "rb") as f:
    te_info["hg38"]["te2int"] = pickle.load(f)
    te_info["hg38"]["int2te"] = {i: j for j,
                                 i in te_info["hg38"]["te2int"].items()}

with open(filemap("TE-copy-hg38.pkl"), "rb") as f:
    te_info["hg38"]["te_info"] = pickle.load(f)

with open(filemap("te2int-mm10.pkl"), "rb") as f:
    te_info["mm10"]["te2int"] = pickle.load(f)
    te_info["mm10"]["int2te"] = {i: j for j,
                                 i in te_info["mm10"]["te2int"].items()}

with open(filemap("TE-copy-mm10.pkl"), "rb") as f:
    te_info["mm10"]["te_info"] = pickle.load(f)

combination_gids = {}

with open(filemap("hg38_combination_gids.json")) as f:
    combination_gids["hg38"] = json.load(f)

with open(filemap("mm10_combination_gids.json")) as f:
    combination_gids["mm10"] = json.load(f)

# END In-memory data

class PGSQL:
    @staticmethod
    def connect(dbname):
        return psycopg2.connect(f"dbname={'crisprtehg38' if form['ga'] == 'hg38' else 'crisprtemm10'} user=postgres port=5432")

    @staticmethod
    def select_gid_gseq(cursor, gseq: str) -> int:
        cursor.execute(f"SELECT gid FROM table2 WHERE gseq = '{gseq}'")
        result = cursor.fetchone()
        return None if not result or len(result) == 0 else result[0]

    @staticmethod
    def select_gid_table1_tedup(cursor, te_class: int, te_dup: int) -> list:
        cursor.execute(f"SELECT gid FROM table1 WHERE te_class = {te_class} AND te_dup = {te_dup}")

        result = list(map(lambda x: x[0], cursor.fetchall()))
        return result

    @staticmethod
    def select_gid_table1_teclass(cursor, te_class: int):
        cursor.execute(f"SELECT gid FROM table1 WHERE te_class = {te_class}")
        result = list(map(lambda x: x[0], cursor.fetchall()))
        return result

    @staticmethod
    def select_gid_seqinfo_table1_tedup(cursor, te_class: int, te_dup: int):
        cursor.execute(f"SELECT gid, pos, gscore_moreno, gscore_azimuth, pam, upstream, downstream FROM table1 WHERE te_class = {te_class} AND te_dup = {te_dup}")

        result = list(cursor.fetchall())
        return result

    @staticmethod
    def select_gid_seqinfo_table1_teclass(cursor, te_class: int):
        cursor.execute(f"SELECT gid,pos,gscore_moreno,gscore_azimuth,pam,upstream,downstream,te_dup FROM table1 WHERE te_class = {te_class}")

        result = list(cursor.fetchall())
        return result

    @staticmethod
    def select_gid_anno_table1_gid(cursor, gid: int):
        cursor.execute(f"SELECT gid, anno_class, te_class, te_dup FROM table1 WHERE gid = {gid}")

        result = list(cursor.fetchall())
        return result

    @staticmethod
    def select_gid_anno_table1_gids(cursor, gids: list):
        cursor.execute(f"SELECT gid, anno_class, te_class, te_dup FROM table1 WHERE gid IN ({','.join(map(str, gids))})")

        result = list(cursor.fetchall())
        return result

    @staticmethod
    def select_anno_table1_gid(cursor, gid: int):
        cursor.execute(f"SELECT anno_class, te_class, te_dup FROM table1 WHERE gid = {gid}")

        result = list(cursor.fetchall())
        return result

    @staticmethod
    def select_allinfo_table1_gid(cursor, gid: int):
        cursor.execute(f"SELECT pos, gid, gscore_moreno,gscore_azimuth, pam, upstream, downstream, anno_class, te_class, te_dup FROM table1 WHERE gid = {gid}")

        result = list(cursor.fetchall())
        return result

    @staticmethod
    def select_allinfo_table1_gids(cursor, gids: List[int]):
        cursor.execute(f"SELECT pos, gid, gscore_moreno,gscore_azimuth, pam, upstream, downstream, anno_class, te_class, te_dup FROM table1 WHERE gid IN ({','.join(list(map(str, gids)))})")

        result = list(cursor.fetchall())
        return result

    @staticmethod
    def select_bed_table1_gid(cursor, gid: int):
        cursor.execute(f"SELECT pos FROM table1 WHERE gid = {gid}")
        result = list(cursor.fetchall())
        result = list(map(lambda x:x[0].split("_"), result))
        result = list(map(lambda x:x[:3] + ['.','.'] + [x[-1]], result))
        return pybedtools.BedTool.from_dataframe(pd.DataFrame(result))

    @staticmethod
    def select_table2_gid(cursor, gid: int):
        cursor.execute(f"SELECT gid, gseq, mm1, mm2, mm3 FROM table2 WHERE gid = {gid}")

        result = list(cursor.fetchall()[0])
        return result

    @staticmethod
    def select_gseq_table2_gid(cursor, gid: int):
        cursor.execute(f"SELECT gseq FROM table2 WHERE gid = {gid}")
        return cursor.fetchall()[0][0]

    @staticmethod
    def select_gseq_table2_gids(cursor, gids: list):
        cursor.execute(f"SELECT gid, gseq FROM table2 WHERE gid IN ({','.join(list(map(str, gids)))})")

        return list(map(lambda x: [x[0], x[1]], cursor.fetchall()))

    @staticmethod
    def select_gseq_table2_where_mm(cursor, gid: int):
        cursor.execute(
            "SELECT gid,mm1,mm2,mm3 FROM table2 WHERE mm1 is not null".format(gid))
        result = cursor.fetchall()
        return result

    @staticmethod
    def select_gseqs_table2_where_mm(cursor, gids: int):
        cursor.execute(f"SELECT gid, gseq FROM table2 WHERE gid IN ({','.join(list(map(str, gids)))}) AND mm1 is not null")

        result = cursor.fetchall()
        return result

    @staticmethod
    def select_gtf_by_genename(cursor, tabname, genename):
        cursor.execute(f"SELECT seqname, source, feature, starting, ending, score, strand, frame, gene_name FROM {tabname} WHERE gene_name='{genename}'")

        result = cursor.fetchall()
        return list(map(lambda x: list(x)[:8] + [{"gene_name": x[8]}], result))

    @staticmethod
    def select_gtf_by_ensembl_geneid(cursor, geneid, tabname="pcg"):
        cursor.execute(f"SELECT seqname, source, feature, starting, ending, score, strand, frame, gene_name FROM {tabname} WHERE gene_id='{geneid}'")

        result = cursor.fetchall()
        return list(map(lambda x: list(x)[:8] + [{"gene_name": x[8]}], result))

    @staticmethod
    def select_gtf_by_tename(cursor, tename, tabname="te"):
        if "dup" in tename:
            cursor.execute(f"SELECT seqname, source, feature, starting, ending, score, strand, frame, gene_id, transcript_id FROM {tabname} WHERE transcript_id = '{tename}'")

        else:
            cursor.execute(f"SELECT seqname, source, feature, starting, ending, score, strand, frame, gene_id, transcript_id FROM {tabname} WHERE gene_id = '{tename}'")

        result = cursor.fetchall()
        return list(map(lambda x: list(x)[:8] + [{"gene_name": x[8]}], result))

    @staticmethod
    def select_gtf_by_region(cursor, tabname, chrom, starting, ending):
        cursor.execute(f"SELECT seqname, source, feature, starting, ending, score, strand, frame, gene_name FROM {tabname} WHERE seqname='{chrom}' AND starting > {starting} AND ending < {ending}")

        result = cursor.fetchall()
        return list(map(lambda x: dict(zip(["seqname", "source", "feature", "start", "end", "score", "strand", "frame", "name"], x)), result))

    @staticmethod
    def select_te_gtf_by_region(cursor, tabname, chrom, starting, ending):
        cursor.execute(f"SELECT seqname, source, feature, starting, ending, score, strand, frame, gene_id, transcript_id FROM {tabname} WHERE seqname='{chrom}' AND starting > {starting} AND ending < {ending}")

        result = cursor.fetchall()
        return list(map(lambda x: dict(zip(["seqname", "source", "feature", "start", "end", "score", "strand", "frame", "gene_id", "name"], x)), result))

    @staticmethod
    def select_gtf_by_tename(cursor, te_name, return_bed):
        cursor.execute(f"SELECT seqname, source, feature, starting, ending, score, strand, frame, gene_id, transcript_id FROM te WHERE gene_id = '{te_name}'")

        result = cursor.fetchall()
        for i in range(len(result)):
            result[i] = list(result[i])
            for j in range(len(result[i])):
                if type(result[i][j]) == memoryview:
                    result[i][j] = result[i][j].tobytes().decode()
        if return_bed:
            return pybedtools.BedTool.from_dataframe(pd.DataFrame(list(map(lambda x: [x[0], x[3], x[4], '.', x[7], x[6], x[-1]], result))))
        return list(map(lambda x:list(x)[:8] + [{"gene_name":x[8]}], result))

    @staticmethod
    def select_te_ttf_by_te(cursor, tid):
        cursor.execute("SELECT name, copies, source, dfamId, classification, clades, description, length, consensus FROM ttf WHERE tid = {}".format(tid))
        result = cursor.fetchall()
        result = dict(zip(["name", "copies", "source", "dfamid", "classification", "clades", "description", "length", "consensus"], result[0]))
        for k,v in result.items():
            if type(v) == memoryview:
                result[k] = v.tobytes().decode()
        return result

    @staticmethod
    def select_te_dtf_by_te(cursor, tid, copy = None):
        if copy != None:
            cursor.execute("SELECT optimal_alignment_score, query_begin, query_end, target_begin, target_end, target_sequence, cigar FROM dtf WHERE tid = {} and dup = {}".format(tid, copy))
            result = cursor.fetchall()
            if len(result) == 0:
                return None 
            result = dict(zip(["optimal_alignment_score", "query_begin", "query_end", "target_begin", "target_end", "target_sequence", "cigar"], result[0]))
            for k,v in result.items():
                if type(v) == memoryview:
                    result[k] = v.tobytes().decode()
        else:
            cursor.execute("SELECT optimal_alignment_score, query_begin, query_end, target_begin, target_end, target_sequence, cigar FROM dtf WHERE tid = {}".format(tid))
            result = cursor.fetchall()
            result = list(map(lambda x:dict(zip(["optimal_alignment_score", "query_begin", "query_end", "target_begin", "target_end", "target_sequence", "cigar"], x)) if len(x) > 0 else None, result))
            for i in range(len(result)):
                for k,v in result[i].items():
                    if type(v) == memoryview:
                        result[i][k] = v.tobytes().decode()
        return result


def parseGidAnno(ga, annotation):
    annotation = list(annotation)
    if annotation[2] == -1:
        annotation[2] = "intron"
    elif annotation[2] == -2:
        annotation[2] = "exon"
    elif annotation[2] == -3:
        annotation[2] = "intergenic"
    elif annotation[2] == -4:
        annotation[2] = "promoter-TSS"
    else:
        annotation[2] = te_info[ga]["int2te"][annotation[2]]
    return [annotation[0], annotation[1].tobytes().decode(), annotation[2], annotation[3]]

def parseGidAnno2(anno):
    if anno.startswith("TE"):
        return anno.split("(")[1].split(";")[0]
    else:
        return anno.split("(")[0]

def FLATTEN(x): return [i for s in x for i in s]


def GROUPBY2(iterable):
    result = {}
    for i in iterable:
        if i[-1] in result.keys():
            result[i[-1]] = result[i[-1]] + [i[:-1]]
        else:
            result[i[-1]] = [i[:-1]]
    return result


off_target_penalty = {"diffTE": .06, "sameTE": .03, 'exon': .09,
                      'intron': .09, 'promoter-TSS': .09, "intergenic": .0, }
mm_penalty = {'0': .004, '3': .001, '2': .002, '1': .003}

@cached(60)
def query_mismatch_annotation(pcursor, **args):
    mm_max = args["mm_max"] if "mm_max" in args.keys() else 3
    int2te = te_info[args["ga"]]["int2te"]
    if 'bydup' in args:
        targetTE_gids = PGSQL.select_gid_seqinfo_table1_tedup(
            pcursor, te_class=args['te_class'], te_dup=args['te_dup'])
    elif 'byclass' in args:
        targetTE_gids = PGSQL.select_gid_seqinfo_table1_teclass(
            pcursor, te_class=args['te_class'])
        targetTE_gids = GROUPBY2(targetTE_gids)
        for i in targetTE_gids.keys():
            targetTE_gids[i] = targetTE_gids[i][:2]
        targetTE_gids = FLATTEN(targetTE_gids.values())
        targetTE_gids = sorted(targetTE_gids, key=lambda x: x[2])[:100]
    result = {}
    sequences = list(map(lambda x: x[5][-4:] + PGSQL.select_gseq_table2_gid(pcursor, x[0]) + x[4] + x[6][:3], targetTE_gids))
    # azimuth_score = [0] * len(sequences)
    azimuth_score = calculate_azimuth(sequences)

    for i,azimuth in zip(targetTE_gids, azimuth_score):
        result[i[0]] = {}
        result[i[0]]['pos'] = i[1]
        result[i[0]]['moreno'] = i[2]
        result[i[0]]['azimuth'] = azimuth
        result[i[0]]['pam'] = i[4]
        result[i[0]]['upstream'] = i[5]
        result[i[0]]['downstream'] = i[6]

    targetTE_gids = [i[0] for i in targetTE_gids]
    pbar = tqdm.tqdm(total=len(targetTE_gids))
    for targetTE_gid in targetTE_gids:
        mm0_annotate, mm1_annotate, mm2_annotate, mm3_annotate = [], [], [], []
        mm_result = PGSQL.select_table2_gid(pcursor, targetTE_gid)
        mm0_annotate = PGSQL.select_gid_anno_table1_gid(pcursor, targetTE_gid)
        for i in range(len(mm0_annotate)):
            mm0_annotate[i] = list(mm0_annotate[i])
            mm0_annotate[i][1] = mm0_annotate[i][1].tobytes().decode()
            annotation = mm0_annotate[i][2]
            if annotation == -1:
                mm0_annotate[i][2] = "intron"
            elif annotation == -2:
                mm0_annotate[i][2] = "exon"
            elif annotation == -3:
                mm0_annotate[i][2] = "intergenic"
            elif annotation == -4:
                mm0_annotate[i][2] = "promoter-TSS"
            else:
                mm0_annotate[i][2] = int2te[mm0_annotate[i][2]]
        if mm_result[2] and mm_max > 0:
            mm1_annotate = PGSQL.select_gid_anno_table1_gids(
                pcursor, list(set(mm_result[2])))
            for i in range(len(mm1_annotate)):
                mm1_annotate[i] = list(mm1_annotate[i])
                mm1_annotate[i][1] = mm1_annotate[i][1].tobytes().decode()
                annotation = mm1_annotate[i][2]
                if annotation == -1:
                    mm1_annotate[i][2] = "intron"
                elif annotation == -2:
                    mm1_annotate[i][2] = "exon"
                elif annotation == -3:
                    mm1_annotate[i][2] = "intergenic"
                elif annotation == -4:
                    mm1_annotate[i][2] = "promoter-TSS"
                else:
                    mm1_annotate[i][2] = int2te[mm1_annotate[i][2]]
        if mm_result[3] and mm_max > 1:
            mm2_annotate = PGSQL.select_gid_anno_table1_gids(
                pcursor, list(set(mm_result[3])))
            for i in range(len(mm2_annotate)):
                mm2_annotate[i] = list(mm2_annotate[i])
                mm2_annotate[i][1] = mm2_annotate[i][1].tobytes().decode()
                annotation = mm2_annotate[i][2]
                if annotation == -1:
                    mm2_annotate[i][2] = "intron"
                elif annotation == -2:
                    mm2_annotate[i][2] = "exon"
                elif annotation == -3:
                    mm2_annotate[i][2] = "intergenic"
                elif annotation == -4:
                    mm2_annotate[i][2] = "promoter-TSS"
                else:
                    mm2_annotate[i][2] = int2te[mm2_annotate[i][2]]
        if mm_result[4] and mm_max > 2:
            mm3_annotate = PGSQL.select_gid_anno_table1_gids(
                pcursor, list(set(mm_result[4])))
            for i in range(len(mm3_annotate)):
                mm3_annotate[i] = list(mm3_annotate[i])
                mm3_annotate[i][1] = mm3_annotate[i][1].tobytes().decode()
                annotation = mm3_annotate[i][2]
                if annotation == -1:
                    mm3_annotate[i][2] = "intron"
                elif annotation == -2:
                    mm3_annotate[i][2] = "exon"
                elif annotation == -3:
                    mm3_annotate[i][2] = "intergenic"
                elif annotation == -4:
                    mm3_annotate[i][2] = "promoter-TSS"
                else:
                    mm3_annotate[i][2] = int2te[mm3_annotate[i][2]]
        result[targetTE_gid]['gseq'] = mm_result[1]
        result[targetTE_gid]['mm0'] = mm0_annotate
        result[targetTE_gid]['mm1'] = mm1_annotate
        result[targetTE_gid]['mm2'] = mm2_annotate
        result[targetTE_gid]['mm3'] = mm3_annotate
        pbar.update()
    pbar.close()
    return result


def map_annotation(te1, te2, x):
    anno = x[0].split('(')[0]
    if te2:
        if anno == "TE":
            if x[2] == te2:
                return (0, "sameTEdup", x[0].split('(')[1].split(';')[0])
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
        result[gid]["mm0"], result[gid]["mm1"], result[gid]["mm2"], result[gid]["mm3"] = {
        }, {}, {}, {}
        result[gid]["mm0"]["brief"], result[gid]["mm1"]["brief"], result[gid]["mm2"]["brief"], result[gid]["mm3"]["brief"] = [], [], [], []
        result[gid]["mm0"]["detail"], result[gid]["mm1"]["detail"], result[gid]["mm2"]["detail"], result[gid]["mm3"]["detail"] = [], [], [], []

        for i in ["mm0", "mm1", "mm2", "mm3"]:
            annotation = list(map(lambda x: (x[1], x[2], x[3]), mm[i]))
            annotation = list(map(par_func, annotation))
            score += sum(list(map(lambda x: x[0],
                                  annotation))) * mm_penalty['0']
            result[gid][i]["brief"] = list(map(lambda x: x[1], annotation))
            result[gid][i]["detail"] = list(map(lambda x: x[2], annotation))

            result[gid][i]["brief"] = dict(Counter(result[gid][i]["brief"]))
            result[gid][i]["detail"] = dict(Counter(result[gid][i]["detail"]))

            gids_mms[gid][i] = {}

            gids_mms[gid][i]["class"] = result[gid][i]["brief"]
            gids_mms[gid][i]["raw"] = result[gid][i]["detail"]
            gids_mms[gid]["score"] = score

    return gids_mms
