from django.shortcuts import render, get_object_or_404
from django.http import HttpResponse
from django.contrib.auth.decorators import login_required
from django.template import loader
from django.conf import settings
from django import forms
from django.views.decorators.csrf import csrf_exempt
from django.http import JsonResponse
from django.views.static import serve

import os
import json

from .apiutils import PGSQL, te_info, combination_gids, parseGidAnno, parseGidAnno2, query_mismatch_annotation, calculate_offtarget_score
from typing import Counter, Optional, Union, List
import json 
from .pako import Pako, Js
import numpy as np
from scipy.signal import savgol_filter
import pandas as pd
import matplotlib.pyplot as plt
import os
import psycopg2

@csrf_exempt
def v1(request):
    if request.method == 'POST':
        form = json.loads(request.body)
        print(form['key'])
        if form['key'] == 'getGidByTedup':
            ga = form['ga']
            te_dup = form['te_dup']
            if '_dup' in te_dup:
                te_class, te_dup = te_dup.split("_dup")
                te_class = te_info[ga]["te2int"][te_class]
                te_int = int(te_dup)
            else:
                te_class = te_info[ga]["te2int"][te_dup]
                te_dup = 0
            pconn = psycopg2.connect(f"dbname={'crisprtehg38' if form['ga'] == 'hg38' else 'crisprtemm10'} user=postgres port=5432")
            pcursor = pconn.cursor()
            result = PGSQL.select_gid_seqinfo_table1_tedup(pcursor, te_class, te_dup)
            result = list(map(lambda x:dict(zip(["Gid", "Location", "Moreno", "Azimuth", "Pam", "Upstream", "Downstream"], x)), result))
            pconn.close()
            return JsonResponse(result)
        elif form['key'] == 'getGidAnno':
            gid = form['gid']
            pconn = psycopg2.connect(f"dbname={'crisprtehg38' if form['ga'] == 'hg38' else 'crisprtemm10'} user=postgres port=5432")
            pcursor = pconn.cursor()
            if type(gid) == int:
                result = PGSQL.select_gid_anno_table1_gid(pcursor, gid)
            else:
                result = PGSQL.select_gid_anno_table1_gids(pcursor, gid)
            result = list(map(lambda x:dict(zip(["Gid", "AnnoClass", "AnnoInfo", "Dup"], parseGidAnno(ga, x))), result))
            pconn.close()
            return JsonResponse(result)
        elif form['key'] == 'getGidInfo':
            gid = form['gid']
            pconn = psycopg2.connect(f"dbname={'crisprtehg38' if form['ga'] == 'hg38' else 'crisprtemm10'} user=postgres port=5432")
            pcursor = pconn.cursor()
            result = PGSQL.select_allinfo_table1_gid(pcursor, gid)
            result = list(map(lambda x:dict(zip(["Location", "Gid", "Moreno", "Azimuth", "Pam", "Upstream", "Downstream", "AnnoClass", "AnnoInfo", "Dup"], parseGidAnno(ga, x))), result))
            pconn.close()
            return JsonResponse(result)
        elif form['key'] == 'getGseqMismatch':
            gid = form['gid']
            pconn = psycopg2.connect(f"dbname={'crisprtehg38' if form['ga'] == 'hg38' else 'crisprtemm10'} user=postgres port=5432")
            pcursor = pconn.cursor()
            result = PGSQL.select_table2_gid(pcursor, gid)
            result = dict(zip(["Gid", "Gseq", "Mismatch1", "Mismatch2", "Mismatch3"], result))
            pconn.close()
            return JsonResponse(result)
        elif form['key'] == 'getGseq':
            gid = form['gid']
            with_mismatch = form['with_mismatch']
            pconn = psycopg2.connect(f"dbname={'crisprtehg38' if form['ga'] == 'hg38' else 'crisprtemm10'} user=postgres port=5432")
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
                return JsonResponse(result)
        elif form['key'] == 'getMismatchBedGid':
            pconn = PGSQL.connect(dbname = "crisprte{}".format(form['ga']))
            pcursor = pconn.cursor()
            result = []
            gid,gseq,mm1,mm2,mm3=PGSQL.select_table2_gid(pcursor, form['gid'])
            result.append(
                pd.DataFrame(
                    list(
                        map(
                            lambda x : [gseq] + x[0].split("_") [:-1] + ['.','.'] + x[0].split("_") [-1:] + [x[7].tobytes().decode() ],
                            PGSQL.select_allinfo_table1_gid(pcursor, gid)
                        )
                    )
                ).to_csv(sep = "\t", header = None, index = False)
            )
            gseqs = dict(zip(mm1, list(map(lambda x: x[1], PGSQL.select_gseq_table2_gids(pcursor, mm1)))))
            result.append(
                pd.DataFrame(
                    list(
                        map(
                            lambda x : [gseqs[x[1]]] + x[0].split("_") [:-1] + ['.','.'] + x[0].split("_") [-1:] + [x[7].tobytes().decode() ],
                            PGSQL.select_allinfo_table1_gids(pcursor, mm1)
                        )
                    )
                ).to_csv(sep = "\t", header = None, index = False)
            ) 
            gseqs = dict(zip(mm2, dict(zip(mm1, list(map(lambda x: x[1], PGSQL.select_gseq_table2_gids(pcursor, mm2)))))))
            result.append(
                pd.DataFrame(
                    list(
                        map(
                            lambda x : [gseqs[x[1]]] + x[0].split("_") [:-1] + ['.','.'] + x[0].split("_") [-1:] + [x[7].tobytes().decode() ],
                            PGSQL.select_allinfo_table1_gids(pcursor, mm2)
                        )
                    )
                ).to_csv(sep = "\t", header = None, index = False)
            ) 
            gseqs = dict(zip(mm3, dict(zip(mm1, list(map(lambda x: x[1], PGSQL.select_gseq_table2_gids(pcursor, mm3)))))))
            result.append(
                pd.DataFrame(
                    list(
                        map(
                            lambda x : [gseqs[x[1]]] + x[0].split("_") [:-1] + ['.','.'] + x[0].split("_") [-1:] + [x[7].tobytes().decode() ],
                            PGSQL.select_allinfo_table1_gids(pcursor, mm3)
                        )
                    )
                ).to_csv(sep = "\t", header = None, index = False)
            )
            pconn.close()
            return JsonResponse(result)
        
            
            
@csrf_exempt
def v2(request):
    if request.method == 'POST':
        form = json.loads(request.body)
        print(form['key'])
        if form['key'] == 'getGRNAByTedup':
            te_dup = form['te_dup']
            if '_dup' in te_dup:
                te_class, te_dup = te_dup.split("_dup")
                te_class = te_info[ga]["te2int"][te_class]
                te_int = int(te_dup)
            else:
                te_class = te_info[ga]["te2int"][te_dup]
                te_dup = 0
            pconn = psycopg2.connect(f"dbname={'crisprtehg38' if form['ga'] == 'hg38' else 'crisprtemm10'} user=postgres port=5432")
            pcursor = pconn.cursor()
            result = query_mismatch_annotation(pcursor, te_class = te_class, te_dup = te_dup, bydup = True, ga = ga)
            pconn.close()
            return HttpResponse(Js.btoa(Pako.deflate(json.dumps(result))), content_type='application/octet-stream')
        elif form['key'] == 'getGid':
            pconn = psycopg2.connect(f"dbname={'crisprtehg38' if form['ga'] == 'hg38' else 'crisprtemm10'} user=postgres port=5432")
            pcursor = pconn.cursor()
            gid = PGSQL.select_gid_gseq(pcursor, gseq)
            pconn.close()
            return JsonResponse({'gid':gid})


@csrf_exempt
def v3(request):
    if request.method == 'POST':
        form = json.loads(request.body)
        print(form)
        if form['key'] == 'validCombinationTE':
            ga = form['ga']
            return HttpResponse(
                Js.btoa(
                    Pako.deflate(json.dumps(list(combination_gids[ga].keys())))
                ), 
                content_type='application/octet-stream'
            )
        elif form['key'] == 'getGRNACombinationByTeclass':
            ga = form['ga']
            te_class = form['te_class']
            pconn = psycopg2.connect(f"dbname={'crisprtehg38' if form['ga'] == 'hg38' else 'crisprtemm10'} user=postgres port=5432")
            pcursor = pconn.cursor()
            if te_class not in combination_gids[ga].keys():
                response = {
                    "err": 0x1,
                    "message": "invalid TE class"
                }
                return HttpResponse(Js.btoa(Pako.deflate(json.dumps(resp))), content_type='application/octet-stream')
            result = combination_gids[ga][te_class]
            gseqs = dict(PGSQL.select_gseq_table2_gids(pcursor, result["gids"]))
            resp = {}
            for k,v in result["comb3"].items():
                gids = json.loads(k)
                resp[json.dumps(list(map(lambda x:gseqs[x], gids)))] = v
            return HttpResponse(Js.btoa(Pako.deflate(json.dumps(resp))), content_type='application/octet-stream')
        elif form['key'] == 'getGRNAByTedup':
            ga = form['ga']
            te_dup = form['te_dup']
            if '_dup' in te_dup:
                te_class_, te_dup = te_dup.split("_dup")
                if te_class_ not in te_info[ga]["te2int"].keys():
                    response = {
                        "err": 0x1,
                        "message": "invalid TE class"
                    }
                    return HttpResponse(Js.btoa(Pako.deflate(json.dumps(response))), content_type='application/octet-stream')
                te_class = te_info[ga]["te2int"][te_class_]
                te_int = int(te_dup)
            else:
                te_class_ = te_dup
                if te_class_ not in te_info[ga]["te2int"].keys():
                    response = {
                        "err": 0x1,
                        "message": "invalid TE class"
                    }
                    return HttpResponse(Js.btoa(Pako.deflate(json.dumps(response))),content_type='application/octet-stream')
                te_class = te_info[ga]["te2int"][te_class_]
                te_dup = 0
            pconn = psycopg2.connect(f"dbname={'crisprtehg38' if form['ga'] == 'hg38' else 'crisprtemm10'} user=postgres port=5432")
            pcursor = pconn.cursor()
            result = query_mismatch_annotation(pcursor, te_class = te_class, te_dup = te_dup, bydup = True, ga = ga)
            result = calculate_offtarget_score(result, te_class_, te_class_ + "_dup" + str(te_dup))
            pconn.close()
            return HttpResponse(Js.btoa(Pako.deflate(json.dumps(result))),content_type='application/octet-stream')
        elif form['key'] == 'getMismatchBedGseq':
            result = {}
            # Single gseq
            gseq = form.get('gseq', None)
            gseqs = form.get('gseqs', None)
            if gseq and not gseqs:
                pconn = psycopg2.connect(f"dbname={'crisprtehg38' if form['ga'] == 'hg38' else 'crisprtemm10'} user=postgres port=5432")
                pcursor = pconn.cursor()
                gid = PGSQL.select_gid_gseq(pcursor, gseq)
                if not gid:
                    response = {
                        "err": 0x1,
                        "message": "invalid TE class"
                    }
                    return Response(Js.btoa(Pako.deflate(json.dumps(response))))
                gid,gseq,mm1,mm2,mm3=PGSQL.select_table2_gid(pcursor, gid)
                result['mm0'] = pd.DataFrame(
                    list(
                        map(
                            lambda x : [gseq] + x[0].split("_") [:-1] + ['.','.'] + x[0].split("_") [-1:] + [x[7].tobytes().decode() ],
                            PGSQL.select_allinfo_table1_gid(pcursor, gid)
                        )
                    )
                ).to_csv(sep = "\t", header = None, index = False)
                gseqs = dict(zip(mm1, list(map(lambda x: x[1], PGSQL.select_gseq_table2_gids(pcursor, mm1)))))
                result['mm1'] = pd.DataFrame(
                    list(
                        map(
                            lambda x : [gseqs[x[1]]] + x[0].split("_") [:-1] + ['.','.'] + x[0].split("_") [-1:] + [x[7].tobytes().decode() ],
                            PGSQL.select_allinfo_table1_gids(pcursor, mm1)
                        )
                    )
                ).to_csv(sep = "\t", header = None, index = False)
                gseqs = dict(zip(mm2, list(map(lambda x: x[1], PGSQL.select_gseq_table2_gids(pcursor, mm2)))))
                result['mm2'] = pd.DataFrame(
                    list(
                        map(
                            lambda x : [gseqs[x[1]]] + x[0].split("_") [:-1] + ['.','.'] + x[0].split("_") [-1:] + [x[7].tobytes().decode() ],
                            PGSQL.select_allinfo_table1_gids(pcursor, mm2)
                        )
                    )
                ).to_csv(sep = "\t", header = None, index = False)
                gseqs = dict(zip(mm3, list(map(lambda x: x[1], PGSQL.select_gseq_table2_gids(pcursor, mm3)))))
                result['mm3'] = pd.DataFrame(
                    list(
                        map(
                            lambda x : [gseqs[x[1]]] + x[0].split("_") [:-1] + ['.','.'] + x[0].split("_") [-1:] + [x[7].tobytes().decode() ],
                            PGSQL.select_allinfo_table1_gids(pcursor, mm3)
                        )
                    )
                ).to_csv(sep = "\t", header = None, index = False)
                
                pconn.close()
                return HttpResponse(Js.btoa(Pako.deflate(json.dumps(result))),content_type='application/octet-stream')
            # Multiple gseqs
            elif gseqs and not gseq:
                pconn = psycopg2.connect(f"dbname={'crisprtehg38' if form['ga'] == 'hg38' else 'crisprtemm10'} user=postgres port=5432")
                pcursor = pconn.cursor()
                for gseq in gseqs:
                    gid = PGSQL.select_gid_gseq(pcursor, gseq)
                    if not gid:
                        continue
                    result[gseq] = {}
                    df = pd.DataFrame(list(map(lambda x: [gseq] + x[0].split("_")[:-1] + ['.','.'] + x[0].split("_")[-1:] + [x[7].tobytes().decode()],PGSQL.select_allinfo_table1_gid(pcursor, gid))))
                    result[gseq]["mm0"] = {
                        "target": df.to_csv(sep="\t",header=None, index=False),
                        "raw": dict(Counter(list(map(parseGidAnno2, df.iloc[:,7]))))
                    }
                    _, _, mm1, mm2, mm3 = PGSQL.select_table2_gid(pcursor, gid)
                    if len(mm1) > 0:
                        gseqs = dict(zip(mm1, list(map(lambda x: x[1], PGSQL.select_gseq_table2_gids(pcursor, mm1)))))
                        df_mm1 = pd.DataFrame(
                            list(
                                map(
                                    lambda x : [gseqs[x[1]]] + x[0].split("_") [:-1] + ['.','.'] + x[0].split("_") [-1:] + [x[7].tobytes().decode() ],
                                    PGSQL.select_allinfo_table1_gids(pcursor, mm1)
                                )
                            )
                        )
                        result[gseq]["mm1"] = {
                            "target": df_mm1.to_csv(sep="\t",header=None, index=False),
                            "raw": dict(Counter(list(map(parseGidAnno2, df_mm1.iloc[:,7]))))
                        }
                    else:
                        result[gseq]["mm1"] = {
                            "target": '',
                            "raw": ''
                        }
                    if len(mm2) > 0:
                        gseqs = dict(zip(mm2, list(map(lambda x: x[1], PGSQL.select_gseq_table2_gids(pcursor, mm2)))))
                        df_mm2 = pd.DataFrame(
                            list(
                                map(
                                    lambda x : [gseqs[x[1]]] + x[0].split("_") [:-1] + ['.','.'] + x[0].split("_") [-1:] + [x[7].tobytes().decode() ],
                                    PGSQL.select_allinfo_table1_gids(pcursor, mm2)
                                )
                            )
                        )
                        result[gseq]["mm2"] = {
                            "target": df_mm2.to_csv(sep="\t",header=None, index=False),
                            "raw": dict(Counter(list(map(parseGidAnno2, df_mm2.iloc[:,7]))))
                        }
                    else:
                        result[gseq]["mm2"] = {
                            "target": '',
                            "raw": ''
                        }

                    if len(mm3) > 0:
                        gseqs = dict(zip(mm3, list(map(lambda x: x[1], PGSQL.select_gseq_table2_gids(pcursor, mm3)))))
                        df_mm3 = pd.DataFrame(
                            list(
                                map(
                                    lambda x : [gseqs[x[1]]] + x[0].split("_") [:-1] + ['.','.'] + x[0].split("_") [-1:] + [x[7].tobytes().decode() ],
                                    PGSQL.select_allinfo_table1_gids(pcursor, mm3)
                                )
                            )
                        )
                        result[gseq]["mm3"] = {
                            "target": df_mm3.to_csv(sep="\t",header=None, index=False),
                            "raw": dict(Counter(list(map(parseGidAnno2, df_mm3.iloc[:,7]))))
                        }
                    else:
                        result[gseq]["mm3"] = {
                            "target": '',
                            "raw": ''
                        }
                pconn.close()
                return HttpResponse(Js.btoa(Pako.deflate(json.dumps(result))), content_type='application/octet-stream')
        elif form['key'] == 'getOnTargetTE':
            ga = form['ga']
            gid = form['gid']
            te = form['te']
            response = form['response']
            if "{}-{}.svg".format(gid, te_info[ga]["te2int"][te]) in os.listdir("/root/CRISPRTE/website-dev/static/img/svg/curves") and response == 2:
                return "{}-{}.svg".format(gid, te_info[ga]["te2int"][te])
            pconn = psycopg2.connect(f"dbname={'crisprtehg38' if form['ga'] == 'hg38' else 'crisprtemm10'} user=postgres port=5432")
            pcursor = pconn.cursor()
            target_consensus_length = PGSQL.select_te_ttf_by_te(pcursor, te_info[ga]["te2int"][te])['length']
            bed1 = PGSQL.select_bed_table1_gid(pcursor, gid)
            bed2 = PGSQL.select_gtf_by_tename(pcursor, te, True)
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
            plt.savefig("/root/CRISPRTE/website-dev/static/img/svg/curves/{}-{}.svg".format(gid, te_info[ga]["te2int"][te]))
            if response == 1:
                return HttpResponse(Js.btoa(Pako.deflate(json.dumps(bed3.to_csv()))), content_type='application/octet-stream')
            elif response == 2:
                return "{}-{}.svg".format(gid, te_info[ga]["te2int"][te])
            else:
                return None
        elif form['key'] == 'getGtfByRegion':
            types = form['type']
            chrom, start, end = form['chrom'], form['start'], form['end']
            pconn = psycopg2.connect(f"dbname={'crisprtehg38' if form['ga'] == 'hg38' else 'crisprtemm10'} user=postgres port=5432")
            pcursor = pconn.cursor()
            result = ''
            if types == "te":
                result = PGSQL.select_te_gtf_by_region(pcursor, "te", chrom, start, end)
                for i in range(len(result)):
                    for key in result[i].keys():
                        if type(result[i][key]) == memoryview:
                            result[i][key] = result[i][key].tobytes().decode()
                    result[i]['feature'] = 'te'
            elif types == "pcg":
                result = PGSQL.select_gtf_by_region(pcursor, "pcg", chrom, start, end)
                for i in range(len(result)):
                    for key in result[i].keys():
                        if type(result[i][key]) == memoryview:
                            result[i][key] = result[i][key].tobytes().decode()
            elif types == "all":
                result_gene = PGSQL.select_gtf_by_region(pcursor, "pcg", chrom, start, end)
                result_te = PGSQL.select_te_gtf_by_region(pcursor, "te", chrom, start, end) 
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
            return HttpResponse(Js.btoa(Pako.deflate(json.dumps(result))), content_type='application/octet-stream')
        elif form['key'] == 'getTtfByTE':
            te_class = form['te_class']
            pconn = psycopg2.connect(f"dbname={'crisprtehg38' if form['ga'] == 'hg38' else 'crisprtemm10'} user=postgres port=5432")
            pcursor = pconn.cursor()
            result = PGSQL.select_te_ttf_by_te(pcursor, te_info[ga]["te2int"][te_class])
            return HttpResponse(Js.btoa(Pako.deflate(json.dumps(result))), content_type='application/octet-stream')
        elif form['key'] == 'getDtfByTE':
            te_class = form['te_class']
            te_copy = form['te_copy']
            pconn = psycopg2.connect(f"dbname={'crisprtehg38' if form['ga'] == 'hg38' else 'crisprtemm10'} user=postgres port=5432")
            pcursor = pconn.cursor()
            result = PGSQL.select_te_dtf_by_te(pcursor, te_info[ga]["te2int"][te_class], te_copy)
            return HttpResponse(Js.btoa(Pako.deflate(json.dumps(result))), content_type='application/octet-stream')
