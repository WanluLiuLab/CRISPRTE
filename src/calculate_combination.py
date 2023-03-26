import numpy as np
import pickle
import tqdm
from itertools import groupby
from queue import PriorityQueue
from copy import deepcopy
from collections import Counter 
import pandas as pd
import ssw 

from src.dbutils import *

class Combination:
    def __init__(self, teid, tecpyn, cursor = None, consensus = None):
        self.gids = []
        self.te_targets = set()
        self.offtargets = {
            'promoter-TSS': 0,
            'exon': 0,
            'intron': 0,
            'intergenic': 0
        }
        self.teid = teid 
        self.tecpyn = tecpyn
        self.gseq = {}
        self.consensus = consensus

    def append(self, annotation):
        self.gids.append(annotation[0])
        self.gseq[annotation[0]] = SQLite.select_gseq_table2_gid(cursor, annotation[0])
        self.te_targets = self.te_targets.union(annotation[1]["TE"])
        for i in annotation[1].keys():
            if i != "TE":
                self.offtargets[i] += annotation[1][i]

    def combination_score(self):
        on_target = len(list(filter(lambda x:x[0] == self.teid, self.te_targets)))
        on_target_perc = on_target / self.tecpyn * 100
        off_target_te = len(list(filter(lambda x:x[0] != self.teid, self.te_targets)))
        off_target_penalty = .04 * self.offtargets["promoter-TSS"] + .03 * self.offtargets["exon"] + .02 * self.offtargets["intron"] + .01 * self.offtargets["intergenic"]
        return on_target_perc - off_target_te * 0.001 - off_target_penalty * 0.0001

    def tostring(self):
        on_target = len(list(filter(lambda x:x[0] == self.teid, self.te_targets)))
        on_target_perc = on_target / self.tecpyn * 100
        return "Coverage: {}({}) of {}, score {}".format(on_target, on_target_perc, self.tecpyn, self.combination_score())

    def totuple(self):
        on_target = len(list(filter(lambda x:x[0] == self.teid, self.te_targets)))
        on_target_perc = on_target / self.tecpyn * 100
        off_target_te = list(map(lambda y: y[0], filter(lambda x:x[0] != self.teid, self.te_targets)))
        off_target_te = dict(Counter(list(map(lambda x:int2te[x], off_target_te))))
        return (self.gids, on_target, on_target_perc, self.tecpyn, self.combination_score(), off_target_te, self.offtargets)

    def gseq(self, cursor):
        return dict(zip(self.gids, list(map(lambda x: SQLite.select_gseq_table2_gid(cursor, x), self.gids))))

    def align2consensus(self):
        raise NotImplementedError("align2consensus not implemented")


##################  
# SQLite version #
##################

parser = argparse.ArgumentParser()
parser.add_argument('-t', '--te', type=str,
                    help='TE name')
args = parser.parse_args()

TE = args.te

with open("/home/zje610/workspace/labW/Snowxue/refData/hg38_TE/te2int-hg38.pkl","rb") as f:
    te2int = pickle.load(f)
    int2te = {i:j for j, i in te2int.items()}

with open("/home/zje610/workspace/labW/Snowxue/refData/hg38_TE/TE-copy-hg38.pkl","rb") as f:
    te_info = pickle.load(f)

with open("/home/zje610/workspace/labW/Snowxue/refData/mm10_TE/te2int-mm10.pkl", "rb") as f: 
    te2int = pickle.load(f)
    int2te = {i:j for j, i in te2int.items()}

with open("/home/zje610/workspace/labW/Snowxue/refData/mm10_TE/TE-copy-mm10.pkl","rb") as f:
    te_info = pickle.load(f)

    
result = sql.select_gid_seqinfo_table1_teclass(cursor, TE)
gids = np.unique(list(map(lambda x:x[0], result)))


def annotation_compile(annotations):
    result = {"TE":[], "promoter-TSS":0, "exon":0, "intron":0, "intergenic":0}
    for annotation in annotations:
        if annotation[0].split("(")[0] == "TE":
            try:
                if len(annotation[2].split("_dup")) <2:
                    result["TE"].append((te2int[annotation[2]], 0))
                result["TE"].append((te2int[annotation[2].split("_dup")[0]], int(annotation[2].split("_dup")[1])))
            except:
                pass
        else:
            result[annotation[0].split("(")[0]] += 1
    return result

result_annotation = {}

for i in tqdm.trange(len(gids)):
    tmp =  SQLite.select_anno_table1_gid(cursor,gids[i])
    result_annotation[gids[i]] = annotation_compile(tmp)

normalize = lambda v:v / np.sqrt(np.sum(v**2)) * 100

def distribution(lst, n):
    result = {}
    for k,g in groupby(sorted(lst),key=lambda x:x//n):
        #  print('{}-{}:{}'.format(k*n,(k+1)*n-1,len(list(g))))
        l = len(list(g))
        result['{}-{}'.format(k*n,(k+1)*n)] = l
    return result


def gScore(annotation, teid, tecpyn):
    on_target = len(np.unique(list(filter(lambda x:x[0] == teid, annotation["TE"]))))
    on_target_perc = on_target / tecpyn * 1000
    off_target_te = len(list(filter(lambda x:x[0] != teid, annotation["TE"])))
    off_target_penalty = .04 * annotation["promoter-TSS"] + .03 * annotation["exon"] + .02 * annotation["intron"] + .01 * annotation["intergenic"]
    return on_target_perc - off_target_te * 0.001 - off_target_penalty * 0.0001

def greedyGidCombination(layer, teid, branch = 10, dropouts = 5, nrec = 3, _currec = None):
    
    if not nrec:
        return list(map(lambda x:x[1], layer))
    assert(len(layer) > 0)
    next_layer = []
    if not _currec:
        for pq, comb in layer:
            cpq = PriorityQueue()
            cpq.queue = deepcopy(pq.queue)
            inc = []
            ccombs = []
            for i in range(branch):
                ccomb = deepcopy(comb)
                if not len(cpq.queue):
                    return 
                gid = cpq.get()[1]
                ccomb.append((gid,result_annotation[gid]))
                inc.append(ccomb.combination_score() - comb.combination_score())
                ccombs.append(ccomb)
            idxs = np.argsort(-np.array(inc))[:branch - dropouts]
            for idx in idxs:
                next_layer.append((cpq, ccombs[idx]))
        return greedyGidCombination(next_layer, teid, branch, dropouts, nrec - 1, 1)

    for pq, comb in layer:
        cpq = PriorityQueue()
        for item in pq.queue:
            annotation = set(list(map(lambda y: y[1], filter(lambda x: x[0] == teid, result_annotation[item[1]]['TE']))))
            difference = len(annotation.difference(comb.te_targets))
            cpq.put( ((-difference/100 + item[0]), item[1]))
        # cpq.queue = deepcopy(pq.queue)
        inc = []
        ccombs = []
        for i in range(branch):
            ccomb = deepcopy(comb)
            gid = cpq.get()[1]
            ccomb.append((gid, result_annotation[gid]))
            inc.append(ccomb.combination_score() - comb.combination_score())
            ccombs.append(ccomb)
        idxs = np.argsort(-np.array(inc))[:branch - dropouts]
        for idx in idxs:
            next_layer.append((cpq, ccombs[idx]))
    return greedyGidCombination(next_layer, teid, branch, dropouts, nrec - 1, 1)

def calculate_combination(result_annotation, TE):
    score = list(map(lambda x:gScore(x, te2int[TE], te_info[TE]), result_annotation.values()))
    # norm_score = normalize(np.array(score))
    pq = PriorityQueue()
    for sc, gid in zip(score, result_annotation.keys()):
        pq.put((-sc, gid))
    result = greedyGidCombination( [(pq,  Combination( te2int[TE], te_info[TE]))], te2int[TE], 10, 5, 3)
    if not result:
        return None
    return pd.DataFrame(list(map(lambda x:x.totuple(), result)))



######################
# END SQLite version #
######################


result_seq = list(map(lambda x:SQLite.select_gseq_table2_gid(cursor,x[0]),result))
gRNA = list(zip(result_seq, result, result_annotation))
gRNA = sorted(gRNA, key=lambda x: -len(x[2]))

print("Found {} on {}".format(len(gRNA), TE))

def unique_by_first(lst):
   s = set()
   out = []
   for i in lst:
      if i[0] not in s:
         out.append(i)
         s.add(i[0])
   return out

gRNA = unique_by_first(gRNA)

popped = []
for i in range(len(gRNA)):
    for j in range(i):
        if set(gRNA[i][2]).issubset(set(gRNA[j][2])):
            popped.append(i)
            continue

gRNA_new = []
for i,j in enumerate(gRNA):
    if i not in popped:
        j = list(j)
        n = []
        for x in j[2]:
            try:
                n.append(int(x.split("dup")[-1]))
            except:
                n.append(0)
        j[2] = n
        gRNA_new.append(j)
        
print("Total {} unique gRNA on {}".format(len(gRNA_new), TE))


LTR7Y_combination_gids = [i[0][1][0] for i in LTR7Y_combination_gRNA]

s = set(gRNA_new[0][2])
index = set([0])
for i in range(20):
    l = []
    for i in range(len(gRNA_new)):
        if i not in index:
            l.append(len(s.union(gRNA_new[i][2])))
        else:
            l.append(0)
    ind = np.argmax(l)
    index.add(ind)
    s = s.union(gRNA_new[ind][2])
    print(len(s) / 257,end="\t")
    print(index)