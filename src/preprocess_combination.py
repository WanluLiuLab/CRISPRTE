import numpy as np
import pickle
from dbutils import *
import tqdm
import argparse 
from dbutils import *

##################
# SQLite version #
##################

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-t', '--te', type=str,
                        help='TE name')
    parser.add_argument('-s', '--species', options=["hg38", "mm10"], type=str)
    args = parser.parse_args()


    TE = args.te

    if args.species == 'hg38':
        with open("/home/zje610/workspace/labW/Snowxue/refData/hg38_TE/TE-copy.pkl","rb") as f:
            te_info = pickle.load(f)
            te2int = dict(zip(list(map(lambda x:"-".join(x.split("_")), te_info.keys())),range(len(te_info))))
        db_path='./Homo_sapiens.GRCh38.97.dna.primary_assembly.gRNA.db'
    else:
        with open("/home/zje610/workspace/labW/Snowxue/refData/mm10_TE/te2intmm10.pkl", "rb") as f:
            te2int = pickle.load(f)
        db_path='./Mus_musculus.GRCm38.97.dna.primary_assembly.gRNA.db'

    
    sql = SQLite()
    conn = sql.connect(db_path)
    cursor = conn.cursor()
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

    with open("./pkl/{}/{}_gid_allinfo.pkl".format(args.species, args.te), "wb+") as f:
        pickle.dump(result_annotation,f)

