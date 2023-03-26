import gzip
import pybedtools
import pandas as pd

def reduce_func(x):
    priorityDict = {'promoter-TSS':4, 'exon':3, 'TE':2, 'intron':1, 'intergenic':0}
    x = list(x)
    tempvalue = []
    for i in x:
        tempkey = i.split('(')[0]
        tempvalue.append(priorityDict[tempkey])    
    return x[tempvalue.index(max(tempvalue))]

fullAnnotationBed = pybedtools.BedTool('./mm10_fullAnnotation.bed')
fullAnnotationBedsrt = fullAnnotationBed.sort()


def grna_annotate(gRNAbed):
    gRNAbedsrt = gRNAbed.sort()
    gRNAbedsrtdf = gRNAbedsrt.to_dataframe(names=[str(i) for i in range(1,10)])
    del(gRNAbedsrt)
    
    gRNAUniqBed = pybedtools.BedTool.from_dataframe(gRNAbedsrtdf.drop_duplicates(subset=['1','2','3','6']))
    gRNAannotate = fullAnnotationBedsrt.intersect(gRNAUniqBed, wb=True)
    gRNAannotateDf = gRNAannotate.to_dataframe(names=[str(i) for i in range(1,17)])
    gRNAannotateSimple = pd.DataFrame()
    gRNAannotateSimple['pos'] = gRNAannotateDf['7'].map(str) + '_' + gRNAannotateDf['8'].map(str) + '_' + gRNAannotateDf['9'].map(str) + '_' + gRNAannotateDf['12'].map(str)
    gRNAannotateSimple['gid'] = gRNAannotateDf['11']
    gRNAannotateSimple['anno_class'] = gRNAannotateDf['4']
    gRNAannotateSimple['pam'] = gRNAannotateDf['10']
    gRNAannotateSimple['upstream'] = gRNAannotateDf['13']
    gRNAannotateSimple['downstream'] = gRNAannotateDf['14']
    ## duplicate because of exon
    gRNAannotateSimpleUniq = gRNAannotateSimple.drop_duplicates()
    ### gRNA may annotate as several types
    #### priority promoter-TSS > exon > TE > intron > intergenic
    gRNAannotateSimpleUniqReduce = gRNAannotateSimpleUniq.groupby(["pos","gid","pam","upstream","downstream"]).agg({"anno_class": reduce_func})
    gRNAannotateSimpleUniqFeature = gRNAannotateSimpleUniqReduce['anno_class'].tolist()
    #  TE((GAATG)n;(GAATG)n;Satellite 
    gRNAaanotateTE_class = []
    gRNAaanotateTE_dup = []
    for i in gRNAannotateSimpleUniqFeature:
        if 'TE' in i:
            if '((' in i:
                tempAnnotate = i.split('(')
                tempAnnotateTE_class = '(' + tempAnnotate[2]
                tempAnnotateTE_class = tempAnnotateTE_class[:-1]
                tempAnnotateTE_dup = '(' + tempAnnotate[3].split(';')[0]
            else:
                tempAnnotate = i.split('(')[1]
                # print(tempAnnotate)
                tempAnnotateTE_class = tempAnnotate.split(';')[0]
                # print(tempAnnotateTE_class)
                tempAnnotateTE_dup = tempAnnotate.split(';')[1]
        else:
            tempAnnotateTE_class = ''
            tempAnnotateTE_dup = ''
        gRNAaanotateTE_class.append(tempAnnotateTE_class)
        gRNAaanotateTE_dup.append(tempAnnotateTE_dup)

    gRNAannotateSimpleUniqReduce['TE_class'] =  gRNAaanotateTE_class
    gRNAannotateSimpleUniqReduce['TE_dup'] = gRNAaanotateTE_dup
    gRNAannotateSimpleUniqReduce = gRNAannotateSimpleUniqReduce.reset_index()
    # print(gRNAannotateSimpleUniqReduce.shape)
    # print(gRNAbedsrtdf.drop_duplicates(subset=['1','2','3','6']).shape)
    return gRNAannotateSimpleUniqReduce