# -*- coding: utf-8 -*-
#!/usr/bin/env python
# Author: Ziwei Xue
#
########## pasrse Data ##############
def parseData(v, k):
    result = {}
    mmsameTEcount = [ v[k]['sameTE'] if 'sameTE' in v[k].keys() else 0][0]
    mmdiffTEcount = [ v[k]['diffTE'] if 'sameTE' in v[k].keys() else 0][0]
    mmintergeniccount = [ v[k]['intergenic'] if 'intergenic' in v[k].keys() else 0][0]
    mmgeniccount = [ v[k]['exon'] if 'exon' in v[k].keys() else 0][0] + [ v[k]['intron'] if 'intron' in v[k].keys() else 0][0]  + [ v[k]['promoter-TSS'] if 'promoter-TSS' in v[k].keys() else 0][0]
    result['sameTEcount'] = mmsameTEcount
    result['diffTEcount'] = mmdiffTEcount
    result['intergeniccount'] = mmintergeniccount
    result['geniccount'] = mmgeniccount
    return result

def formatData(data):
    formatdata = {}
    for k,v in data.items():
        tempDict = {}
        tempDict['gseq'] = v['gseq']
        tempDict['pos'] = '_'.join(v['pos'].split('_')[:3])
        tempDict['strand'] = v['pos'].split('_')[3]
        tempDict['OT_score'] = v['score']
        tempDict['initial_score'] = v['gscore']
        tempMM0,tempMM1,tempMM2,tempMM3 = v['mm0'], v['mm1'], v['mm2'], v['mm3']
        tempMM0count = [ sum(tempMM0.values()) -1 if (sum(tempMM0.values()) -1)  > 0 else 0 ][0]
        tempMM1count = sum(tempMM1.values())
        tempMM2count = sum(tempMM2.values())
        tempMM3count = sum(tempMM3.values())
        tempDict['mm0count'],tempDict['mm1count'], tempDict['mm2count'],tempDict['mm3count'],= tempMM0count, tempMM1count, tempMM2count, tempMM3count
        tempDict['OT_count'] =  tempMM0count + tempMM1count + tempMM2count+ tempMM3count
        tempDict['sameTEcount'] = parseData(v, 'mm0')['sameTEcount'] + parseData(v, 'mm1')['sameTEcount'] + parseData(v, 'mm2')['sameTEcount'] + parseData(v, 'mm3')['sameTEcount']
        tempDict['diffTEcount'] = parseData(v, 'mm0')['diffTEcount'] + parseData(v, 'mm1')['diffTEcount'] + parseData(v, 'mm2')['diffTEcount'] + parseData(v, 'mm3')['diffTEcount']
        tempDict['geniccount'] = parseData(v, 'mm0')['geniccount'] + parseData(v, 'mm1')['geniccount'] + parseData(v, 'mm2')['geniccount'] + parseData(v, 'mm3')['geniccount']
        tempDict['intergeniccount'] = parseData(v, 'mm0')['intergeniccount'] + parseData(v, 'mm1')['intergeniccount'] + parseData(v, 'mm2')['intergeniccount'] + parseData(v, 'mm3')['intergeniccount']
        formatdata[k] = tempDict
    return formatdata

s2n = {'A':0, 'C':1, 'G':2, 'T':3}
n2s = ['A','C','G','T']

def SymbolToNumber(symbol: str) -> int:
	"""
	@brief: should be ordered by lexicographical order  {'A':0, 'C':1, 'G':2, 'T':3}
	"""
	return s2n[symbol]

def NumberToSymbol(num: int) -> str:
	"""
	@brief: reverse function of SymbolToNumber
	"""
	return n2s[num]

def PatternToNumber(pattern: str) -> int:
	if (len(pattern) == 0):
		return 0
	else:
		return 4 * PatternToNumber(pattern[:-1]) + SymbolToNumber(pattern[-1])

def NumberToPattern(index: int, k: int) -> str:
	"""
	k should be the length of output string
	"""
	if k == 1:
		return NumberToSymbol(index)
	prefix_index = index // 4
	remainder = index % 4
	symbol = NumberToSymbol(remainder)
	prefix_pattern = NumberToPattern(prefix_index, k - 1)
	return prefix_pattern + symbol