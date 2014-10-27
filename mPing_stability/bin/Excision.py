#!/opt/Python/2.7.3/bin/python
import sys
from collections import defaultdict
from numpy import *
import re
import os
import argparse
from Bio import SeqIO
import subprocess

def usage():
    test="name"
    message='''
python Excision.py --input ../input/10092013.mpings.gff
    '''
    print message

def fasta_id(fastafile):
    fastaid = defaultdict(str)
    for record in SeqIO.parse(fastafile,"fasta"):
        fastaid[record.id] = 1
    return fastaid

'''
Chr1    RelocaTE        mPing   10903901        10903903
'''
def readmPing(nonref):
    data = defaultdict(lambda: int)
    #data.update(readfile(ref, 0))
    #print len(data.keys())
    data.update(readfile(nonref, 1))
    #print len(data.keys())
    #data.update(readfile(shared, 2)) 
    #print len(data.keys())
    return data

def readfile(infile, flag):
    data = defaultdict(lambda: int)
    with open (infile, 'r') as filehd:
        for line in filehd:
            line = line.rstrip()
            if len(line) > 2: 
                unit = re.split(r'\t',line)
                mping = '%s:%s-%s' %(unit[0], unit[3], unit[4])
                data[mping] = flag
    return data

'''
Chr1    RIL10_0 RelocaTE        1132975 1132977 .       .       .       ID=mPing_1;Strain=RIL10_0;TSD=TAA;
'''
def readtable(infile):
    data = defaultdict(lambda: defaultdict(lambda: int()))
    p = re.compile(r'Strain=RIL(\d+)\_\d+;')
    with open (infile, 'r') as filehd:
        for line in filehd:
            line = line.rstrip()
            if len(line) > 2: 
                unit = re.split(r'\t',line)
                mping = '%s:%s-%s' %(unit[0], unit[3], unit[4])
                m = p.search(unit[8])
                strain = m.groups(0)[0] if m else 'NA'
                #print strain, mping, line
                data[strain][mping] = 1
    return data


'''
Chr1    RIL10_0 RelocaTE        1132975 1132977 .       .       .       ID=mPing_1;Strain=RIL10_0;TSD=TAA;
'''
def mping_frequency(infile):
    data = defaultdict(lambda: defaultdict(lambda: int))
    rils = defaultdict(int)
    inf  = defaultdict()
    p = re.compile(r'Strain=RIL(\d+)\_\d+;')
    with open (infile, 'r') as filehd:
        for line in filehd:
            line = line.rstrip()
            if len(line) > 2: 
                unit = re.split(r'\t',line)
                mping = '%s:%s-%s' %(unit[0], unit[3], unit[4])
                m = p.search(unit[8])
                strain = m.groups(0)[0] if m else 'NA'
                data[mping][strain] = 1
                rils[strain] =1
                inf[mping] = [unit[0], unit[3], unit[4]]
    total = len(rils.keys())
    #print total
    mping_frq = defaultdict(lambda: int)
    for m in data.keys():
        count = len(data[m].keys())
        frq = float(count)/total
        #print '%s\t%s\t%s\t%s\t%s\t%s\t%s' %(inf[m][0], inf[m][1], inf[m][2], m, '+', str(count), str(frq))
        mping_frq[m] = frq
    return mping_frq


#Convert BIN MAP
'''
""      "GN1"   "GN10"  "GN100" "GN101" "GN102" "GN103" "GN104" "GN105" "GN106" "GN107" "GN108
"0100222046"    1       1       0       0       0       0       1       1       1       1
"0100500860"    1       1       0       0       0       0       1       1       1       1
'''

def convert_MAP(infile):
    rils = []
    data = defaultdict(lambda : defaultdict(lambda: defaultdict(lambda: int)))
    with open (infile, 'r') as filehd:
        for line in filehd:
            line = line.rstrip()
            line = line.replace('"', '')
            if line[0:1].isdigit():
                unit = re.split(r'\t',line)
                #print '%s\t%s\t%s' %(chrs, int(unit[0][2:]), str(chr_start[chrs]))
                chrs = 'Chr%s' %(str(int(unit[0][0:2])))
                pos = int(unit[0][2:])
                for i in range(1,len(unit)):
                    #print i, rils[i], chrs, pos
                    ril = rils[i]
                    data[ril][chrs][pos] = unit[i]
            else:
                unit = re.split(r'\t',line)
                unit[0] = 'Position'
                #print unit[0], unit[1], unit[2]
                rils.extend(unit)
    return data

def binarySearch(data, val):
    highIndex = len(data)-1
    lowIndex = 0
    while highIndex > lowIndex:
            index = (highIndex + lowIndex) / 2
            sub = int(data[index])
            #print highIndex, index, lowIndex, sub, val
            if data[lowIndex] == val:
                    return [lowIndex, lowIndex]
            elif sub == val:
                    return [index, index]
            elif data[highIndex] == val:
                    return [highIndex, highIndex]
            elif sub > val:
                    if highIndex == index:
                            return sorted([highIndex, lowIndex])
                    highIndex = index
            else:
                    if lowIndex == index:
                            return sorted([highIndex, lowIndex])
                    lowIndex = index
    return sorted([highIndex, lowIndex])
 

def findbin(start, binmap, ril, chrs, transition):
    #for i in sorted(binmap[ril][chrs].keys(), key=int):
    #    print i
    array = []
    array.extend(sorted(binmap[ril][chrs].keys(), key=int))
    #for i in array:
        #print i
    if int(start) > int(array[-1]):
        return 0
    index = binarySearch(array, int(start))
    #print 'BINsearch', len(array), start, index[0], index[1], array[index[0]], array[index[1]],binmap[ril][chrs][array[index[1]]]
    if transition[ril][chrs].has_key(array[index[1]]):
        return 0
    return array[index[1]]

def genotyping(ril, mping, binmap, transition):
    ril = 'GN%s' %(ril)
    p = re.compile(r'(\w+):(\d+)\-(\d+)')
    m = p.search(mping)
    chrs = ''
    start = 0
    end   = 0
    if m:
        chrs  = m.groups(0)[0]
        start = m.groups(0)[1]
        end   = m.groups(0)[2]
    #print ril, chrs, start, len(binmap[ril][chrs].keys())
    pos = findbin(start, binmap, ril, chrs, transition)
    #print ril,chrs,start,pos
    genotype = binmap[ril][chrs][pos] if pos > 0 else '3'
    return genotype

def validmap(binmap):
    last = 0
    data = defaultdict(lambda : defaultdict(lambda: defaultdict(lambda: int))) 
    for ril in sorted(binmap.keys()):
        for chrs in sorted(binmap[ril].keys()):
            for pos in sorted(binmap[ril][chrs].keys(), key=int):
                if binmap[ril][chrs][pos] != binmap[ril][chrs][last]:
                    data[ril][chrs][last] = 1
                    data[ril][chrs][pos] = 1
                    last = pos
                    #print ril, chrs, pos
    return data

def bamcheck(bam, mping):
    reg = re.compile(r'(Chr\d+):(\d+)\-(\d+)')
    match = reg.search(mping)
    start = int(match.groups(0)[1])
    end   = int(match.groups(0)[2])
    cmd = 'samtools view %s %s' %(bam, mping)
    ofile = open('bamcheck.txt', 'a') 
    print >> ofile, bam, mping
    print >> ofile, cmd
    out = subprocess.check_output(cmd, shell=True)
    lines = re.split(r'\n', out)
    covered = 0
    clipped  = 0
    if len(lines) < 2:
        return 2
    pattern = re.compile(r'([0-9]+)([A-Z]+)')

    for line in lines:
        print >> ofile, line
        unit = re.split(r'\t', line)
        if len(unit) < 2:
            continue
        matches = []
        for (base, match) in re.findall(pattern, unit[5]):
            if match == 'I' or match == 'D':
                continue
            #print >> ofile, base, match
            matches.append([base, match])
            #print >> ofile, matches[0]
        #print >> ofile, len(matches)
        if len(matches) == 1 and matches[0][1] == 'M': # perfect match
            if int(unit[3]) < start - 10 and int(unit[3]) + int(matches[0][0]) > end + 10: # read cover mping insertion site
                covered += 1
        elif len(matches) == 2: # may have soft clip
            if matches[0][1] == 'S' and int(matches[0][0]) > 10 and matches[1][1] == 'M' and int(matches[1][0]) > 10: # clip at start
                if int(unit[3]) > start - 5 and int(unit[3]) < end + 5: # read at mping insertion site, so the clip should be due to mping insertion
                    covered += 1
                    clipped += 1
                elif int(unit[3]) <= start - 5 and int(unit[3]) + int(matches[1][0]) >= end + 5: # read cover mping insertion site
                    covered += 1
            elif matches[1][1] == 'S' and int(matches[1][0]) > 10 and matches[0][1] == 'M' and int(matches[0][0]) > 10: # clip at end
                if int(unit[3]) + int(matches[0][0]) > start - 5 and int(unit[3]) + int(matches[0][0]) < end + 5: # read at mping insertion site, so the clip should be due to mping insertion
                    covered += 1
                    clipped += 1
                elif int(unit[3]) <= start - 5 and int(unit[3]) + int(matches[0][0]) >= end + 5: # read cover mping insertion site
                    covered += 1
        elif len(matches) == 3 and matches[0][1] == 'S' and matches[1][1] == 'M' and matches[2][1] == 'S': # may have soft clip, but the other end of reads have clip too
            if int(unit[3]) > start - 5 and int(unit[3]) < end + 5 and int(matches[0][0]) > 10: # read at mping insertion site, so the clip should be on the left if due to mping
                covered += 1
                clipped += 1
            elif int(unit[3]) + int(matches[1][0]) > start - 5 and int(unit[3]) + int(matches[1][0]) < end + 5: # read start before mping insertion site, but clipped at mping
                covered += 1
                clipped += 1
            elif int(unit[3]) <= start - 5 and int(unit[3]) + int(matches[1][0]) >= end + 5: # read cover the mping insertion site
                covered += 1
        print >> ofile, covered, clipped
    print >> ofile, covered, clipped
    rate = float(clipped/covered) if covered > 0 else 0
    if rate > 0.2:
        return 1
    else:
        return 0
 
def excision(mPing_ancestor, mPing_rils, mPing_frq):
    mPing_excision = defaultdict(lambda : defaultdict(lambda : int()))
    binmap = convert_MAP('../input/MPR.geno.bin.uniq')
    transition = validmap(binmap)
    for mping in mPing_ancestor.keys():
        print mping
        if not mPing_frq.has_key(mping) or mPing_frq[mping] < 0.05:
            continue
            print 'not in rils'
        print 'in rils'
        for ril in mPing_rils.keys():
            #print 'TT',ril,len(binmap[ril].keys())
            if not binmap.has_key('GN%s' %(ril)):
                continue
            genotype = genotyping(ril, mping, binmap, transition)
            bam = '/rhome/cjinfeng/BigData/00.RD/RILs/QTL_pipe/input/fastq/RILs_ALL_bam/GN%s.bam' %(ril)
            flag = 2
            #print mping, ril, genotype
            #if mping == 'Chr1:1715117-1715119' and int(ril) == 10:
                #print 'EX', mping, ril, genotype, mPing_rils[ril][mping], mPing_ancestor[mping]
            if (int(genotype) == 0 and int(mPing_ancestor[mping]) != 1): # ril has genotype of reference and mping has genotype of reference or both
                #print 'S1', mPing_rils[ril].has_key(mping)
                #if os.path.isfile(bam):
                #    flag = bamcheck(bam, mping)
                if not mPing_rils[ril].has_key(mping):
                    if os.path.isfile(bam):
                        flag = bamcheck(bam, mping)
                        if flag == 0: ## bam check showed no mping insertion in this ril
                            mPing_excision[mping][ril] = 1
                    #print mPing_excision[mping]
            elif (int(genotype) == 1 and int(mPing_ancestor[mping]) != 0): # ril has genotype of nonref and mping has genotype of nonref or both
                #print 'S2', mPing_rils[ril].has_key(mping)
                #if os.path.isfile(bam):
                #    flag = bamcheck(bam, mping)
                if not mPing_rils[ril].has_key(mping):
                    if os.path.isfile(bam):
                        flag = bamcheck(bam, mping)
                        if flag == 0: ## bam check showed no mping insertion in this ril
                            mPing_excision[mping][ril] = 1
                    #print mPing_excision[mping]
    return mPing_excision

def preplot(frq, ancestor, excision):
    for m in excision.keys():
        print '>', m, len(excision[m].keys())
        for i in sorted(excision[m].keys()):
            print i

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input')
    parser.add_argument('-o', '--output')
    parser.add_argument('-v', dest='verbose', action='store_true')
    args = parser.parse_args()
    try:
        len(args.input) > 0
    except:
        usage()
        sys.exit(2)

    mPing_ancestor = readmPing('../input/HEG4.mping.non-ref.gff')
    mPing_rils     = readtable(args.input)
    mPing_frq      = mping_frequency(args.input)
    mPing_excision = excision(mPing_ancestor, mPing_rils, mPing_frq)
    preplot(mPing_frq, mPing_ancestor, mPing_excision)

if __name__ == '__main__':
    main()

