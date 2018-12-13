from DeepSeqG9 import *
from ProcessBenign import *
import pysam
import math
import numpy as np
#PDB2Uniprot = "../dat/test.map"
#PDB2Uniprot = "../dat/pdbsws_res.txt"
#ins.Reformat(PDB2Uniprot, "../dat/PDB2SWISS.tsv")
#PDB2Uniprot = "../dat/PDB2SWISS.sort.tsv"
#ins.ProcessPDB2Uniprot(PDB2Uniprot)

ExonTab = "../dat/gencodeV19.selected.Transcripts.bed"
hg19 = "/share/data/resources/hg19/references/hg19.fasta"
allmisFil = "../dat/revel_all_chromosomes.txt.gz"
#allmisFil = "../dat/SYN/ALL.MIS.sort.txt.gz"
#allsynFil = "../dat/SYN/ALL.SYN.sort.txt.gz"

#allmisFil = "../dat/SYN/tmp/test.mis.sort.txt.gz"
#allsynFil = "../dat/SYN/tmp/test.sort.txt.gz"

gnomAD_Cov_Fil = "/share/data/resources/ANNOVAR_DATA/gnomAD/WES/Coverage/gnomad.exomes.coverage.txt.gz"

N = 1e7

def LoadMutable():
    #reader = csv.reader(open("../dat/rate.txt", 'rb'), delimiter="\t")
    fin = open("../dat/rate.txt", 'rb')
    res = {}
    for l in fin:
        ref, alt, rate = l.strip().split()
        res[ref+">"+alt] = float(rate)
    return res

def GenerateMisRefFil():
    Dir = "/home/jw3514/DeepseqG9/Jiayao/dat/DBNSFP/"
    Chrs = map(str, range(1,23))
    Chrs.extend(["X", "Y"])
    res = {}
    for CHR in Chrs:
        res[CHR] = None

def LoadTrans():
    fin = open(ExonTab, 'rb')
    res = {}
    for l in fin:
        row = l.strip().split("\t")
        gene = row[6]
        trans = row[3].split("_")[0]
        exonid = row[3].split("_")[2]
        Chr = row[0].strip("chr")
        start, end = int(row[1]), int(row[2])
        if gene not in res:
            res[gene] = []
        res[gene].append([trans, exonid, Chr, start, end])
    return res
    
def TestGene(Gene, GeneExons, HG19, ALLMIS, mutable):
    Rate = 0
    for trans, exonid, Chr, start, end in GeneExons[Gene]:
        for v in ALLMIS.fetch(Chr, start-1, end):
            v = v.split()
            Chr, Pos, Ref, Alt = v[0], int(v[1]), v[2], v[3]
            Trinuc = HG19.fetch(Chr, Pos-2, Pos+1)
            TrinucCG = Trinuc[0] + Alt + Trinuc[2]
            rate = mutable[Trinuc+">"+TrinucCG]
            Rate += rate
    print Rate

    pass

def AggregateMutationRate(Chr, start, end, HG19, ALLMIS, mutable):
    Rate = 0
    for v in ALLMIS.fetch(Chr, start-1, end):
        v = v.split()
        Chr, Pos, Ref, Alt = v[0], int(v[1]), v[2], v[3]
        Trinuc = HG19.fetch(Chr, Pos-2, Pos+1)
        TrinucCG = Trinuc[0] + Alt + Trinuc[2]
        rate = mutable[Trinuc+">"+TrinucCG]
        Rate += rate
    return Rate


def AggregateMutationRate2(Chr, start, end, TabixFil):
    Rate = 0
    for v in TabixFil.fetch(Chr, start-1, end):
        v = v.split()
        #Chr, Pos, Rate = v[0], int(v[1]), v[2], v[3]
        rate = float(v[12])
        Rate += rate
    return Rate

def ExpectedCount(Beta0, Beta1, Rate):
    return Beta0 + Beta1 * Rate * N

def ExpectedCountCovCorrection(Cov):
    cov = getDepth(Chr, start, end, Cov)
    cor = CorrectCount(ExpC, cov)

def GetDepth(Chr, start, end, Cov):
    covs = Cov.fetch(Chr, start, end)
    Medians = []
    for l in covs:
        row = l.split()
        median = float(row[3])
        Medians.append(median)
    if Medians == []:
        return 50
    else:
        return np.mean(Medians)

def CorrectCount(ExpC, Cov):
    if Cov >= 50:
        return ExpC
    elif Cov >1:
        return ExpC * (0.089 + 0.217 * math.log(Cov))
    else:
        return 0.089 * ExpC

def Zscore(obs, exp):
    chi2 = (obs-exp)**2/exp
    if obs > exp:
        return -1 * math.sqrt(chi2)
    else:
        return math.sqrt(chi2)


def ComputeRatePerExon(Tab, HG19, ALLMIS, mutable):
    fin = open(Tab, 'rb')
    fout = csv.writer(open("ExonMuteRate.tsv", 'wb'), delimiter="\t")
    res = {}
    for l in fin:
        row = l.strip().split("\t")
        Chr = row[0].strip("chr")
        start, end = int(row[1]), int(row[2])
        Rate =  AggregateMutationRate(Chr, start, end, HG19, ALLMIS, mutable)
        row.append(Rate)
        fout.writerow(row)

def ComputeObsNumPerExon(Tab, gnomAD, DoC):
    reader = csv.reader(open(Tab, 'rb'), delimiter="\t")
    writer = csv.writer(open("Regioins.addObsDepth.tsv", 'wb'), delimiter="\t")
    for row in reader:
        Chr, Start, End = row[0].strip("chr"), int(row[1]), int(row[2])
        syn_obs, mis_obs = GetObsSynMis(gnomAD, Chr, Start, End)
        Depth = GetDepth(Chr, Start, End, DoC)
        row.extend([Depth, syn_obs, mis_obs])
        writer.writerow(row)

def GetObsSynMis(gnomAD, Chr, Start, End):
    syn, mis = 0, 0
    for v in gnomAD.fetch(Chr, Start-1, End):
        llist = v.split("\t")
        INFO_string = llist[7]
        info_dict = getINFO(INFO_string)
        CSQ = info_dict["CSQ"][0].split("|")
        CSQ = dict(zip(CSQ_FMT, CSQ))
        if CSQ["Consequence"] == "synonymous_variant":
            syn += 1
        elif CSQ["Consequence"] == "missense_variant":
            mis += 1
    return syn, mis

# Things to compute: 
# SYN_Rate MIS_Rate
def ComputeRatePerExon2(Table, SYN_ALL, MIS_ALL):
    fin = open(Table, 'rb')
    fout = csv.writer(open("ExonMuteRate.tsv", 'wb'), delimiter="\t")
    res = {}
    for l in fin:
        row = l.strip().split("\t")
        Chr = row[0].strip("chr")
        start, end = int(row[1]), int(row[2])
        SYNRate =  AggregateMutationRate2(Chr, start, end, SYN_ALL)
        MISRate =  AggregateMutationRate2(Chr, start, end, MIS_ALL)
        row.extend([SYNRate, MISRate])
        fout.writerow(row)

def ComputeExpPerExon(HG19, ALLMIS, mutable):
    fin = open(ExonTab, 'rb')
    fout = csv.writer(open("ExonMuteRate.tsv", 'wb'), delimiter="\t")
    res = {}
    for l in fin:
        row = l.strip().split("\t")
        Chr = row[0].strip("chr")
        start, end = int(row[1]), int(row[2])
        Rate =  AggregateMutationRate(Chr, start, end, HG19, ALLMIS, mutable)
        row.append(Rate)
        fout.writerow(row)

def ComputeOBSPerExon(HG19, ALLMIS, mutable):
    fin = open(ExonTab, 'rb')
    fout = csv.writer(open("ExonMuteRate.tsv", 'wb'), delimiter="\t")
    res = {}
    for l in fin:
        row = l.strip().split("\t")
        Chr = row[0].strip("chr")
        start, end = int(row[1]), int(row[2])
        Rate =  AggregateMutationRate(Chr, start, end, HG19, ALLMIS, mutable)
        row.append(Rate)
        fout.writerow(row)

def main():
    mutable = LoadMutable()
    HG19 = pysam.FastaFile(hg19)
    GeneExons = LoadTrans()
    ALLMIS = pysam.TabixFile(allmisFil)
    #ALLSYN = pysam.TabixFile(allsynFil)
    Cov = pysam.TabixFile(gnomAD_Cov_Fil)
    gnomAD = pysam.TabixFile("gnomad.Exon.Rare.HighConf.vcf.gz")
    #print HG19.fetch('1', 67000041, 67000051)
    #print ALLMIS.fetch('1', 67000041, 67000051)
    #TestGene("KLHL17", GeneExons, HG19, ALLMIS, mutable)
    Tab = "Regioins.addObsDepth.tsv"
    ComputeRatePerExon(Tab, HG19, ALLMIS, mutable) 
    #ComputeRatePerExon2(Table, ALLSYN, ALLMIS) 
    #ComputeObsNumPerExon(ExonTab, gnomAD, Cov)


main()
