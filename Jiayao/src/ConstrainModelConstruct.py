from DeepSeqG9 import *
import pysam
#PDB2Uniprot = "../dat/test.map"
#PDB2Uniprot = "../dat/pdbsws_res.txt"
#ins.Reformat(PDB2Uniprot, "../dat/PDB2SWISS.tsv")
#PDB2Uniprot = "../dat/PDB2SWISS.sort.tsv"
#ins.ProcessPDB2Uniprot(PDB2Uniprot)

ExonTab = "../dat/gencodeV19.selected.Transcripts.bed"
hg19 = "/share/data/resources/hg19/references/hg19.fasta"
#allmisFil = "../dat/ALL.Mis.txt.gz"
allmisFil = "../dat/revel_all_chromosomes.txt.gz"
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

def ComputeRatePerExon(HG19, ALLMIS, mutable):
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
    #GeneExons = LoadTrans()
    ALLMIS = pysam.TabixFile(allmisFil)
    #print HG19.fetch('1', 67000041, 67000051)
    #print ALLMIS.fetch('1', 67000041, 67000051)
    #TestGene("KLHL17", GeneExons, HG19, ALLMIS, mutable)
    ComputeRatePerExon(HG19, ALLMIS, mutable) 


main()
