#from Mutations import *
import gzip

gnomad_file = "/Users/jiayao/Work/DeepSeq/dat/gnomAD/test.vcf.gz"
#gnomad_file = "/Users/jiayao/Work/DeepSeq/dat/gnomAD/haha.vcf.gz"

head = ['genename', 'assembly', 'chr', 'pos', 'ref', 'alt', 'ref_aa', 'alt_aa', 'swissprot_pos', 'ensp', 'swissprot',
        'dbsnp', 'source', 'source_id', 'target', 'descr', 'var_id', 'swissprot_var_id', 'af']
print zip(range(19), head)

CSQ_FMT = 'Allele|Consequence|IMPACT|SYMBOL|Gene|Feature_type|Feature|BIOTYPE|EXON|INTRON|HGVSc|HGVSp|cDNA_position|CDS_position|Protein_position|Amino_acids|Codons|Existing_variation|ALLELE_NUM|DISTANCE|STRAND|FLAGS|VARIANT_CLASS|MINIMISED|SYMBOL_SOURCE|HGNC_ID|CANONICAL|TSL|APPRIS|CCDS|ENSP|SWISSPROT|TREMBL|UNIPARC|GENE_PHENO|SIFT|PolyPhen|DOMAINS|HGVS_OFFSET|GMAF|AFR_MAF|AMR_MAF|EAS_MAF|EUR_MAF|SAS_MAF|AA_MAF|EA_MAF|ExAC_MAF|ExAC_Adj_MAF|ExAC_AFR_MAF|ExAC_AMR_MAF|ExAC_EAS_MAF|ExAC_FIN_MAF|ExAC_NFE_MAF|ExAC_OTH_MAF|ExAC_SAS_MAF|CLIN_SIG|SOMATIC|PHENO|PUBMED|MOTIF_NAME|MOTIF_POS|HIGH_INF_POS|MOTIF_SCORE_CHANGE|LoF|LoF_filter|LoF_flags|LoF_info'
CSQ_FMT = CSQ_FMT.split("|")
print CSQ_FMT
def getINFO(info_string):
    infolist = info_string.split(';')
    infodict = {}
    for kv in infolist:
        kv = kv.split('=')
        if len(kv) == 2:
            k, v = kv
            if k in infodict:
                infodict[k].append(v)
            else:
                infodict[k] = [v]
    return infodict
def searchCANONICAL(CSQLIST):
    for csq in CSQLIST:
        csq = csq.split("|")
        CSQ = dict(zip(CSQ_FMT, csq))
        if CSQ["CANONICAL"] == "YES":
            #return csq
            return "YES"
    return None
def displayCSQ(head, _dict):
    for item in head:
        print "{}:{}\t".format(item, _dict[item]),

def processBenign():
    hand = gzip.open(gnomad_file, 'rb')
    outfil = csv.writer(open("benigh.tsv", 'wb'), delimiter="\t")
    outfil.writerow(head)
    for l in hand:
        if l.startswith("#"):
            continue
        llist = l.strip().split("\t")
        Chr, Pos, ID, ref, alt = llist[0:5]
        if len(ref) != 1:
            continue
        INFO_string = llist[7]
        info_dict = getINFO(INFO_string)
        #print info_dict
        #CSQ = searchCANONICAL(info_dict["CSQ"])
        CSQ = info_dict["CSQ"][0].split("|")
        CSQ = dict(zip(CSQ_FMT, CSQ))
        if CSQ["Consequence"] != "missense_variant":
            continue
        AA = CSQ["Amino_acids"].split("/")
        res = [""]*19
        res[0] = info_dict["Gene.refGene"][0]
        res[1] = 'b37'
        res[2] = Chr
        res[3] = Pos
        res[4] = ref
        res[5] = alt.split(",")[0]
        res[6] = AA[0]
        res[7] = AA[1]
        res[8] = CSQ["Protein_position"]
        res[9] = CSQ["ENSP"]
        res[10] = CSQ["SWISSPROT"]
        res[11] = info_dict["avsnp147"][0]
        res[12] = "gnomAD_exome"
        res[13] = "-"
        res[14] = "-"
        res[15] = "-"
        res[16] = "{}_{}_{}_{}".format(Chr,Pos,ref,res[5])
        res[17] = "{}_{}_{}_{}".format(res[10],res[8],res[6],res[7])
        res[18] = info_dict["gnomAD_exome_ALL"][0]
        outfil.writerow(res)

def processObs():
    hand = gzip.open("gnomAD.Exon.CDS.vcf.gz", 'rb')
    out = open("gnomad.Exon.Rare.HighConf.vcf", 'wb')
    for l in hand:
        if l.startswith("#"):
            out.write(l)
        else:
            llist = l.strip().split("\t")
            if llist[6] != "PASS":
                continue
            INFO_string = llist[7]
            info_dict = getINFO(INFO_string)
            AC = int(info_dict["AC"][0])
            if AC > 251 or AC == 0:
                continue
            CSQ = info_dict["CSQ"][0].split("|")
            CSQ = dict(zip(CSQ_FMT, CSQ))
            if CSQ["Consequence"] not in ["synonymous_variant", "missense_variant"]:
                continue
            out.write(l)


def main():
    processObs()

#main()
