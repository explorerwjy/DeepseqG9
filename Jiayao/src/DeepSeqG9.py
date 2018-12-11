import pandas as pd
import csv
import numpy
import gzip
#import matplotlib as mpl
#import matplotlib.pyplot as plt
import math
import sys

class ManipulateData:
    def __init__(self):
        self.std_header = ['genename', 'assembly', 'chr', 'pos', 'ref', 'alt', 'ref_aa', 'alt_aa', 'swissprot_pos', 'ensp', 'swissprot', 'dbsnp', 'source', 'source_id', 'target', 'descr', 'var_id', 'swissprot_var_id', 'af']
        return 
    def mergeTable(self, GeneFile, hgmd, clinvar, cancer, AF = 1e-4):
        out = csv.writer(open("hgmd_clinvar_cancer_combined_pathogenic_variants.tsv", 'wb'), delimiter = "\t")
        out.writerow(self.std_header)
        self.hgmd_genes = self.reduce_var_gene(hgmd)
        self.cancer_genes = self.reduce_var_gene(cancer)
        self.clinvar_genes = self.reduce_var_gene(clinvar) 
        reader = csv.reader(open(GeneFile), delimiter="\t")
        header = reader.next()
        res = {}
        for row in reader:
            gene = row[1]
            gene_vars = {}
            for dat_name, dataset in zip(['hgmd', 'clinvar', 'cancer'], [self.hgmd_genes, self.clinvar_genes, self.cancer_genes]):
                if gene not in dataset:
                    continue
                for _var, row in dataset[gene].items():
                    if row[14] != '1':
                        continue
                    if float(row[18]) > AF:
                        continue
                    if _var not in gene_vars:
                        gene_vars[_var] = row
                    else:
                        row_tmp = gene_vars[_var][:]
                        #print row_tmp
                        row_tmp[12] = row_tmp[12] + ";" + row[12]
                        row_tmp[13] = row_tmp[13] + ";" + row[13]
                        gene_vars[_var] = row_tmp
            for K in sorted(gene_vars.keys()):
                out.writerow(gene_vars[K])
    def reduce_var_gene(self, _file):
        res = {}
        reader = csv.reader(open(_file), delimiter="\t")
        header = reader.next()
        idx_var_id = header.index("var_id")
        for row in reader:
            gene = row[0]
            var_id = row[idx_var_id]
            if gene not in res:
                res[gene] = {}
            res[gene][var_id] = row
        return res
    def MakeBenighVars(self, srcFil, PathFil):
        Genes = sorted(list(set(pd.read_csv(PathFil, delimiter="\t")["genename"].values)))
        print len(Genes)
    def ProcessPDB2Uniprot(self, mapfile):
        PrList = []
        last_PDB, last_Uni = None, None
        tmp_pr, tmp_pdb = Gene2PDB("Fake"), myPDB("fake","fake") 
        #reader = csv.reader(open(mapfile, 'rb'), delimiter="\t")
        reader = open(mapfile, 'rb')
        for l in reader:
            row = l.strip().split()
            try:
                PDB, Chain, trash, AA, PDB_Coord, Uni, aa, Uni_Coord = row
            except:
                continue
            ID = PDB + "." + Chain
            if last_Uni != Uni:
                #print last_Uni, Uni
                tmp_pr.PDBs.append(tmp_pdb)
                tmp_pr.PDBs = tmp_pr.PDBs[1:]
                PrList.append(tmp_pr)
                last_Uni = Uni
                tmp_pr = Gene2PDB(Uni)
                tmp_pdb = myPDB(PDB, Chain)
                tmp_pdb.start = int(Uni_Coord)
                #print Uni, ID
            elif last_Uni == Uni and ID != last_PDB:
                #print PDB, last_PDB
                tmp_pr.PDBs.append(tmp_pdb)
                last_PDB = ID 
                tmp_pdb = myPDB(PDB, Chain)
                tmp_pdb.start = int(Uni_Coord)
            else:
                tmp_pdb.end = int(Uni_Coord)
        tmp_pr.PDBs.append(tmp_pdb)
        tmp_pr.PDBs = tmp_pr.PDBs[1:]
        PrList.append(tmp_pr)
        PrList = PrList[1:]
        for gene in PrList:
            #print gene.Uni
            gene.process()
            #gene.display(mode=1)
            #print
            gene.display(mode=2)
    def Reformat(self, mapfile, outfile):
        PrList = []
        out = csv.writer(open(outfile, 'wb'), delimiter="\t")
        last_PDB, last_Uni = None, None
        tmp_pr, tmp_pdb = Gene2PDB("Fake"), myPDB("fake","fake") 
        #reader = csv.reader(open(mapfile, 'rb'), delimiter="\t")
        reader = open(mapfile, 'rb')
        for l in reader:
            row = l.strip().split()
            try:
                PDB, Chain, trash, AA, PDB_Coord, Uni, aa, Uni_Coord = row
                out.writerow([PDB, Chain, trash, AA, PDB_Coord, Uni, aa, Uni_Coord])
            except:
                continue
                
class myPDB:
    def __init__(self, PDB, chain):
        self.PDB = PDB
        self.chain = chain
        self.ID = PDB+"."+chain
        self.start = 1 
        self.end = 1 
class Gene2PDB:
    def __init__(self, Uni):
        self.Uni = Uni
        self.PDBs = []
    def display(self, mode=1):
        if mode == 1:
            for pdb in self.PDBs:
                print self.Uni, pdb.ID, pdb.start, pdb.end
        elif mode == 2:
            for pdb in self.finalist:
                #print self.Uni, pdb.ID, pdb.start, pdb.end
                print "%s\t%s\t%s\t%s"%(self.Uni, pdb.ID, pdb.start, pdb.end)
    def process(self):
        self.finalist = []
        try:
            segments = [self.PDBs[0]]
        except IndexError:
            #print self.Uni, print self.PDBs
            sys.stderr.write("{}\n".format(self.Uni))
            return
        for PDB in self.PDBs[1:]:
            keep = True
            for seg in segments:
                #print seg.ID, seg.start, seg.end, PDB.ID, PDB.start, PDB.end
                if PDB.end < seg.start or PDB.start > seg.end:
                    keep = True
                elif seg.start == PDB.start and seg.end == PDB.end:
                    segments.remove(seg)
                    keep = False
                elif seg.start <= PDB.start and seg.end >= PDB.end:
                    keep = False
                elif (PDB.start <= seg.start and PDB.end > seg.end) or (PDB.start < seg.start and PDB.end >= seg.end):
                    #print seg.ID + " To remove"
                    segments.remove(seg)
                    #segments.append(PDB)
                    keep = True
                else:
                    if (PDB.end - PDB.start) > (seg.end - seg.start):
                        keep = True
                        segments.remove(seg)
                    else:
                        keep = False
            if keep:
                segments.append(PDB)
        self.finalist = segments
        return

ins = ManipulateData()

