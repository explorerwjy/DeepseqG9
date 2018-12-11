from DeepSeqG9 import *
#PDB2Uniprot = "../dat/test.map"
#PDB2Uniprot = "../dat/pdbsws_res.txt"
#ins.Reformat(PDB2Uniprot, "../dat/PDB2SWISS.tsv")
PDB2Uniprot = "../dat/PDB2SWISS.sort.tsv"
ins.ProcessPDB2Uniprot(PDB2Uniprot)
