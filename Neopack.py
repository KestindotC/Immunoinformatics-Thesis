# Neoantigen Prediction Script
import sys
import time
import subprocess
import multiprocessing as mp
from pyensembl import EnsemblRelease

# release 77 uses human reference genome GRCh38

db = EnsemblRelease(87)


desc = '''
        ===========================================================================
        Neopack are written for calling netMHCpan3.0 with following filtering setup.
        HLA alleles were iterated for each input samples.

        Required input:
        - fastafiles (should be *.fasta files, in the same work directory)
        - hlalist (python list)

        e.g. Neopack.start("~/test.fasta",["HLA-A01:01","HLA-C06:07"])
        
        ** Warning: Now netMHCpan realared para and path is hardcoded.
           Please Check to publish

        Author : Kestin
        ===========================================================================
        '''
    
print(desc)

def netmhcpancalling(args):
    ## infilename must have extension ".fasta"
    ## the hlaSting is a comma seperated list of MHC alleles. Remove the *s from the allele names.""
    infilename,hla = args
    netpanexe = "/tempwork173/kestin/tool/netMHCpan-3.0/netMHCpan"
    attr = hla[:-3]
    outfilename = infilename[:-6] + attr + ".txt"
    errfilename = infilename[:-6] + attr + ".err"

    errfile = open(errfilename, "w") # open err and std files
    outfile = open(outfilename, "w")
    
    cmdlist = [netpanexe,
               "-a",hla,
               "-l","9",
               infilename]

    print("netpan start running...")
    execode = subprocess.call(cmdlist, stdout = outfile, stderr = errfile)
    errfile.close()  #close files
    outfile.close()

    print("netpan finish, coumsuming time: " + time.time() + "sec.")
    netmhcpanreportfilter(outfilename)



def netmhcpanreportfilter(outfile):
    filtertmp = []
    with open(outfile,'r') as tmpreport:
        for line in tmpreport:
            tmpstr_list=[]
            if "<= SB" in line or "<= WB" in line:
                tmpstr_list = line.split(" ")
                tmpstr_list = filter(None,tmpstr_list)
                filtertmp.append([tmpstr_list[9],pyensembl(tmpstr_list[10]),tmpstr_list[12]])
    
    with open(outfile[:-4] + "_filtered.txt","w") as out:
		for item in filtertmp:
			print >> out, "\t".join(item)




def pyensembl(protein_id):
    geneinfo = db.gene_by_protein_id(protein_id)
    gene_symbl = str(geneinfo).split(",")[1].split("=")[1]
    return gene_symbl




## Main progress

def start(fasta,hla = "HLA-A01:01"):
    args = ((fasta,hlas) for hlas in hla)
    pool = mp.Pool(processes=6)
    pool.map(netmhcpancalling,args)


if __name__ == '__main__':
    hla = ["HLA-A01:01","HLA-B15:01"]
    start("~/test.fasta",hla)
    print("finish testing.")
