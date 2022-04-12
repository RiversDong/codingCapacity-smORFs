from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import os

def merged(inFile, output):
    #fasta = ["Euk_Ara_coding.fasta", "Euk_Hum_coding.fasta", "Euk_Mou_coding.fasta"]
    fasta = ["Euk_Ara_noncoding.fasta", "Euk_Hum_noncoding.fasta", "Euk_Mou_noncoding.fasta"]
    path = "/home/chuand/small_orf/data/data_sets_fromJiang/Euk/"

    records_list = []
    for i in fasta:
        i_file = os.path.join(path, i)
        records = SeqIO.parse(i_file, "fasta")
        for j in records:
            j_id = str(j.id)
            j_seq = str(j.seq)
            rec = SeqRecord(Seq(j_seq), id=j_id, description="")
            records_list.append(rec)
    #SeqIO.write(records_list, "/home/chuand/small_orf/data/mergerd_eukaryotes.fasta", "fasta")
    SeqIO.write(records_list, "/home/chuand/small_orf/data/mergerd_eukaryotes_noncoding.fasta", "fasta")

def splitFile(infile,outfile):
    file1 = open(outfile+"100", "w")
    file2 = open(outfile+"200", "w")
    file3 = open(outfile+"300", "w")
    records = SeqIO.parse(infile, "fasta")
    for i in records:
        iid = str(i.id)
        iseq = str(i.seq)
        iLen = len(iseq)
        if 1<= iLen <=100:
            file1.write(">"+iid+"\n")
            file1.write(iseq + "\n")
        if 100 <= iLen <= 200:
            file2.write(">"+iid+"\n")
            file2.write(iseq + "\n")
        if 200 <= iLen <= 300:
            file3.write(">"+iid+"\n")
            file3.write(iseq + "\n")
    file1.close()
    file2.close()
    file3.close()

if __name__ == "__main__":
    merged(["Euk_Ara_coding.fasta", "Euk_Hum_coding.fasta", "Euk_Mou_coding.fasta"],\
            "/home/chuand/small_orf/data/mergerd_eukaryotes_coding.fasta")    
    merged(["Euk_Ara_noncoding.fasta", "Euk_Hum_noncoding.fasta", "Euk_Mou_noncoding.fasta"],\
            "/home/chuand/small_orf/data/mergerd_eukaryotes_noncoding.fasta")    

    splitFile("/home/chuand/small_orf/data/mergerd_eukaryotes.fasta",\
            "/home/chuand/small_orf/data/splitFasta/split_eu_coding")
    splitFile("/home/chuand/small_orf/data/trainAndTest/test/Pro_coding.fasta",\
            "/home/chuand/small_orf/data/splitFasta/split_pro_coding")

    splitFile("/home/chuand/small_orf/data/mergerd_eukaryotes_noncoding.fasta",\
            "/home/chuand/small_orf/data/splitFasta/split_eu_noncoding")
    splitFile("/home/chuand/small_orf/data/trainAndTest/test/Pro_nocoding.fasta",\
            "/home/chuand/small_orf/data/splitFasta/split_pro_noncoding")

