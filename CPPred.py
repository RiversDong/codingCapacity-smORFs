import sys
import os

infile = sys.argv[1]
out = sys.argv[2]

cmd = "python /home/chuand/small_orf/bin1/CPPred-sORF/CPPred-sORF.py -i {0} -hex /home/chuand/small_orf/bin1/CPPred-sORF/Hexamer/Integrated_Hexamer.txv -r /home/chuand/small_orf/bin1/CPPred-sORF/Model/Integrated.range -mol /home/chuand/small_orf/bin1/CPPred-sORF/Model/Integrated.model -spe Integrated -o {1}".format(infile, out)
os.system(cmd)
