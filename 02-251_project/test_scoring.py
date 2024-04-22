import subprocess
from Bio import AlignIO

# Define input and output file paths
input_fasta = r'"2022_03_07_first_3.fasta"'
output_aln = r'"aligned.aln"'

import os
import time
start = time.time()
def wait_until(somepredicate, timeout, output_file, period=0.25):
  mustend = time.time() + timeout
  while time.time() < mustend:
    if somepredicate(output_file): return True
    time.sleep(period)
  return False

def aligner_done(file_path):
    with open(file_path, "r") as f:
        return len(f.readlines()) != 0 # This would give length of files.

# clustal_command = [r'..\ClustalW2\clustalw2.exe', input_fasta, '-outfile='+output_aln, '-gapext='+'"5.5"']
# print(" ".join(clustal_command))
# os.system(" ".join(clustal_command))
# time.sleep(0.5)
# print(aligner_done(output_aln[1:-1]))
# wait_until(aligner_done, 100000000, output_aln[1:-1], period=0.1)
# print(time.time()-start)
# k=3
print('Remove-Item -Path "aligned.aln"')
print(r'Remove-Item -Path "aligned.aln"')
os.system('Remove-Item -Path "aligned.aln"')



