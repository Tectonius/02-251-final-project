import Bio
from Bio import SeqIO
from Bio import AlignIO

def seq_file_reader(filename):
    '''
    Input: file name of a sequence (no matter alignment or sequence set).
    Output: a list of strings of that sequence/alignment.

    Note: file should be in the same folder as code file.
    '''
    if filename[-4:] == '.aln':
        seqs = AlignIO.read(filename, 'clustal')
    else:
        try:
            seqs = SeqIO.parse(filename, 'fasta')
        except:
            print('File is neither fasta nor aln file!!')

    sequences = []
    for rec in seqs:
        sequences.append(str(rec.seq))
    return sequences
