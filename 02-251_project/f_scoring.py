from Bio import AlignIO
import blosum as bl
import os
import time
import math


def wait_until(somepredicate, timeout, output_file, period=0.25):
    mustend = time.time() + timeout
    while time.time() < mustend:
        if somepredicate(output_file):
            return True
        time.sleep(period)
    return False


def aligner_done(file_path):
    try:
        with open(file_path, "r") as f:
            return len(f.readlines()) != 0
    except Exception as E:
        return False


def f_scoring(path_to_sequence_file, path_to_aligner, params, ref_alignment, aligner):
    '''
    path_to_sequence_file: path to the unaligned sequence in fasta format
    path_to_clustal: path to the aligner
    parameters: a dictionary of numbers that is put into MSA software
    ref_alignment: list of strings representing the reference alignment sequences
    aligner: a string--name of msa aligner used
    s_matrix: Dictionary that represents the scoring matrix for SP Score
    gop: gap openning penalty
    gep: gap extension penalty
    Output: the SP score difference between MSA output by the softwares and the reference.
    '''
    s_matrix = bl.BLOSUM(62)
    gop = 1
    gep = 1
    # gapopen, gapext, maxdiv, transweight
    parameters = {'gapopen': params[0], 'gapext': params[1], 'maxdiv': params[2], 'transweight': params[3]}
    if aligner == 'clustalw':
        cwAln = path_to_aligner + ' -infile=' + path_to_sequence_file
        cwAln += ' -outfile="aligned.aln"'
        for key in parameters.keys():
            cwAln += f' -{key}="{parameters[key]}"'
        try:
            os.remove('aligned.aln')
        except Exception as E:
            pass
        os.system(cwAln)
        time.sleep(0.5)
        wait_until(aligner_done, 100000000, "aligned.aln", period=0.1)
        aln = AlignIO.read('aligned.aln', 'clustal')
    elif aligner == 'dialigntx':
        daAln = path_to_aligner + ' ..' + path_to_sequence_file
        daAln += ' -outfile="aligned"'
        for key in parameters.keys():
            daAln += f' -{key}="{parameters[key]}"'
        os.system(r'Remove-Item -Path "aligned.aln"')  # if dialigntx even makes .aln files
        os.system(daAln)
        time.sleep(0.5)
        wait_until(aligner_done, 100000000, "aligned.aln", period=0.1)
        aln = AlignIO.read('aligned.aln', 'clustal')
    alignment = []
    for rec in aln:
        alignment.append(rec.seq)
    ref_alignment = [ref_alignment[19], ref_alignment[16], ref_alignment[7]]
    sp_diff = SP_Score_Diff(alignment, ref_alignment, s_matrix, gop, gep)
    return sp_diff / len(ref_alignment[0])

def run_aligner(path_to_sequence_file, path_to_aligner, params, output_alignment_name, aligner):
    # gapopen, gapext, maxdiv, transweight
    parameters = {'gapopen': params[0], 'gapext': params[1], 'maxdiv': params[2], 'transweight': params[3]}
    if aligner == 'clustalw':
        cwAln = path_to_aligner + ' -infile=' + path_to_sequence_file
        cwAln += f' -outfile="{output_alignment_name}.aln"'
        for key in parameters.keys():
            cwAln += f' -{key}="{parameters[key]}"'
        try:
            os.remove('aligned.aln')
        except Exception as E:
            pass
        os.system(cwAln)
        time.sleep(0.5)
        wait_until(aligner_done, 100000000, "aligned.aln", period=0.1)
    elif aligner == 'dialigntx':
        daAln = path_to_aligner + ' ..' + path_to_sequence_file
        daAln += f' -outfile="{output_alignment_name}"'
        for key in parameters.keys():
            daAln += f' -{key}="{parameters[key]}"'
        os.system(r'Remove-Item -Path "aligned.aln"')  # if dialigntx even makes .aln files
        os.system(daAln)
        time.sleep(0.5)
        wait_until(aligner_done, 100000000, "aligned.aln", period=0.1)



def SP_Score_Diff(alignment, ref_alignment, s_matrix, gop, gep):
    if (len(alignment) != len(ref_alignment)):
        return "Produce extra or less strings in alignment"
    l = len(alignment[0])
    result = 0
    for i in range(len(alignment)):
        for j in range(i + 1, len(alignment)):
            score = 0
            score_r = 0
            si = alignment[i]
            sj = alignment[j]
            si_r = ref_alignment[i]
            sj_r = ref_alignment[j]
            si_p = ''
            sj_p = ''
            si_r_p = ''
            sj_r_p = ''
            # Compute pair score for alignment
            for k in range(len(si)):
                if not (si[k] == '-' and sj[k] == '-'):
                    si_p += si[k]
                    sj_p += sj[k]
                if si[k] != '-' and sj[k] != '-':
                    s = s_matrix[si[k].upper()][sj[k].upper()]
                    if math.isinf(s): s = 0
                    score += s
            for ig in computeGapIntervals(si_p).union(computeGapIntervals(sj_p)):
                score += gCost(ig, gop, gep)
            # Compute pair score for ref_alignment
            for k in range(len(si_r)):
                if not (si_r[k] == '-' and sj_r[k] == '-'):
                    si_r_p += si_r[k]
                    sj_r_p += sj_r[k]
                if si_r[k] != '-' and sj_r[k] != '-':
                    s = s_matrix[si_r[k].upper()][sj_r[k].upper()]
                    if math.isinf(s): s = 0
                    score_r += s
            for ig in computeGapIntervals(si_r_p).union(computeGapIntervals(sj_r_p)):
                score_r += gCost(ig, gop, gep)
            result += (score_r - score)
    return abs(result)


def computeGapIntervals(sequence):
    gap_list = set()
    cur_gap = []
    for i in range(len(sequence)):
        if sequence[i] == '-' and len(cur_gap) == 0:
            cur_gap.append(i)
        elif sequence[i] != '-' and len(cur_gap) == 1:
            gap_list.add((cur_gap[0], i))
            cur_gap = []
    if len(cur_gap) != 0:
        gap_list.add((cur_gap[0], len(sequence)))
    return gap_list


def gCost(ig, gop, gep):
    return gop + gep * (ig[1] - ig[0])

def test_f_scoring(path_to_sequence_file, path_to_aligner, ref_alignment, aligner):
    '''
    path_to_sequence_file: path to the unaligned sequence in fasta format
    path_to_clustal: path to the aligner
    parameters: a dictionary of numbers that is put into MSA software
    ref_alignment: list of strings representing the reference alignment sequences
    aligner: a string--name of msa aligner used
    s_matrix: Dictionary that represents the scoring matrix for SP Score
    gop: gap openning penalty
    gep: gap extension penalty
    Output: the SP score difference between MSA output by the softwares and the reference.
    '''
    s_matrix = bl.BLOSUM(62)
    gop = 1
    gep = 1
    # gapopen, gapext, maxdiv, transweight
    if aligner == 'clustalw':
        cwAln = path_to_aligner + ' -infile=' + path_to_sequence_file
        cwAln += ' -outfile="test_aligned.aln"'
        try:
            os.remove('test_aligned.aln')
        except Exception as E:
            pass
        os.system(cwAln)
        time.sleep(0.5)
        wait_until(aligner_done, 100000000, "test_aligned.aln", period=0.1)
        aln = AlignIO.read('test_aligned.aln', 'clustal')
    alignment = []
    for rec in aln:
        alignment.append(rec.seq)
    ref_alignment = [ref_alignment[19], ref_alignment[16], ref_alignment[7]]
    sp_diff = SP_Score_Diff(alignment, ref_alignment, s_matrix, gop, gep)
    return sp_diff / len(ref_alignment[0])

aligner = "clustalw"
path_to_aligner = r'..\ClustalW2\clustalw2.exe'
path_to_sequence_file = r'"2022_03_07_first_3.fasta"'
ref_alignment = r'"2022_03_07_A.fasta"'

print(test_f_scoring(path_to_sequence_file, path_to_aligner, ref_alignment, aligner))