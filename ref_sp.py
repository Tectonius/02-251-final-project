import Bio
from Bio import AlignIO
from seq_file_reader import seq_file_reader


def reference_based_SP_Score(ref_alignment_name, alignment_name):
    '''
    This function assess the percentage of reference alignment recovered by the
    alignment of interst.

    ref_alignment: a list of strings containing the reference alignment.
    alignment: a list of strings containing the alignment of interest.

    output: A decimal between 0 and 1.
    '''
    ref_alignment = seq_file_reader(ref_alignment_name)

    alignment = seq_file_reader(alignment_name)

    if len(ref_alignment) != len(alignment):
        print('Have different number of rows!!!')

    ref_len = len(ref_alignment[0])
    ali_len = len(alignment[0])
    rows = len(ref_alignment)
    col = 0
    total_pairs = 0
    total_recovered = 0

    while col < ref_len and col < ali_len:
        recorder = dict()
        for i in range(rows):
            for j in range(i, rows):
                pair = (ref_alignment[i][col], ref_alignment[j][col])
                if pair in recorder:
                    recorder[pair] += 1
                elif (pair[1], pair[0]) in recorder:
                    recorder[(pair[1], pair[0])] += 1
                else:
                    recorder[pair] = 1
                total_pairs += 1
        for i in range(rows):
            for j in range(i, rows):
                pair = (alignment[i][col], alignment[j][col])
                if pair in recorder:
                    recorder[pair] -= 1
                    if recorder[pair] == 0: recorder.pop(pair)
                    total_recovered += 1
                elif (pair[1], pair[0]) in recorder:
                    recorder[(pair[1], pair[0])] -= 1
                    if recorder[(pair[1], pair[0])] == 0:
                        recorder.pop((pair[1], pair[0]))
                    total_recovered += 1
        col += 1

    return total_recovered / total_pairs
