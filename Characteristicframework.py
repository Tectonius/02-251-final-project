import numpy as np
from collections import Counter


def calculate_aa_percentages(sequences):
    aa_categories = {
        'positive': ['R', 'H', 'K'],
        'negative': ['D', 'E'],
        'polar': ['S', 'T', 'N', 'Q'],
        'special': ['C', 'U', 'G', 'P'],
        'hydrophobic': ['A', 'V', 'I', 'L', 'M', 'F', 'Y', 'W']
    }
    t_length = sum(len(seq) for seq in sequences)
    aa_counts = Counter("".join(sequences))

    percentages = {}
    for category, aas in aa_categories.items():
        category_count = sum(aa_counts[aa] for aa in aas)
        percentages[category] = (category_count / t_length) * 100

    return percentages


# def kimura_distances(sequences):

# with open("dna_sequences.txt", "w") as file:
# for sequence in sequences:
# file.write(sequence + "\n")

def kimura_distance(sequences):
    distances = []

    for i, seq1 in enumerate(sequences):
        for seq2 in sequences[i + 1:]:
            d_sites = sum(a != b for a, b in zip(seq1, seq2) if a != '-' and b != '-')
            t_sites = sum(1 for a, b in zip(seq1, seq2) if a != '-' and b != '-')
            p = d_sites / t_sites

            D = -np.log(1 - p - 0.2 * p ** 2)
            distances.append(D)

    mean_distance = np.mean(distances)
    std_distance = np.std(distances)

    return [mean_distance, std_distance]


def extract_characteristics(sequences):
    a1 = len(sequences)
    a2 = np.mean([len(seq) for seq in sequences])
    a3 = np.std([len(seq) for seq in sequences])

    b1, b2 = kimura_distance(sequences)

    cs = calculate_aa_percentages(sequences)

    return [a1, a2, a3, b1, b2, cs['positive'], cs['negative'], cs['polar'], cs['special'], cs['hydrophobic']]
