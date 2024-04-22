import os
from Bio import AlignIO
from sklearn.decomposition import PCA
from sklearn.neighbors import NearestNeighbors
from Characteristicframework import extract_characteristics


def read_data_and_extract_features(folder_path):
    data = []
    labels = []
    for label in os.listdir(folder_path):
        label_path = os.path.join(folder_path, label)
        if os.path.isdir(label_path):
            for file_name in os.listdir(label_path):
                file_path = os.path.join(label_path, file_name)
                # Read sequences from FASTA file
                alignment = AlignIO.read(file_path, "fasta")
                # Extract features using Characteristicframework.py
                features = extract_characteristics(alignment)
                data.append(features)
                labels.append(label)
    return data, labels


def find_nearest_neighbor(input_data, data):
    # Perform PCA
    pca = PCA(n_components=2)
    pca.fit(data)
    reduced_data = pca.transform(data)

    # Reduce input data to 2 components
    input_reduced = pca.transform(input_data)

    # Find nearest neighbor
    nn = NearestNeighbors(n_neighbors=1)
    nn.fit(reduced_data)
    nearest_neighbor_index = nn.kneighbors(input_reduced, return_distance=False)[0][0]

    return nearest_neighbor_index


def main(input_folder_path, input_sequence):
    # Read data and extract features
    data, labels = read_data_and_extract_features(input_folder_path)

    # Process input sequence using Characteristicframework.py
    input_features = extract_characteristics(input_sequence)

    # Find nearest neighbor
    nearest_neighbor_index = find_nearest_neighbor(input_features, data)
    nearest_neighbor_label = labels[nearest_neighbor_index]

    return nearest_neighbor_label


# Example usage
input_folder_path = ""
input_sequence = "path/to/input_sequence.fasta"
nearest_label = main(input_folder_path, input_sequence)
print("Nearest neighbor label:", nearest_label)
