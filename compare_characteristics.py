import os
from sklearn.decomposition import PCA
from sklearn.neighbors import NearestNeighbors
from Characteristicframework import extract_characteristics
from seq_file_reader import seq_file_reader
from sklearn.preprocessing import StandardScaler


def read_data_and_extract_features(folder_path):
    data = []
    labels = []
    for label in os.listdir(folder_path):
        label_path = os.path.join(folder_path, label)
        if os.path.isdir(label_path):
            for file_name in os.listdir(label_path):
                file_path = os.path.join(label_path, file_name)
                alignment = seq_file_reader(file_path)

                features = extract_characteristics(alignment)
                data.append(features)
                labels.append(label)
    return data, labels


def find_nearest_neighbor(input_data, data):
    scaler = StandardScaler()
    scaler.fit(data)

    data_scaled = scaler.transform(data)
    input_data_scaled = scaler.transform([input_data])

    # Perform PCA
    pca = PCA(n_components=2)
    pca.fit(data_scaled)
    reduced_data = pca.transform(data_scaled)

    input_reduced = pca.transform(input_data_scaled)

    nn = NearestNeighbors(n_neighbors=5)
    nn.fit(reduced_data)
    nearest_neighbor_index = nn.kneighbors(input_reduced, return_distance=False)[0][0]

    return nearest_neighbor_index


def closest_label(input_sequence):
    # Read data and extract features
    data, labels = read_data_and_extract_features("./PCA Seqs")

    # Process input sequence using Characteristicframework.py
    input_features = extract_characteristics(input_sequence)

    # Find nearest neighbor
    nearest_neighbor_index = find_nearest_neighbor(input_features, data)
    nearest_neighbor_label = labels[nearest_neighbor_index]

    return nearest_neighbor_label


# Example usage
input_folder_path = "./PCA Seqs"
input_sequence = "./sup_002.fasta"
nearest_label = closest_label( input_sequence)
print("Nearest neighbor label:", nearest_label)
