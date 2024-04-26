import os
import pickle
from sklearn.decomposition import PCA
from sklearn.neighbors import NearestNeighbors
from Characteristicframework import extract_characteristics
from seq_file_reader import seq_file_reader
from sklearn.preprocessing import StandardScaler
import matplotlib.pyplot as plt
import numpy as np


def read_data_and_extract_features(folder_path):
    data = []
    labels = []
    for idx, label in enumerate(os.listdir(folder_path)):  # enumerate to get index for labels
        label_path = os.path.join(folder_path, label)
        if os.path.isdir(label_path):
            for file_name in os.listdir(label_path):
                file_path = os.path.join(label_path, file_name)
                alignment = seq_file_reader(file_path)
                features = extract_characteristics(alignment)
                data.append(features)
                labels.append(label_path.split("/")[-1])  # Use index as label instead of string

    # Pickle the extracted data
    with open('final_extracted_data.pkl', 'wb') as f:
        pickle.dump((data, labels), f)

    return data, labels


def load_extracted_data():
    # Unpickle the extracted data
    with open('final_extracted_data.pkl', 'rb') as f:
        data, labels = pickle.load(f)
    return data, labels


def find_nearest_neighbors(input_data, data, n_neighbors=5):
    scaler = StandardScaler()
    scaler.fit(data)

    data_scaled = scaler.transform(data)
    input_data_scaled = scaler.transform([input_data])

    # Perform PCA
    pca = PCA(n_components=2)
    pca.fit(data_scaled)
    reduced_data = pca.transform(data_scaled)

    input_reduced = pca.transform(input_data_scaled)

    nn = NearestNeighbors(n_neighbors=n_neighbors)
    nn.fit(reduced_data)
    nearest_neighbor_indices = nn.kneighbors(input_reduced, return_distance=False)[0]

    return nearest_neighbor_indices, reduced_data, input_reduced

def most_common(lst):
    return max(set(lst), key=lst.count)

def closest_label(input_sequence):
    # Check if extracted data has been pickled
    data, labels = load_extracted_data()

    # Process input sequence using Characteristicframework.py
    input_features = extract_characteristics(input_sequence)

    # Find nearest neighbors and get reduced data for visualization
    nearest_neighbor_indices, reduced_data, input_reduced = find_nearest_neighbors(input_features, data)
    nearest_neighbor_labels = [labels[i] for i in nearest_neighbor_indices]

    # # Plot PCA
    # plt.figure(figsize=(8, 6))
    # for i, label in enumerate(set(labels)):
    #     indices = np.where(np.array(labels) == label)
    #     plt.scatter(reduced_data[indices, 0], reduced_data[indices, 1], label=f'Label {label}', alpha=0.7)
    # plt.scatter(reduced_data[nearest_neighbor_indices, 0], reduced_data[nearest_neighbor_indices, 1],
    #             marker='x', color='red', label='Nearest Neighbors')
    # plt.scatter(input_reduced[0, 0], input_reduced[0, 1], marker='o', color='black', label='Input Sequence')
    # plt.xlabel('Principal Component 1')
    # plt.ylabel('Principal Component 2')
    # plt.title('PCA Visualization')
    # plt.legend()
    # plt.grid(True)
    # plt.show()

    return most_common(nearest_neighbor_labels)


# Example usage
# input_folder_path = "./PCA Seqs"
# input_sequence = "./PCA Seqs/ox_seqs/12s103"
# nearest_labels = closest_label(input_folder_path, input_sequence)
# print("Nearest neighbor labels:", nearest_labels)
