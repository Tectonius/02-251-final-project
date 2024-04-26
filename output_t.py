from final_wrapper import run_clustal_test
import os
import pickle

input_folder_path = r".\PCA Seqs"

tests_path = r".\PCA Seqs\sabre_seqs"

refs_path = "Sup Refs"

score_pairs = []

# run_clustal_test(file1, output1_name, output2_name, ref_alignment_path)

for file in os.listdir(tests_path):
    if '.dnd' not in file:
        ref = os.path.join(refs_path, file)
        inp_path = os.path.join(tests_path, file)
        score_pairs.append(run_clustal_test(inp_path, ref))
with open('scored_pairs.pkl', 'wb') as f:
    pickle.dump(score_pairs, f)