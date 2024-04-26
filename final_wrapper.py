import pca_t
import f_scoring
import argparse
from ref_sp import reference_based_SP_Score
import os

opt_param_dict = {'covid_seqs': [78.94887853, 65.3188784, 95.17580108, 0.99500936],
                  'ox_seqs': [9.05710228, 36.24392932, 40.41719035, 0.48147639],
                  'sabre_seqs': [9.5614024, 38.05434494, 46.38328272, 0.25971668]}

path_to_aligner = r'..\ClustalW2\clustalw2.exe'


def run_clustal(file, output_name):
    closest_label = pca_t.closest_label(file)
    params = opt_param_dict[closest_label]
    f_scoring.run_aligner(file, path_to_aligner, params, output_name)




def run_clustal_test(file1, ref_alignment_path):
    closest_label = pca_t.closest_label(file1)
    params = opt_param_dict[closest_label]
    f_scoring.run_aligner(file1, path_to_aligner, params, 't1')
    f_scoring.run_aligner(file1, path_to_aligner, {}, 't2')
    opt, reg = reference_based_SP_Score(ref_alignment_path, 't1.aln'), reference_based_SP_Score(ref_alignment_path, 't2.aln')
    try:
        os.remove('t1.aln')
        os.remove('t2.aln')
    except:
        pass
    return opt, reg

def main():
    parser = argparse.ArgumentParser(description="Run Clustal function from the command line")
    parser.add_argument("file", help="Input file")
    parser.add_argument("output_name", help="Output file name")
    args = parser.parse_args()
    run_clustal(args.file, args.output_name)


if __name__ == "__main__":
    main()
