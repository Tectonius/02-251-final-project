import compare_characteristics
import f_scoring
import argparse


opt_param_dict = {'covid_seqs': [], 'ox_seqs': [], 'sabre_seqs': []}

path_to_aligner = r'..\ClustalW2\clustalw2.exe'

def run_clustal(file, output_name):
    closest_label = compare_characteristics.closest_label(file)
    params = opt_param_dict[closest_label]
    f_scoring.run_aligner(file, path_to_aligner, params, output_name)

def main():
    parser = argparse.ArgumentParser(description="Run Clustal function from the command line")
    parser.add_argument("file", help="Input file")
    parser.add_argument("output_name", help="Output file name")
    args = parser.parse_args()

    run_clustal(args.file, args.output_name)

if __name__ == "__main__":
    main()

