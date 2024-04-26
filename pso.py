import pyswarms as ps
import numpy as np
from f_scoring import f_scoring
import time
from f_scoring import run_aligner


def wait_until(somepredicate, timeout, output_file, period=0.25):
    mustend = time.time() + timeout
    while time.time() < mustend:
        if somepredicate(output_file):
            return True
        time.sleep(period)
    return False


def aligner_done(file_path):
    with open(file_path, "r") as f:
        return len(f.readlines()) != 0


def append_to_file(text, file_name):
    try:
        with open(file_name, 'a') as file:
            file.write(text + '\n')
    except FileNotFoundError:
        with open(file_name, 'w') as file:
            file.write(text + '\n')


def optimize_clustal(path_to_aligner, path_to_sequence_file, ref_alignment, logname, final_output_name):
    def clustalscore(params):
        # input: params is a numpy.ndarray
        # set of inputs of shape :code:`(n_particles, dimensions)`
        #
        # output:  numpy.ndarray of size :code:`(n_particles, )`
        output_array = []
        global iter_
        start = time.time()
        for particle in params:
            output_array.append(f_scoring(path_to_sequence_file, path_to_aligner, particle, ref_alignment))
        append_to_file(
            f"###\niteration: {iter_}\ntotal time: {time.time() - true_start}\niteration time:{time.time() - start}\nparameters:\n{params}\noutput scores:\n{output_array}\n###\n",
            f"{logname}.txt")
        iter_ += 1
        global best_params
        global best_params_score
        if best_params_score is None:
            best_params = params[output_array.index(min(output_array))]
            best_params_score = output_array[output_array.index(min(output_array))]
        elif output_array[output_array.index(min(output_array))] < best_params_score:
            best_params = params[output_array.index(min(output_array))]
        return np.array(output_array)

    def optimize_params(output_alignment_names):
        # clustalw param order:
        # gapopen, gapext, maxdiv, transweight
        bounds = (np.array([0, 0, 0, 0]), np.array([100, 100, 100, 1]))
        options = {'c1': 0.5,
                   'c2': 0.5,
                   'w': 0.9}
        optimizer = ps.single.GlobalBestPSO(n_particles=25, dimensions=4, options=options, bounds=bounds)
        best_position, _ = optimizer.optimize(clustalscore, iters=50)
        run_aligner(path_to_sequence_file, path_to_aligner, best_params, output_alignment_names)
        return best_position

    print(f"covid params: {optimize_params(final_output_name)}")


def run_pso(path_to_aligner, path_to_sequence_file_list, ref_alignment_list, log_name_list, final_output_name_list):
    for i in range(len(path_to_sequence_file_list)):
        optimize_clustal(path_to_aligner, path_to_sequence_file_list[i], ref_alignment_list[i], log_name_list[i],
                         final_output_name_list[i])


path_to_aligner = r'..\ClustalW2\clustalw2.exe'
path_to_sequence_file_list = ['"1a0cA_1a0dA.fasta"', '"1aab_ref1.fasta"', '"BB11001.fasta"', '"sup_002.fasta"']
ref_alignment_list = ['1a0cA_1a0dA_ref_align.fasta', '1aab_ref1_ref_align.fasta', 'BB11001_ref_align.fasta',
                      'sup_002_ref_align.fasta']  # shouldn't include quotes
log_name_list = ["1a0cA_1a0dA_log", "1aab_ref1_log", "BB11001_log", "sup_002_log"]
final_output_name_list = ["1a0cA_1a0dA_final", "1aab_ref1_final", "BB11001_final", "sup_002_final"]
iter_ = 0
true_start = time.time()
best_params = None
best_params_score = None

run_pso(path_to_aligner, path_to_sequence_file_list, ref_alignment_list, log_name_list, final_output_name_list)
