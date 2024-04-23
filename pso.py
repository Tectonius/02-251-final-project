import pyswarms as ps
import numpy as np
from f_scoring import f_scoring
import time
from f_scoring import run_aligner

path_to_aligner = r'..\ClustalW2\clustalw2.exe'
path_to_sequence_file = r'"ox_8_seqs.fasta"'
ref_alignment = r'"ox_8_ref_align.fasta"'
iter = 0
true_start = time.time()
best_params = None
best_params_score = None

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

def optimize_clustal(path_to_aligner, path_to_sequence_file, ref_alignment, logname):
    def clustalscore(params):
        # input: params is a numpy.ndarray
        # set of inputs of shape :code:`(n_particles, dimensions)`
        #
        # output:  numpy.ndarray of size :code:`(n_particles, )`
        output_array = []
        global iter
        start = time.time()
        for particle in params:
            output_array.append(f_scoring(path_to_sequence_file, path_to_aligner, particle, ref_alignment))
        append_to_file(f"###\niteration: {iter}\ntotal time: {time.time() - true_start}\niteration time:{time.time() - start}\nparameters:\n{params}\noutput scores:\n{output_array}\n###\n", f"{logname}.txt")
        iter += 1
        global best_params
        global best_params_score
        if best_params_score == None:
            best_params = params[output_array.index(min(output_array))]
            best_params_score = output_array[output_array.index(min(output_array))]
        elif output_array[output_array.index(min(output_array))] < best_params_score:
            best_params = params[output_array.index(min(output_array))]
        return np.array(output_array)
    def optimize_params(output_alignment_names):
        # clustalw param order:
        # gapopen, gapext, maxdiv, transweight
        bounds = (np.array([0, 0, 0, 0]),np.array([100, 100, 100, 1]))
        options = {'c1': 0.5,
                   'c2': 0.5,
                   'w': 0.9}
        optimizer = ps.single.GlobalBestPSO(n_particles=5, dimensions=4, options=options, bounds=bounds)
        best_position, _ = optimizer.optimize(clustalscore, iters=10)
        run_aligner(path_to_sequence_file, path_to_aligner, best_params, output_alignment_names)
        return best_position
    print(f"covid params: {optimize_params('ox_final_alignment')}")

optimize_clustal(path_to_aligner, path_to_sequence_file, ref_alignment)