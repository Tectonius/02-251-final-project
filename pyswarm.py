import pyswarms as ps
import numpy as np
from f_scoring import f_scoring
import time

aligner = "clustalw"
path_to_aligner = r'..\ClustalW2\clustalw2.exe'
path_to_sequence_file = r'"2022_03_07_first_3.fasta"'
ref_alignment = r'"2022_03_07_first_3_A.fasta"'
iter = 0
true_start = time.time()

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

def optimize_clustal(aligner, path_to_aligner, path_to_sequence_file, ref_alignment):
    def clustalscore(params):
        # input: params is a numpy.ndarray
        # set of inputs of shape :code:`(n_particles, dimensions)`
        #
        # output:  numpy.ndarray of size :code:`(n_particles, )`
        output_array = []
        global iter
        start = time.time()
        for particle in params:
            output_array.append(f_scoring(path_to_sequence_file, path_to_aligner, particle, ref_alignment, aligner))
        append_to_file(f"###\niteration: {iter}\ntotal time: {time.time() - true_start}\niteration time:{time.time() - start}\nparameters:\n{params}\noutput scores:\n{output_array}\n###\n", "log_2022_03_07_first_3.txt")
        iter += 1
        return np.array(output_array)
    def optimize_params(aligner):
        if aligner == "clustalw":
            # clustalw param order:
            # gapopen, gapext, maxdiv, transweight
            bounds = (np.array([0, 0, 0, 1]),np.array([100, 100, 100, 1]))
            options = {'c1': [1, 2, 3],
                       'c2': [1, 2, 3],
                       'w': [2, 3, 5],
                       'k': [5, 10, 15],
                       'p': 1}
            optimizer = ps.single.GlobalBestPSO(n_particles=10, dimensions=4, options=options, bounds=bounds)
            best_position, _ = optimizer.optimize(clustalscore, iters=50)
            return best_position
    print(f"covid params: {optimize_params('clustalw')}")




