import sys
from run_experiments import run_experiments
from generate_graphs import generate_graphs

''' argv[1] = experiments folder
    argv[2] = results folder '''
if __name__ == '__main__':
    run_experiments(sys.argv[1],sys.argv[2])
    generate_graphs(sys.argv[2])
