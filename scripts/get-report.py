import os
import subprocess
import repair
import gcis


'''
    Run GCIS in order to get report.log statistics.
    Receives an experiment folder and stores the resulting reports
    in an output folder.
'''
def run_gcis(experiments_folder,output_folder):
    files = [os.path.join(f,experiments_folder) for f in os.listdir(experiments_folder)]
    for f in files:
        ouf = os.path.join(output_folder)
        gcis.decompress_gc_is_statistics(f,)


if __name__ == "__main__":
    pass