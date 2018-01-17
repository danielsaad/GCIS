import os
import sys


def run_extract(input_folder_path,output_folder_path):
    input_file_paths = os.listdir(input_folder_path)
    for f in input_file_paths:
        print("Compressing",os.path.basename(f):
        process = Popen(['../bin/gc-is-codec', '-c', input,output])
        process.communicate()

#argv[1] contains the path to the experiments folder
#argv[2] contains the path to the output folder

if __name__ == "__main__":
    run_extract(sys.argv[1],sys.argv[2])
