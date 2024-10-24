import os
import random
import sys


number_of_substrings = 10000
substring_size = [1, 5, 10, 25, 50, 100]


'''
    Generates number_of_substrings random interval queries of 
    substring_size length
'''


def generate_extract_input(input_file_path, output_file_path):
    # Open input file
    with open(input_file_path, 'r') as input_file:
        # Get input size in bytes
        input_file.seek(0, os.SEEK_END)
        input_size = input_file.tell()
        # for each substring length
        for ss in substring_size:
            with open(output_file_path + '.' + str(ss)+'_extract', 'w') as extract_file:
                print(number_of_substrings,ss, file=extract_file)
                intervals = [(x, x+ss-1)
                             for x in [random.randint(0, input_size-ss-1) for k in range(number_of_substrings)]]
                for i in intervals:
                    print(i[0], i[1], file=extract_file)


''' 
    argv[1] = Input File
    argv[2] = Output File 
'''
if __name__ == "__main__":
    random.seed(42)
    generate_extract_input(sys.argv[1], sys.argv[2])
