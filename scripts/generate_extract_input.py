import os
import random
import sys


number_of_substrings = 1000
substring_size = [1, 10, 100, 1000, 10000]


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
            with open(output_file_path + '.' + str(ss)+'_query', 'w') as extract_file:
                for k in range(0, number_of_substrings):
                    # Generates number_of_substrings random queries of
                    # length ss and outputs it in output_file
                    beg = random.randint(0, input_size-ss-1)
                    end = beg+ss-1
                    print(beg, end, file=extract_file)


''' 
    argv[1] = Input File
    argv[2] = Output File 
'''
if __name__ == "__main__":
    generate_extract_input(sys.argv[1], sys.argv[2])
