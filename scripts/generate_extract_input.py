import os
import random
import sys
number_of_substrings = 1000
substring_size = [1,10,100,1000,10000]

def generate_extract_input(input_file_path):
    with open(input_file_path,'r') as input_file:
        input_file.seek(0, os.SEEK_END)
        input_size = input_file.tell()
        for ss in substring_size:
            with open(sys.argv[1]+'.'+str(ss)+'_query','w') as extract_file:
                for k in range(0,number_of_substrings):
                    beg = random.randint(0,input_size-ss-1)
                    end = beg+ss-1
                    print(beg,end,file=extract_file)



''' argv[1] = Input File '''
if __name__ == "__main__":
    generate_extract_input(sys.argv[1])
