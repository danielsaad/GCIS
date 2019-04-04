import os
import sys
import gcis
import tempfile
"""
    Produces SACA data considering: 
    1) Yuta SAIS
    2) GCIS suffix array construction underdecoding

    The results are aggregated with the following:
    1) GCIS decode + Yuta SAIS
    2) SA constrution directly from GCIS decompression
    3) Yuta SAIS from plain text
"""

def parse_time(s):
    pass
def build_gcis_dictionaries(input_folder,output_folder):
    print("Creating output folder",output_folder)
    os.makedirs(output_folder,exist_ok=True)
    
    input_files = [os.path.join(input_folder,f) for f in os.listdir(input_folder)]
    output_files = [os.path.join(output_folder,f+'.gcis') for f in os.listdir(input_folder)]
    gcis.compress_gc_is_parallel(input_files,output_files)


def decode_and_build(dictionary_folder):
    
    time_list = []

    # obtain .gcis file list
    files = [os.path.join(dictionary_folder,f) 
        for f in os.listdir(dictionary_folder)
        if os.path.splitext(f)[0] == '.gcis']

    '''
        Create temporary file, build the suffix
        array from decompression and store it on
        the temporary file then delete temporary file.

        Stdout is captured and parsed to obtain total time.
    '''    
    for f in files:
        # Create temporary file
        _,tmp_file_path = tempfile.mkstemp()
        out,err = gcis.decompress_gc_is(f,tmp_file_path)
        time_list.append(parse_time(out.decode()))
        # Remove temporary file
        os.remove(tmp_file_path)

    return time_list

def save(input_folder,results_db,results_gs,results_sa):
    pass

"""
    argv[1] = folder containing experiments
    argv[2] = result folder containing dictionary files
"""
if __name__== '__main__':
    print('Input Folder =',sys.argv[1])
    print('Output Folder =',sys.argv[2])

    print('Building GCIS Dictionaries')
    build_gcis_dictionaries(sys.argv[1],sys.argv[2])
    print("Evaluating GCIS decode + Yuta SAIS. ")
    results_db = decode_and_build(sys.argv[2])
    print("Evaluating GCIS SACA.")
    results_gs = gcis_saca(sys.argv[2])
    print("Evaluating Yuta SAIS")
    results_sa = sais(sys.argv[1])

    # Save time spent and experiment names into .csv file
    save(sys.argv[1],results_db,results_gs,results_sa)