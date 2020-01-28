
from subprocess import Popen
from subprocess import PIPE
import os.path
import report2csv

'''
    Compress a single file creating a GCIS dictionary
'''


def compress_gc_is(input, output,codec='-ef'):
    process = Popen(['../bin/gc-is-codec', '-c', input, output,codec])
    p = process.communicate()
    return p


'''
    Compute several GCIS dictionaries from input list to
    output list. If an element of output list already
    exists in filesystem, the compression is not performed
'''


def compress_gc_is_parallel(input_list, output_list):
    process_list = []
    outerr_l = []
    for inf, ouf in zip(input_list, output_list):
        if(os.path.isfile(ouf)):
            print(os.path.basename(ouf), 'already exists: skipping')
        else:
            process_list.append(Popen(['../bin/gc-is-codec', '-c', inf, ouf]))
    for p in process_list:
        outerr_l.append(p.communicate())

    return outerr_l


'''
    Build SA under decompressing using GCIS
'''


def decompress_saca(input, output):
    process = Popen(['../bin/gc-is-codec', '-s', input, output])
    return process.communicate()


'''
    Build SA + LCP under decompressing using GCIS
'''


def decompress_saca_lcp(input, output):
    process = Popen(['../bin/gc-is-codec', '-l', input, output])
    return process.communicate()


'''
    Decode a compressed GCIS file and return the time used on
    decompression (without saving file).
'''


def decompress_gc_is(input, output):
    p = Popen(['../bin/gc-is-codec', '-d', input, output],
              stdout=PIPE, stderr=PIPE)
    out, _ = p.communicate()
    lines = out.decode("utf-8").split('\n')
    for l in lines:
        if(l.startswith('time:')):
            decompress_time = l.split()[1]
    return (float(decompress_time))


def compress_gc_is_statistics(input, output,codec='-ef'):
    process = Popen(['../bin/gc-is-codec-memory', '-c', input, output,codec])
    return process.communicate()


def decompress_gc_is_statistics(input, output):
    process = Popen(['../bin/gc-is-codec-memory', '-d', input, output])
    return process.communicate()


''' 
    Extract substrings from the dictionary file accordling
    to the query file
'''


def gcis_extract(inf, query_file):

    command_gcis = ['../bin/gc-is-codec', '-e', inf, query_file]
    p = Popen(command_gcis,stdout=PIPE)
    return p.communicate()


def gcis_build_index(inf,ouf):
    command = ['../bin/gcis-ef','-i',inf,ouf]
    p = Popen(command,stdout=PIPE,stderr=PIPE);
    return p.communicate();

if __name__ == "__main__":
    pass
