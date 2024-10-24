
from subprocess import Popen
from subprocess import PIPE
import os.path
import report2csv
import util

'''
    Compress a single file creating a GCIS dictionary
'''


gcis_path = '../bin/gcis'
gcis64_path = '../bin/gcis-64'
gcis_statistics_path = '../bin/gcis-memory'
gcis64_statistics_path = '../bin/gcis-64-memory'

SIZE_LIMIT_32_BIT = 2**31 - 1


@util.timing
def compress(input, output, codec='-ef'):
    path = gcis_path if util.calc_file_size(
        input) <= SIZE_LIMIT_32_BIT else gcis64_path
    process = Popen([path, '-c', input, output, codec])
    p = process.communicate()
    return p


@util.timing
def decompress(input, output, codec='-ef'):
    path = gcis_path if util.calc_file_size(
        input) <= SIZE_LIMIT_32_BIT else gcis64_path
    process = Popen([path, '-d', input, output, codec])
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
        path = gcis_path if util.calc_file_size(
            inf) <= SIZE_LIMIT_32_BIT else gcis64_path
        if(os.path.isfile(ouf)):
            print(os.path.basename(ouf), 'already exists: skipping')
        else:
            process_list.append(Popen([path, '-c', inf, ouf]))
    for p in process_list:
        outerr_l.append(p.communicate())

    return outerr_l


'''
    Build SA under decompressing using GCIS
'''


def decompress_saca(input, output, codec='-s8b'):
    path = gcis_path if util.calc_file_size(
        input) <= SIZE_LIMIT_32_BIT else gcis64_path
    process = Popen([path, '-s', input, output, codec])
    return process.communicate()


'''
    Build SA + LCP under decompressing using GCIS
'''


def decompress_saca_lcp(input, output, codec='-s8b'):
    path = gcis_path if util.calc_file_size(
        input) <= SIZE_LIMIT_32_BIT else gcis64_path

    process = Popen([path, '-l', input, output, codec])
    return process.communicate()


'''
    Decode a compressed GCIS file and return the time used on
    decompression (without saving file).
'''


def decompress_gc_is(input, output, codec='-s8b'):
    path = gcis_path if util.calc_file_size(
        input) <= SIZE_LIMIT_32_BIT else gcis64_path
    p = Popen([path, '-d', input, output, codec],
              stdout=PIPE, stderr=PIPE)
    out, _ = p.communicate()
    lines = out.decode("utf-8").split('\n')
    decompress_time = 0
    for l in lines:
        if(l.startswith('time:')):
            decompress_time = l.split()[1]
    return (float(decompress_time))


def compress_gc_is_statistics(input, output, codec='-s8b'):
    path = gcis_statistics_path if util.calc_file_size(
        input) <= SIZE_LIMIT_32_BIT else gcis64_statistics_path
    process = Popen([path, '-c', input, output, codec])
    return process.communicate()


def decompress_gc_is_statistics(input, output, codec='-s8b'):
    path = gcis_statistics_path if util.calc_file_size(
        input) <= SIZE_LIMIT_32_BIT else gcis64_statistics_path
    process = Popen([path, '-d', input, output, codec])
    return process.communicate()


''' 
    Extract substrings from the dictionary file accordling
    to the query file
'''


def gcis_extract(inf, query_file):
    path = gcis_path if util.calc_file_size(
        inf) <= SIZE_LIMIT_32_BIT else gcis64_path

    command_gcis = [path, '-e', inf, query_file]
    p = Popen(command_gcis, stdout=PIPE)
    return p.communicate()


def gcis_self_index_extract(inf, index_file, query_file):
    command_gcis = ['../bin/gcis-ef', '-x', inf, index_file, query_file]
    p = Popen(command_gcis, stdout=PIPE)
    return p.communicate()


def gcis_self_index_locate(inf, index_file, query_file):
    command_gcis = ['../bin/gcis-ef', '-p', inf, index_file, query_file]
    p = Popen(command_gcis, stdout=PIPE)
    return p.communicate()


def gcis_build_index(inf, ouf):
    command = ['../bin/gcis-ef', '-i', inf, ouf]
    p = Popen(command, stdout=PIPE, stderr=PIPE)
    return p.communicate()


def gcis_build_index_statistics(inf, ouf):
    command = ['../bin/gcis-ef-memory', '-i', inf, ouf]
    p = Popen(command, stdout=PIPE, stderr=PIPE)
    return p.communicate()


if __name__ == "__main__":
    pass
