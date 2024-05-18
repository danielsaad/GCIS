import subprocess
from subprocess import Popen
from subprocess import PIPE
import os.path
import report2csv
import util
import config
import shutil

'''
    Compress a single file creating a GCIS dictionary
'''


gcis_path = '../bin/gcis'
gcis64_path = '../bin/gcis-64'
gcis_statistics_path = '../bin/gcis-memory'
gcis64_statistics_path = '../bin/gcis-64-memory'
gcis_alpha_path = '../bin/gcis_alpha'
gcis_repair_hybrid_path = '../bin/gcis_repair_hybrid'
gcis_repair_hybrid_2_path = '../bin/gcis_repair_hybrid_2'
mprof_path = '/home/danielsaad/.local/bin/mprof'
python_path = 'python3.6'

SIZE_LIMIT_32_BIT = 2**31 - 1


@util.timing
def compress(input, output, codec='-s8b'):
    path = gcis_path if util.calc_file_size(
        input) <= SIZE_LIMIT_32_BIT else gcis64_path
    if (output.endswith('.gcis') and path == gcis64_path):
        output += '64'
    print(path, '-c', input, output, codec)
    process = Popen([path, '-c', input, output, codec],
                    stdout=PIPE, stderr=PIPE)
    p = process.communicate()
    return p


@util.timing
def compress_statistics(input, output, codec='-s8b'):
    path = gcis_statistics_path if util.calc_file_size(
        input) <= SIZE_LIMIT_32_BIT else gcis64_path
    if (output.endswith('.gcis') and path == gcis64_statistics_path):
        output += '64'
    print(path, '-c', input, output, codec)
    process = Popen([path, '-c', input, output, codec])
    p = process.communicate()
    print(p)
    return p


@util.timing
def gcis_alpha_compress(input, output):
    path = gcis_alpha_path
    print(path, '-c', input, output)
    process = Popen([path, '-c', input, output],
                    stdout=PIPE, stderr=PIPE)
    p = process.communicate()
    return p


@util.timing
def gcis_repair_hybrid_compress(input, output, level=0):
    path = gcis_repair_hybrid_path
    cmd = [path, '-c', input, output, str(level)]
    print(' '.join(cmd))
    process = Popen(cmd,
                    stdout=PIPE, stderr=PIPE)
    p = process.communicate()
    return p


@util.timing
def gcis_repair_hybrid_2_compress(input, output):
    path = gcis_repair_hybrid_2_path
    cmd = [path, '-c', input, output]
    print(' '.join(cmd))
    process = Popen(cmd,
                    stdout=PIPE, stderr=PIPE)
    p = process.communicate()
    return p


@util.timing
def gcis_repair_hybrid_decompress(input, output, level=0):
    path = gcis_repair_hybrid_path
    cmd = [path, '-d', input, output, str(level)]
    print(' '.join(cmd))
    process = Popen(cmd,
                    stdout=PIPE, stderr=PIPE)
    p = process.communicate()
    return p


@util.timing
def gcis_alpha_decompress(input, output):
    path = gcis_alpha_path
    print(path, '-d', input, output)
    process = Popen([path, '-d', input, output],
                    stdout=PIPE, stderr=PIPE)
    p = process.communicate()
    return p


def compress_mprof(input, output, mprof_folder, codec='-s8b'):
    path = gcis_path if util.calc_file_size(
        input) <= SIZE_LIMIT_32_BIT else gcis64_path
    if (output.endswith('.gcis') and path == gcis64_path):
        output += '64'
    process = Popen([config.memusg_path, path, '-c', input, output, codec])
    p = process.communicate()
    util.safe_move('memory_profile.csv', os.path.join(mprof_folder, os.path.basename(
        input)+'-gcis-compress.csv'))
    return p


@util.timing
def decompress(input, output, codec='-s8b'):
    path = gcis_path if input.endswith('.gcis') else gcis64_path
    print(' '.join([path, '-d', input, output, codec]))
    process = Popen([path, '-d', input, output, codec])
    p = process.communicate()
    return p


def decompress_mprof(input, output, mprof_folder, codec='-s8b'):
    path = gcis_path if input.endswith('.gcis') else gcis64_path
    process = Popen([config.memusg_path, path, '-d', input, output, codec])
    p = process.communicate()
    util.safe_move('memory_profile.csv', os.path.join(mprof_folder, os.path.splitext(os.path.basename(
        input))[0]+'-gcis-decompress.csv'))
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
        if (os.path.isfile(ouf)):
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
    path = gcis_path if input.endswith('gcis') else gcis64_path
    process = Popen([path, '-s', input, output, codec])
    return process.communicate()


'''
    Build SA + LCP under decompressing using GCIS
'''


def decompress_saca_lcp(input, output, codec='-s8b'):
    path = gcis_path if input.endswith('gcis') else gcis64_path

    process = Popen([path, '-l', input, output, codec])
    return process.communicate()


'''
    Decode a compressed GCIS file and return the time used on
    decompression (without saving file).
'''


def decompress_gc_is(input, output, codec='-s8b'):
    path = gcis_path if input.endswith('gcis') else gcis64_path
    p = Popen([path, '-d', input, output, codec],
              stdout=PIPE, stderr=PIPE)
    out, _ = p.communicate()
    lines = out.decode("utf-8").split('\n')
    decompress_time = 0
    for l in lines:
        if (l.startswith('time:')):
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


def gcis_extract(inf, query_file, codec='-ef'):
    path = gcis_path if inf.endswith('gcis') else gcis64_path
    command_gcis = [path, '-e', inf, query_file, codec]
    p = subprocess.run(command_gcis, stdout=PIPE, stderr=PIPE)
    return p


''' Parse gcis_extract stdout and returns the time in microseconds '''


def parse_extract(p):
    stdout = p.stdout.decode('utf-8')
    result = 'ERROR'
    for line in stdout.split('\n'):
        if line.startswith('Mean time'):
            result = line.split(': ')[-1]
    if (p.stderr):
        print(p.stderr)
    return result


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
