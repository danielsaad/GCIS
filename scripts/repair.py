import subprocess
from subprocess import Popen
from subprocess import PIPE
import util
import os.path
import config
import shutil

repair_path = '../bin/repair'
despair_path = '../bin/despair'

repair_memory_path = '../bin/repair-memory'
despair_memory_path = '../bin/despair-memory'

repair_navarro_path = '../external/repair-navarro/repair'
despair_navarro_path = '../external/repair-navarro/despair'

irepair_navarro_path = '../external/repair-navarro/irepair'
idespair_navarro_path = '../external/repair-navarro/idespair'


@util.timing
def compress_repair(input, output_folder):
    process = subprocess.run(
        [repair_path, '-i', input], stdout=PIPE, stderr=PIPE)
    input_basename = os.path.basename(input)
    input_folder = os.path.dirname(input)
    util.safe_move(os.path.join(input_folder, input_basename+'.prel'),
                   os.path.join(output_folder, input_basename+'.prel'))
    util.safe_move(os.path.join(input_folder, input_basename+'.seq'),
                   os.path.join(output_folder, input_basename+'.seq'))


def compress_mprof(input, output_folder, mprof_folder):
    process = subprocess.run(
        [config.memusg_path, repair_path, '-i', input], stdout=PIPE, stderr=PIPE)
    input_basename = os.path.basename(input)
    input_folder = os.path.dirname(input)
    util.safe_move(os.path.join(input_folder, input_basename+'.prel'),
                   os.path.join(output_folder, input_basename+'.prel'))
    util.safe_move(os.path.join(input_folder, input_basename+'.seq'),
                   os.path.join(output_folder, input_basename+'.seq'))
    util.safe_move('memory_profile.csv', os.path.join(mprof_folder, os.path.basename(
        input)+'-repair-compress.csv'))


@util.timing
def decompress_repair(input, output_folder):
    input_basename = os.path.basename(input)
    process = subprocess.run(
        [despair_path, '-i', input], stdout=PIPE, stderr=PIPE)
    util.safe_move(input+'.u', os.path.join(output_folder, input_basename))


def decompress_mprof(input, output_folder, mprof_folder):
    input_basename = os.path.basename(input)
    process = subprocess.run(
        [config.memusg_path, despair_path, '-i', input], stdout=PIPE, stderr=PIPE)
    util.safe_move(input+'.u', os.path.join(output_folder, input_basename))
    util.safe_move('memory_profile.csv', os.path.join(mprof_folder, os.path.basename(
        input)+'-repair-decompress.csv'))


def compress_repair_statistics(input, output):
    process = Popen([repair_memory_path, '-i', input],
                    stdout=PIPE, stderr=PIPE)
    process.communicate()


def decompress_repair_statistics(input, output):
    process = Popen([despair_memory_path, '-i', input],
                    stdout=PIPE, stderr=PIPE)
    process.communicate()


''' 
    Compresses the input with Navarro's version of repair
    and returns the (stdout,stderr) pair 
'''


def compress_repair_navarro(input, output_folder):
    p = subprocess.run([repair_navarro_path, input], stdout=PIPE, stderr=PIPE)
    if (all([os.path.isfile(input+ext) for ext in ['.C', '.R']])):
        for f in [input + ext for ext in ['.C', '.R']]:
            util.safe_move(f, os.path.join(output_folder, os.path.basename(f)))
    else:
        print(F'Repair Navarro compression failed for {input}')
    return p


def decompress_repair_navarro(input_path, output_folder):
    p = subprocess.run([despair_navarro_path, input_path],
                       stdout=PIPE, stderr=PIPE)
    decompressed_file = os.path.splitext(input_path)[0]
    ouf = os.path.join(output_folder, os.path.basename(decompressed_file))
    util.safe_move(decompressed_file, ouf)
    return p


def compress_repair_navarro_int(input, output_folder):
    p = subprocess.run([irepair_navarro_path, input], stdout=PIPE, stderr=PIPE)
    if (all([os.path.isfile(input+ext) for ext in ['.C', '.R']])):
        for f in [input + ext for ext in ['.C', '.R']]:
            util.safe_move(f, os.path.join(output_folder, os.path.basename(f)))
    else:
        print(F'Repair Navarro compression failed for {input}')
    return p


def decompress_repair_navarro_int(input_path, output_folder):
    p = subprocess.run([idespair_navarro_path, input_path],
                       stdout=PIPE, stderr=PIPE)
    # decompressed_file = os.path.splitext(input_path)[0]
    decompressed_file = input_path
    print(f'input_path = {input_path} decompressed_file = {decompressed_file}')
    ouf = os.path.join(output_folder, os.path.basename(decompressed_file))
    util.safe_move(decompressed_file, ouf)
    return p


if __name__ == "__main__":
    pass
