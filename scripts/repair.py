import subprocess
from subprocess import Popen
from subprocess import PIPE
import util
import os.path

repair_path = '../bin/repair'
despair_path = '../bin/despair'

repair_memory_path = '../bin/repair-memory'
despair_memory_path = '../bin/despair-memory'

@util.timing
def compress_repair(input,output_folder):
    process = subprocess.run([repair_path, '-i', input], stdout=PIPE, stderr=PIPE)
    input_basename = os.path.basename(input)
    input_folder = os.path.dirname(input)
    os.rename(os.path.join(input_folder,input_basename+'.prel'),os.path.join(output_folder,input_basename+'.prel'))
    os.rename(os.path.join(input_folder,input_basename+'.seq'),os.path.join(output_folder,input_basename+'.seq'))

@util.timing
def decompress_repair(input,output_folder):
    input_basename = os.path.basename(input)
    process = subprocess.run([despair_path, '-i', input], stdout=PIPE, stderr=PIPE)
    os.rename(input+'.u',os.path.join(output_folder,input_basename))

def compress_repair_statistics(input, output):
    process = Popen([repair_memory_path,'-i', input], stdout=PIPE, stderr=PIPE)
    process.communicate()


def decompress_repair_statistics(input, output):
    process = Popen([despair_memory_path, '-i', input],
                    stdout=PIPE, stderr=PIPE)
    process.communicate()


''' 
    Compresses the input with Navarro's version of repair
    and returns the (stdout,stderr) pair 
'''


def compress_repair_navarro(input):
    process = Popen(['../bin/repair-navarro', input], stdout=PIPE, stderr=PIPE)
    return process.communicate()


def decompress_repair_navarro(input):
    process = Popen(['../bin/despair-navarro', input],
                    stdout=PIPE, stderr=PIPE)
    process.communicate()


if __name__ == "__main__":
    pass
