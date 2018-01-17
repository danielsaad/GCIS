
from subprocess import Popen
import os.path
import report2csv


def compress_gc_is(input,output):
    process = Popen(['../bin/gc-is-codec', '-c', input,output])
    process.communicate()


def decompress_gc_is(input,output):
    process = Popen(['../bin/gc-is-codec', '-d', input,output])
    process.communicate()


def compress_gc_is_statistics(input,output):
    process = Popen(['../bin/gc-is-codec-memory', '-c', input,output])
    process.communicate()



def decompress_gc_is_statistics(input,output):
    process = Popen(['../bin/gc-is-codec-memory', '-d', input,output])
    process.communicate()

if __name__ == "__main__":
    pass
