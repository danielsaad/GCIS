from subprocess import Popen
from subprocess import PIPE
import os.path


def compress_repair(input):
    process = Popen(['../bin/repair','-i', input],stdout=PIPE,stderr=PIPE)
    process.communicate()

def decompress_repair(input):
    process = Popen(['../bin/despair','-i',input],stdout=PIPE,stderr=PIPE)
    process.communicate()

def compress_repair_statistics(input,output):
    process = Popen(['../bin/repair-memory','-i',input],stdout=PIPE,stderr=PIPE)
    process.communicate()


def decompress_repair_statistics(input,output):
    process = Popen(['../bin/despair-memory','-i',input],stdout=PIPE,stderr=PIPE)
    process.communicate()


def compress_repair_navarro(input):
    process = Popen(['../bin/repair-navarro', input],stdout=PIPE,stderr=PIPE)
    process.communicate()

def decompress_repair_navarro(input):
    process = Popen(['../bin/despair-navarro', input],stdout=PIPE,stderr=PIPE)
    process.communicate()


if __name__ == "__main__":
    pass
