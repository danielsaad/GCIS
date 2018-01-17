import os

from subprocess import Popen, PIPE


def compress_7z(input,output):
    process = Popen(['7za', 'a', output, input, "-mmt=1"],stdout=PIPE,stderr=PIPE)
    process.communicate()

def decompress_7z(input,output):
    outputdir = os.path.split(input)[0]
    process = Popen(['7za', 'e', input, "-o" + outputdir,"-mmt=1"],stdout=PIPE,stderr=PIPE)
    process.communicate()


if __name__ == "__main__":
    pass
