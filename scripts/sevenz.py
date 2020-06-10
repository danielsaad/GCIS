import os
import util
from subprocess import Popen, PIPE

@util.timing
def compress(input,output):
    process = Popen(['7za', 'a', output, input, "-mmt=1"],stdout=PIPE,stderr=PIPE)
    process.communicate()

@util.timing
def decompress(input,output):
    outputdir = os.path.dirname(output)
    process = Popen(['7za', 'e', input, "-y -o" + outputdir,"-mmt=1"],stdout=PIPE,stderr=PIPE)
    process.communicate()


if __name__ == "__main__":
    pass
