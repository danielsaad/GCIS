import re
import sys
import os
import glob
import shutil
import util

from subprocess import Popen, PIPE

@util.timing
def compress(input,output):
    with open(output, "w") as f:
        process = Popen(['gzip', '-c', input], stdout=f,stderr=PIPE)
        process.communicate()
        
@util.timing
def decompress(input,output):
    with open(output, "w") as f:
        process = Popen(['gzip', '-dck', input], stdout=f,stderr=PIPE)
        process.communicate()

if __name__ == "__main__":
    pass
