import re
import sys
import os
import glob
import shutil
import time

from subprocess import Popen, PIPE


def compress_gz(input,output):
    with open(output, "w") as f:
        process = Popen(['gzip', '-c', input], stdout=f,stderr=PIPE)
        process.communicate()

def decompress_gz(input,output):
    with open(output, "w") as f:
        process = Popen(['gzip', '-dck', input], stdout=f,stderr=PIPE)
        process.communicate()

if __name__ == "__main__":
    pass
