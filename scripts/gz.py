import re
import sys
import os
import glob
import shutil
import util
import config

from subprocess import Popen, PIPE


@util.timing
def compress(input, output):
    with open(output, "w") as f:
        process = Popen(['gzip', '-c', input], stdout=f, stderr=PIPE)
        process.communicate()


def compress_mprof(input, output, mprof_folder):
    command = "\"" + ' '.join([config.memusg_path, 'gzip', '-c',
                        input, '>', output]) + "\""
    print(command)
    process = Popen(command, shell=True, stderr=PIPE)
    process.communicate()
    util.safe_move('memory_profile.csv', os.path.join(
        mprof_folder, os.path.basename(input)+'-gz-compress.csv'))


@util.timing
def decompress(input, output):
    with open(output, "w") as f:
        process = Popen(['gzip', '-dc', input], stdout=f, stderr=PIPE)
        process.communicate()


def decompress_mprof(input, output, mprof_folder):
    command = "\"" + ' '.join([config.memusg_path, 'gzip', '-dc',
                        input, '>', output])  + "\""
    print(command)
    process = Popen(command, shell=True, stderr=PIPE)
    process.communicate()
    util.safe_move('memory_profile.csv', os.path.join(mprof_folder, os.path.splitext(os.path.basename(
        input))[0]+'-gz-decompress.csv'))


if __name__ == "__main__":
    pass
