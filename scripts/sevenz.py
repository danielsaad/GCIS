import os
import util
import config
import shutil
from subprocess import Popen, PIPE


BUFFER_SIZE_MB = 1000


@util.timing
def compress(input, output):
    process = Popen(['7za', '-md={}m'.format(BUFFER_SIZE_MB), 'a', output, input, "-mmt=1"],
                    stdout=PIPE, stderr=PIPE)
    process.communicate()


def compress_mprof(input, output, mprof_folder):
    process = Popen([config.memusg_path, '7za', '-md={}m'.format(BUFFER_SIZE_MB), 'a', output, input, "-mmt=1"],
                    stdout=PIPE, stderr=PIPE)
    process.communicate()
    util.safe_move('memory_profile.csv', os.path.join(mprof_folder, os.path.basename(
        input)+'-7zip-compress.csv'))


@util.timing
def decompress(input, output):
    outputdir = os.path.dirname(output)
    process = Popen(['7za', 'e', input, '-o' + outputdir, '-y',
                     "-mmt=1"], stdout=PIPE, stderr=PIPE)
    process.communicate()


def decompress_mprof(input, output, mprof_folder):
    outputdir = os.path.dirname(output)
    process = Popen([config.memusg_path, '7za', 'e', input, '-o' + outputdir, '-y',
                     "-mmt=1"], stdout=PIPE, stderr=PIPE)
    process.communicate()
    util.safe_move('memory_profile.csv', os.path.join(mprof_folder, os.path.splitext(os.path.basename(
        input))[0]+'-7zip-decompress.csv'))


if __name__ == "__main__":
    pass
