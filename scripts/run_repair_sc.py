import os
import subprocess
import sys

# rp_bin = '/home/danielsaad/git/libcds2/libcds2/bin/repair_extract_sc'
# corpus_folder = '/home/danielsaad/corpora/pizza-chili-repetitive-bkp'
# output_folder = os.path.join(corpus_folder,'repairsc')


def run(rp_bin, corpus_folder, output_folder, c_sample_rate, delta_sample_rate, ss_rate):
    processes = []
    os.makedirs(output_folder, exist_ok=True)

    filename = [f for f in os.listdir(corpus_folder) if os.path.isfile(os.path.join(corpus_folder,f))]
    for f in filename:
        fp = os.path.join(corpus_folder, f)
        of = os.path.join(
            output_folder, f+'.rpsc-{}-{}-{}'.format(c_sample_rate, delta_sample_rate, ss_rate))
        if (os.path.isfile(of)):
            print(of, 'already present: skipping')
            continue
        command = [rp_bin, '-c', fp, of,
                   c_sample_rate, delta_sample_rate, ss_rate]
        print('Command = ', command)
        p = subprocess.Popen(command)
        processes.append(p)

    for p in processes:
        p.wait()


if __name__ == "__main__":
    # repair binary
    rp_bin = sys.argv[1]
    # corpus folder
    corpus_folder = sys.argv[2]
    # output folder
    output_folder = sys.argv[3]
    # c sample rate parameter
    c_sample_rate = sys.argv[4]
    # delta sampling parameter
    delta_sample_rate = sys.argv[5]
    # superblock sample parameter
    ss_rate = sys.argv[6]
    run(rp_bin, corpus_folder, output_folder,
        c_sample_rate, delta_sample_rate, ss_rate)
