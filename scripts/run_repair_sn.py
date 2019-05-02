import os
import subprocess
import sys


'''
    Creates a RepairSN dictionaries from every input folder
'''


def run(rp_bin, corpus_folder, output_folder, c_sample_rate, delta_sample_rate, ss_rate):
    processes = []
    os.makedirs(output_folder, exist_ok=True)

    filename = [f for f in os.listdir(corpus_folder) if os.path.isfile(
        os.path.join(corpus_folder, f))]
    for f in filename:
        fp = os.path.join(corpus_folder, f)
        of = os.path.join(
            output_folder, f+'.rpsn-{}-{}-{}'.format(c_sample_rate, delta_sample_rate, ss_rate))
        if(os.path.isfile(of)):
            print(of, 'already present: skipping')
            continue
        command = [rp_bin, '-c', fp, of,
                   c_sample_rate, delta_sample_rate, ss_rate]
        print('Command = ', command)
        p = subprocess.Popen(command)
        processes.append(p)

    for p in processes:
        p.wait()


"""
    Using rp_bin, RePair compress every input file from input_list into output files from
    output_list using the parameters c_sample_rate, delta_sample_rate and
    ss_rate
"""


def compress_repair_s_parallel(rp_bin, input_list, output_list, repair_parameters):

    # List of processes
    processes = []
    for (inf, ouf) in zip(input_list, output_list):

        # Check if Dictionary already exists
        if(os.path.isfile(ouf)):
            # If exists, skip
            print(os.path.basename(ouf), 'already present: skipping')
            continue

        # Run RepairSN asynchronously to compute every dictionary
        command = [rp_bin, '-c', inf, ouf,
                   str(repair_parameters[0]), str(repair_parameters[1]), str(repair_parameters[2])]
        print('Command = ', command)
        p = subprocess.Popen(command)
        processes.append(p)

    # Wait for the processses to complete
    for p in processes:
        p.wait()


"""
    Using rp_bin, RePair compress  input file into output file from
    output_list using the parameters c_sample_rate, delta_sample_rate and
    ss_rate
"""


def compress_repair_s(rp_bin, inf, ouf, repair_parameters):

    # Check if Dictionary already exists
    if(os.path.isfile(ouf)):
        # If exists, skip
        print(os.path.basename(ouf), 'already present: skipping')
    else:
        # Run RepairS
        p = repair_parameters
        command = [rp_bin, '-c', inf, ouf,
                   str(p[0]), str(p[1]), str(p[2])]
        print('Command = ', command)
        p = subprocess.Popen(command)
        p.communicate()


def extract_repair_s(rp_bin, inf, query_file,repair_parameters):
    # Run RepairS extraction
    p = repair_parameters
    command = [rp_bin, '-e', inf, query_file,
               str(p[0]), str(p[1]), str(p[2])]
    print('Command = ', command)
    p = subprocess.Popen(command,stdout=subprocess.PIPE)
    return p.communicate()


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
