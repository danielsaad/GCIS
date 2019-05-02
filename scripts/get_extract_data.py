import os
import csv
import subprocess
import itertools
import gcis
import run_repair_sc
import run_repair_sn
import generate_extract_input
import sys

# data_path = '/home/danielsaad/corpora/pizza-chili-repetitive-bkp'
# repair_sc_folder = os.path.join(data_path, 'repairsc')
# repair_sn_folder = os.path.join(data_path, 'repairsn')
# gcis_folder = os.path.join(data_path, 'gcis')
# query_folder = os.path.join(*[data_path, 'statistics', 'query'])

# TODO change hardcoded binaries from repair
repair_sc_binary = '/home/danielsaad/git/libcds2/bin/repair_extract_sc'
repair_sn_binary = '/home/danielsaad/git/libcds2/bin/repair_extract_sn'


sampling_rate_rp_sc = [6, 7, 8, 9, 10]
sampling_rate_rp_sn = [10, 11, 12, 13, 14]
rule_sampling = [0, 1, 2, 3, 4]
superblock_sampling = [0, 5, 8]

rp_sc_parameters = [
    (x, y, z) for x in sampling_rate_rp_sc for y in rule_sampling for z in superblock_sampling]
rp_sn_parameters = [
    (x, y, z) for x in sampling_rate_rp_sn for y in rule_sampling for z in superblock_sampling]

# # print(rp_sc_parameters)
# # print(rp_sn_parameters)


def extract(data_folder, results_folder):

    # Input files
    file_path_list = sorted([os.path.join(data_folder, f) for f in os.listdir(data_folder) if
                             os.path.isfile(os.path.join(data_folder, f))])

    # Size in bytes for each input file
    file_size_list = [os.path.getsize(f) for f in file_path_list]

    # Dictionary folder
    dictionary_folder = os.path.join(results_folder, 'dictionaries')

    # Query folder
    query_folder = os.path.join(results_folder, 'queries')

    # Bits per symbol for each variation
    bps_gcis = []
    bps_rp_sc = []
    bps_rp_sn = []

    # Time for each variation
    time_rp_sc = []
    time_rp_sn = []
    time_gcis = []

    print('Collecting GCIS extract data')
    # Gather GCIS data
    for f, sz in zip(file_path_list, file_size_list):

        filename = os.path.basename(f)
        print(filename, 'size = ', os.path.getsize(f))

        dictionary = os.path.join(dictionary_folder, filename+'.gcis')

        # We are only regarding the extraction of a single symbol
        query_file = os.path.join(query_folder, filename+'.1_query')

        # Get GCIS dictionary size
        file_size_gcis = os.path.getsize(dictionary)

        # Append to the list the size in bits of the gcis file
        bps_gcis.append(8*float(file_size_gcis)/sz)
        out, err = gcis.gcis_extract(dictionary, query_file)
        total_time = out.split(b'\n')[-2].split(b' ')[-1]
        time_gcis.append(float(total_time))

        print('Time GCIS (seconds) = ', total_time.decode())

    print('Collecting RPSC extract data')
    # Gather RPSC data
    for parameter in rp_sc_parameters:
        time = []
        bps = []
        for f, sz in zip(file_path_list, file_size_list):
            filename = os.path.basename(f)
            print('Processing', filename)
            dictionary = os.path.join(
                dictionary_folder, filename+'.rpsc-{}-{}-{}'.format(parameter[0], parameter[1], parameter[2]))

            # We are only regarding the extraction of a single symbol
            query_file = os.path.join(query_folder, filename+'.1_query')

            # Get RepairSC dictionary file
            file_size_rp_sc = os.path.getsize(os.path.join(dictionary))

            # Append RepairSC bps
            bps.append(8*float(file_size_rp_sc)/sz)

            out, err = run_repair_sn.extract_repair_s(
                repair_sc_binary, dictionary, query_file, parameter)
            out = out.decode("utf-8")
            total_time = out.split('\n')[-2].split(' ')[-1]
            print('Time Repair SC (seconds) = ', total_time)
            time.append(float(total_time))

        time_rp_sc.append(time)
        bps_rp_sc.append(bps)

    print('Collecting RPSN extract data')
    # Gather RPSN data
    for parameter in rp_sn_parameters:
        time = []
        bps = []
        for f, sz in zip(file_path_list, file_size_list):
            filename = os.path.basename(f)
            print('Processing', filename)
            dictionary = os.path.join(
                dictionary_folder, filename+'.rpsn-{}-{}-{}'.format(parameter[0], parameter[1], parameter[2]))

            # We are only regarding the extraction of a single symbol
            query_file = os.path.join(query_folder, filename+'.1_query')

            # Get RepairSN dictionary file
            file_size_rp_sn = os.path.getsize(os.path.join(dictionary))

            # Append RepairSC bps
            bps.append(8*float(file_size_rp_sc)/sz)

            out, err = run_repair_sn.extract_repair_s(
                repair_sn_binary, dictionary, query_file, parameter)
            out = out.decode("utf-8")
            total_time = out.split('\n')[-2].split(' ')[-1]
            print('Time Repair SN (seconds) = ', total_time)
            time.append(float(total_time))

        time_rp_sn.append(time)
        bps_rp_sn.append(bps)

    # Store results in csv file
    csv_file = os.path.join(results_folder, 'results_extract.csv')
    with open(csv_file, 'w') as of:
        writer = csv.writer(of)
        header = ['Filename', 'File Size(MB)', 'GCIS bps', 'Access GCIS (us)']
        for parameter in rp_sc_parameters:
            header.append(
                'RePair SC bps {}-{}-{}'.format(parameter[0], parameter[1], parameter[2]))
            header.append(
                'Access RePair SC (us) {}-{}-{}'.format(parameter[0], parameter[1], parameter[2]))
        for parameter in rp_sn_parameters:
            header.append(
                'RePair SN bps {}-{}-{}'.format(parameter[0], parameter[1], parameter[2]))
            header.append(
                'Access RePair SN (us) {}-{}-{}'.format(parameter[0], parameter[1], parameter[2]))

        # Write header
        writer.writerow(header)

        for i in range(len(file_path_list)):
            line = [os.path.basename(file_path_list[i]), "{}".format(
                file_size_list[i]/1000000), "{}".format(bps_gcis[i]), '{}'.format(time_gcis[i]*100)]
            for j in range(len(rp_sc_parameters)):
                line += ['{}'.format(bps_rp_sc[j][i]), '{}'.format(
                    time_rp_sc[j][i]*100)]
            for j in range(len(rp_sn_parameters)):
                line += ['{}'.format(bps_rp_sn[j][i]), '{}'.format(
                    time_rp_sn[j][i]*100)]
            writer.writerow(line)


'''
    Computes GCIS, RepairSC and RepairSN dictionaries files in order
    to perform extraction
'''


def compute_dictionaries(data_folder, dictionary_folder):

    file_path_list = sorted([os.path.join(data_folder, f) for f in os.listdir(data_folder) if
                             os.path.isfile(os.path.join(data_folder, f))])

    # Generate output file lists for every compressor
    gcis_ouf = [os.path.join(dictionary_folder, os.path.basename(
        f)+'.gcis') for f in file_path_list]

    gcis.compress_gc_is_parallel(file_path_list, gcis_ouf)

    # For each RepairSN parameter, compress all files of the corpus
    for p in rp_sn_parameters:
        for f in file_path_list:
            rp_ouf = os.path.join(dictionary_folder, os.path.basename(
                f)+'.rpsn-{}-{}-{}'.format(p[0], p[1], p[2]))
            run_repair_sn.compress_repair_s(repair_sn_binary, f, rp_ouf, p)

    # For each RepairSC parameter, compress all files of the corpus
    for p in rp_sc_parameters:
        for f in file_path_list:
            rp_ouf = os.path.join(dictionary_folder, os.path.basename(
                f)+'.rpsc-{}-{}-{}'.format(p[0], p[1], p[2]))
            run_repair_sn.compress_repair_s(repair_sc_binary, f, rp_ouf, p)


def get_extract_data(data_folder, results_folder):
    print('Creating Dictionary folder')
    dictionary_folder = os.path.join(results_folder, 'dictionaries')
    os.makedirs(dictionary_folder, exist_ok=True)

    # Compute Dictionaries
    compute_dictionaries(data_folder, dictionary_folder)

    # Extract Symbols
    extract(data_folder, results_folder)


'''
    Generates query files for each input folder
    and places them into queries folder
'''


def produce_query_files(data_folder, results_folder):
    file_path_list = sorted([os.path.join(data_folder, f) for f in os.listdir(data_folder) if
                             os.path.isfile(os.path.join(data_folder, f))])
    query_folder = os.path.join(results_folder, 'queries')
    print('Creating queries folder')
    os.makedirs(query_folder, exist_ok=True)

    # For every input
    for f in file_path_list:
        # Get path for output file
        ouf = os.path.join(query_folder, os.path.basename(f))
        # Produce queries in output file
        generate_extract_input.generate_extract_input(f, ouf)


'''
    argv[1] = input_folder
    argv[2] = results folder
'''

if __name__ == "__main__":

    print('Creating Output directory', sys.argv[2])
    os.makedirs(sys.argv[2], exist_ok=True)
    produce_query_files(sys.argv[1], sys.argv[2])
    get_extract_data(sys.argv[1], sys.argv[2])
