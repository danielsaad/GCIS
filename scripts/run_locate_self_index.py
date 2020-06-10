import os
import subprocess
from sys import argv
import generate_locate_input
import gcis
import csv

'''
    Generate queries from input folder and stores them in queries_folder
'''


def generate_queries(input_folder, queries_folder):
    input_files = [os.path.join(input_folder, f)
                   for f in os.listdir(input_folder)]
    query_filenames = [os.path.join(queries_folder, f)
                       for f in os.listdir(input_folder)]

    for (inf, qry) in zip(input_files, query_filenames):
        generate_locate_input.generate_extract_input(inf, qry)


'''
    Computes the index for every file in input_folder and store the
    respective index in the indices_folder
'''


def compute_indices(input_folder, indices_folder):
    print('Computing indices of ', input_folder)
    input_files = [os.path.join(input_folder, f)
                   for f in os.listdir(input_folder)]
    index_files = [os.path.join(indices_folder, f+'.gcisidx')
                   for f in os.listdir(input_folder)]

    for (inf, ouf) in zip(input_files, index_files):
        print('Building index file for', os.path.basename(inf))
        gcis.gcis_build_index(inf, ouf)


'''
    Check if the folders exist and
    take proper actions.
'''


def check_folders(input_folder, indices_folder, queries_folder, output_folder):
    if(not os.path.isdir(input_folder)):
        # The input folder does not exists. Aborting the program.
        print('Error: ', input_folder, 'does not exists')
        exit(0)
    if(not os.path.isdir(indices_folder)):
        ''' The indices were not computed. Compute all indices from the
            input folder '''
        print(indices_folder, 'does not exists, creating.')
        os.makedirs(indices_folder)
        compute_indices(input_folder, indices_folder)
    if(not os.path.isdir(output_folder)):
        # Create the output folder, were the .csv files shall be stored.
        print(output_folder, 'does not exists, creating')
        os.makedirs(output_folder, exist_ok=True)
    if(not os.path.isdir(queries_folder)):
        # Check if the queries folder exists and creates it if necessary.
        print(queries_folder, 'does not exists, creating')
        os.makedirs(queries_folder)
        generate_queries(input_folder, queries_folder)

    for inf in [os.path.join(input_folder, f) for f in os.listdir(input_folder)]:
        basename = os.path.basename(inf)
        ouf = os.path.join(indices_folder, basename+'.gcisidx')
        if(not os.path.isfile(ouf)):
            print('Missing index for ', basename, 'in', output_folder,'building it.')
            gcis.gcis_build_index(inf, ouf)


def parse_file(statistics_ouptut):
    pattern_len = 0
    occ_total = 0
    total_locate_time = 0
    mean_locate_time_per_pattern = 0
    mean_locate_time_per_occ = 0

    text = statistics_ouptut.decode()
    print(text)
    text = text.split('\n')
    pattern_len = text[2].split(' ')[-1]
    occ_total = text[3].split(' ')[-1]
    total_locate_time = text[4].split(' ')[-2]
    mean_locate_time_per_pattern = text[5].split(' ')[-4]
    mean_locate_time_per_occ = text[6].split(' ')[-2]

    return (pattern_len,occ_total,total_locate_time,mean_locate_time_per_pattern,mean_locate_time_per_occ)


def run_extract_self_index(input_folder, indices_folder, queries_folder, output_folder):
    check_folders(input_folder, indices_folder, queries_folder, output_folder)

    input_files = [os.path.join(input_folder, f) for f in os.listdir(
        input_folder) if os.path.isfile(os.path.join(input_folder, f))]

    header = ['pattern length', 'number of occ',
              'total time (ms)', 'mean time per pattern (us)', 'mean time per occ (us)']
    for inf in input_files:
        basename = os.path.basename(inf)
        print('processing', basename)
        ouf = os.path.join(output_folder, basename+'-self_index_locate.csv')
        idxf = os.path.join(indices_folder, basename+'.gcisidx')
        if(not os.path.isfile(idxf)):
            print(idxf,' was not produced. Error. Skipping.')
            continue
        queries_files = [os.path.join(queries_folder, f) for f in os.listdir(
            queries_folder) if f.startswith(basename)]
        if(not queries_files):
            print('There are no query file for ', basename)
            continue

        info = []
        for qf in queries_files:
            out, _ = gcis.gcis_self_index_locate(inf, idxf, qf)
            info.append(parse_file(out))
            print(info)
            # check if profile ouput exists and renames it to the result folder
            if(os.path.isfile('gperf.out')):
                new_gperf_location = os.path.join(
                    output_folder, basename+'-self_index_locate-'+info[-1][0]+'-gperf.out')
                os.rename('gperf.out', new_gperf_location)
                # convert to cachegrind format
                with open(os.path.join(output_folder, basename+'-self_index_locate-'+info[-1][0]+'-gperf.callgrind'), 'w') as callgrind_file:
                    subprocess.run(['pprof', '--callgrind', '../bin/gcis-ef',
                                    new_gperf_location], stdout=callgrind_file)

        info = sorted(info, key=lambda x: int(x[0]))
        with open(ouf, 'w') as csv_file:
            writer = csv.writer(csv_file)
            writer.writerow(header)
            for i in info:
                writer.writerow(i)
        
        # check if profile ouput exists and renames it to the result folder
        if(os.path.isfile('gmon.out')):
            os.rename('gmon.out',os.path.join(output_folder,basename+'-self_index_locate-gmon.out'))



'''
    Prints to the screen the correct way to execute the program
'''


def print_usage(argv):
    print('Usage: ',
          argv[0], '<input folder> <indices folder> <queries_folder> <output folder>')


'''
    argv[1] = Experiments Folder
    argv[2] = Index Folder
    argv[3] = queries_folders
    argv[4] = Output Folder
'''

if __name__ == "__main__":
    if(len(argv) < 5):
        print_usage(argv)
        exit(0)
    run_extract_self_index(argv[1], argv[2], argv[3], argv[4])
