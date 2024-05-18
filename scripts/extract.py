import os
import gcis
import bigrepair
import repair
import tomohiro_slp
import operator
import generate_extract_input
import csv
import util
import statistics

NUMBER_OF_ITERATIONS = 10

tomohiro_variants = ['PlainSlp_32Fblc', 'PlainSlp_FblcFblc',
                     'PlainSlp_IblcFblc', 'PoSlp_Iblc', 'PoSlp_Sd']


def build_extractors(corpora_folder_p, extractor_folder_p, machine):

    if(not os.path.isdir(extractor_folder_p)):
        print('Creating directory', extractor_folder_p)
        os.makedirs(extractor_folder_p, exist_ok=False)
    input_files_p = [os.path.join(corpora_folder_p, f)
                     for f in os.listdir(corpora_folder_p)]

    # Create Indices
    for f in input_files_p:
        gcis_ouf = [os.path.join(extractor_folder_p, os.path.basename(
            f)+ext) for ext in ['.gcis', '.gcis64']]
        repair_ouf = [os.path.join(
            extractor_folder_p, os.path.basename(f)+ext) for ext in ['.C', '.R']]
        already_built_gcis = [os.path.isfile(f) for f in gcis_ouf]
        already_built_repair = [os.path.isfile(f) for f in repair_ouf]
        if(not all(already_built_repair)):
            if(machine != 'JAU'):
                repair.compress_repair_navarro(f, extractor_folder_p)
            else:
                bigrepair.compress(f, extractor_folder_p)
        else:
            print('Skipping, RePair grammar for ', f, ' already build')
        if(not any(already_built_gcis)):
            gcis.compress(f, os.path.join(extractor_folder_p,
                                          os.path.basename(f)+'.gcis'), '-ef')
        else:
            print('Skipping, GCIS extractor for ', f, ' already built')

    # Create SLP Encodings for RePair Grammars
    for f in input_files_p:
        repair_basename_p = os.path.join(
            extractor_folder_p, os.path.basename(f))
        repair_files_p = [os.path.join(
            extractor_folder_p, os.path.basename(f)+ext) for ext in ['.C', '.R']]
        # The Grammar Files are Present
        if(all([os.path.isfile(p) for p in repair_files_p])):
            for var in tomohiro_variants:
                slp_ouf = os.path.join(
                    extractor_folder_p, os.path.basename(f)+'.slp-'+var)
                if(os.path.isfile(slp_ouf)):
                    print('Skipping construction of SLP encoder ',
                          slp_ouf, 'already built')
                else:
                    print('Encoding SLP variant', var,
                          'for', repair_basename_p)
                    format_type = 'NavarroRepair' if machine != 'JAU' else 'Bigrepair'
                    tomohiro_slp.build_slp_encoding(
                        repair_basename_p, slp_ouf, format_type, var)
        else:
            print('Grammar file for ', f, ' does not exists')


def create_extract_input(corpora_folder_p, extractor_folder_p):
    # Create extract queries directory
    os.makedirs(os.path.join(extractor_folder_p,
                             'extract_queries'), exist_ok=True)
    input_files_p = [os.path.join(corpora_folder_p, f)
                     for f in os.listdir(corpora_folder_p)]
    for inf in input_files_p:
        generate_extract_input.generate_extract_input(inf, os.path.join(
            *[extractor_folder_p, 'extract_queries', os.path.basename(inf)]))


def extract(corpora_folder_p, extractor_folder_p):

    extract_queries_p = os.path.join(extractor_folder_p, 'extract_queries')
    extract_results_folder = os.path.join(
        extractor_folder_p, 'extract_results')
    os.makedirs(extract_results_folder, exist_ok=True)
    input_files = [os.path.join(corpora_folder_p, f)
                   for f in os.listdir(corpora_folder_p)]

    extract_data = {}
    for inf in input_files:
        csv_time_fp = open(os.path.join(
            extract_results_folder, os.path.basename(inf)+'-extract-time.csv'), 'w')
        csv_space_fp = open(os.path.join(
            extract_results_folder, os.path.basename(inf)+'-compression-ratio.csv'), 'w')
        csv_time_writer = csv.writer(csv_time_fp)
        csv_space_writer = csv.writer(csv_space_fp)
        query_files_p = [os.path.join(extract_queries_p, f) for f in os.listdir(
            extract_queries_p) if f.startswith(os.path.basename(inf))]
        extract_data = {}
        space_data = {}
        for q_p in query_files_p:
            substring_size = int(os.path.splitext(
                os.path.basename(q_p))[1].split('_')[0][1:])
            extract_data[substring_size] = {
                k: 'ERROR' for k in ['GCIS']+tomohiro_variants}
            space_data = {k: 'ERROR' for k in ['GCIS'] + tomohiro_variants}

            print('Extracting ', os.path.basename(q_p))
            # Extract with Tomohiro's SLP
            for var in tomohiro_variants:
                slp_p = os.path.join(extractor_folder_p,
                                     os.path.basename(inf)+'.slp-'+var)
                if(os.path.isfile(slp_p)):
                    print(var)
                    time_list = []
                    for _ in range(NUMBER_OF_ITERATIONS):
                        p = tomohiro_slp.random_extract(slp_p, var, q_p)
                        t = tomohiro_slp.parse_extract_time(p)
                        time_list.append(t)
                    print(time_list)
                    extract_data[substring_size][var] = 'ERROR' if 'ERROR' in time_list else str(
                        statistics.median(sorted([float(x) for x in time_list])))
                    space_data[var] = util.calc_compression_ratio(slp_p, inf)
                else:
                    print('SLP file ', slp_p, 'is not present')
            # Extract with GCIS
            gcis_p = [os.path.join(extractor_folder_p, x) for x in os.listdir(
                extractor_folder_p) if x.startswith(os.path.basename(inf)+'.gcis')][0]

            if(os.path.isfile(gcis_p)):
                time_list = []
                print('GCIS')
                for _ in range(NUMBER_OF_ITERATIONS):
                    p = gcis.gcis_extract(gcis_p, q_p)
                    t = gcis.parse_extract(p)
                    time_list.append(t)
                space_data['GCIS'] = util.calc_compression_ratio(gcis_p, inf)
            else:
                print('GCIS file', gcis_p, 'is not present')
            print(time_list)
            extract_data[substring_size]['GCIS'] = 'ERROR' if 'ERROR' in time_list else str(
                statistics.median(sorted([float(x) for x in time_list])))
            print(extract_data)
        csv_time_writer.writerow(['n']+['GCIS']+tomohiro_variants)
        csv_space_writer.writerow(['Extractor', 'Compression Ratio'])
        for k in sorted(extract_data.keys()):
            row = [str(k)] + [extract_data[k][variant]
                              for variant in ['GCIS']+tomohiro_variants]
            csv_time_writer.writerow(row)

        for k in ['GCIS'] + tomohiro_variants:
            csv_space_writer.writerow([k, space_data[k]])

        csv_time_fp.close()
        csv_space_fp.close()
        print(extract_data)
        print(space_data)
