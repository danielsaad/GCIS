import os
import csv
import subprocess
import itertools
import run_repair_sc
import run_repair_sn


data_path = '/home/danielsaad/corpora/pizza-chili-repetitive-bkp'
repair_sc_folder = os.path.join(data_path, 'repairsc')
repair_sn_folder = os.path.join(data_path, 'repairsn')
gcis_folder = os.path.join(data_path, 'gcis')
query_folder = os.path.join(*[data_path, 'statistics', 'query'])
repair_sc_binary = '/home/danielsaad/git/libcds2/libcds2/bin/repair_extract_sc'
repair_sn_binary = '/home/danielsaad/git/libcds2/libcds2/bin/repair_extract_sn'


file_path_list = sorted([os.path.join(data_path, f) for f in os.listdir(data_path) if
                         os.path.isfile(os.path.join(data_path, f))])
file_size_list = [os.path.getsize(f) for f in file_path_list]
bps_gcis = []
bps_rp_sc = []
bps_rp_sn = []
time_rp_sc = []
time_rp_sn = []
time_gcis = []

sampling_rate_rp_sc = [6, 7, 8, 9, 10]
sampling_rate_rp_sn = [10, 11, 12, 13, 14]
rule_sampling = [0, 1, 2, 3, 4]
superblock_sampling = [0, 5, 8]

rp_sc_parameters = [
    (x, y, z) for x in sampling_rate_rp_sc for y in rule_sampling for z in superblock_sampling]
rp_sn_parameters = [
    (x, y, z) for x in sampling_rate_rp_sn for y in rule_sampling for z in superblock_sampling]

# print(rp_sc_parameters)
# print(rp_sn_parameters)

# TODO: Run GCIS

print('Collecting GCIS extract data')

# Gather GCIS data
for f, sz in zip(file_path_list, file_size_list):
    filename = os.path.basename(f)
    print(filename, 'size = ', os.path.getsize(f))

    file_size_gcis = os.path.getsize(
        os.path.join(gcis_folder, filename+'.gcisef'))
    bps_gcis.append(8*float(file_size_gcis)/sz)
    command_gcis = ['/home/danielsaad/git/gcis/gcis-ef/bin/gc-is-codec',
                    '-e', os.path.join(gcis_folder, filename+'.gcisef'),
                    os.path.join(query_folder, filename+'.query')]
    print('GCIS command', command_gcis)
    p = subprocess.Popen(command_gcis, stdout=subprocess.PIPE)
    out, err = p.communicate()
    total_time = out.split(b'\n')[-2].split(b' ')[-1]
    time_gcis.append(float(total_time))
    print('Time GCIS (seconds) = ', total_time.decode())


print('Building RPSC dictionaries')

# Build RPSC and RPSN dictionaries
for parameter in rp_sc_parameters:
    # Run RepairSC to make the executables
    run_repair_sc.run(repair_sc_binary,
                      data_path, repair_sc_folder,
                      str(parameter[0]),
                      str(parameter[1]),
                      str(parameter[2]))

print('Building RPSN dictionaries')
for parameter in rp_sn_parameters:
    run_repair_sn.run(repair_sn_binary,
                      data_path, repair_sn_folder,
                      str(parameter[0]),
                      str(parameter[1]),
                      str(parameter[2]))


print('Collecting RPSC extract data')
# Gather RPSC data
for parameter in rp_sc_parameters:
    time = []
    bps = []
    for f, sz in zip(file_path_list, file_size_list):
        filename = os.path.basename(f)
        print('Processing', filename)
        file_size_rp_sc = os.path.getsize(os.path.join(
            repair_sc_folder, filename)+'.rpsc-{}-{}-{}'.format(parameter[0], parameter[1], parameter[2]))
        bps.append(8*float(file_size_rp_sc)/sz)
        repair_sc_dict_filepath = os.path.join(
            repair_sc_folder, filename + '.rpsc-{}-{}-{}'.format(parameter[0], parameter[1], parameter[2]))
        repair_sc_query_filepath = os.path.join(
            query_folder, filename)+'.query'
        command_repair_sc = [repair_sc_binary, '-e',
                             repair_sc_dict_filepath,
                             repair_sc_query_filepath,
                             str(parameter[0]),
                             str(parameter[1]),
                             str(parameter[2])]

        print('Repair SC command = ', ' '.join(command_repair_sc))
        p = subprocess.Popen(command_repair_sc, stdout=subprocess.PIPE)
        out, err = p.communicate()
        out = out.decode("utf-8")
        # print(out,err)
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
        file_size_rp_sn = os.path.getsize(os.path.join(
            repair_sn_folder, filename+'.rpsn-{}-{}-{}'.format(parameter[0], parameter[1], parameter[2])))
        bps.append(8*float(file_size_rp_sn)/sz)
        repair_sn_dict_filepath = os.path.join(
            repair_sn_folder, filename+'.rpsn-{}-{}-{}'.format(parameter[0], parameter[1], parameter[2]))
        repair_sn_query_filepath = os.path.join(
            query_folder, filename)+'.query'
        command_repair_sn = [repair_sn_binary,
                             '-e', repair_sn_dict_filepath,
                             repair_sn_query_filepath,
                             str(parameter[0]),
                             str(parameter[1]),
                             str(parameter[2])]

        print('Repair SN command = ', ' '.join(command_repair_sn))
        p = subprocess.Popen(command_repair_sn, stdout=subprocess.PIPE)
        out, err = p.communicate()
        out = out.decode("utf-8")
        # print(out,err)
        total_time = out.split('\n')[-2].split(' ')[-1]
        print('Time Repair SN (seconds) = ', total_time)
        time.append(float(total_time))
    time_rp_sn.append(time)
    bps_rp_sn.append(bps)

with open('results_extract.csv', 'w') as of:
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
            line.append('{}'.format(bps_rp_sc[j][i]), '{}'.format(
                time_rp_sc[j][i]*100))
        for j in range(len(rp_sn_parameters)):
            line.append('{}'.format(bps_rp_sn[j][i]), '{}'.format(
                time_rp_sn[j][i]*100))
        writer.writerow(line)
