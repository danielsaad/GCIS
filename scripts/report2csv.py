import sys
import os
import re
import csv


rpn_report2csv = True
gcis_report2csv = True

'''
    Takes a log file produced by Navarro's repair
    and translate it to .csv
'''


def report2csv_rpn(report_path, csv_folder_path):
    files = [os.path.join(report_path, f) for f in os.listdir(
        report_path) if os.path.isfile(os.path.join(report_path, f))]
    files = [f for f in files if f.endswith('rpn.log')]
    for f in files:
        csv_file_path = os.path.join(
            csv_folder_path, os.path.basename(f)[:-8]+'-rpn-statistics.csv')
        print("Processing", os.path.basename(f))
        print(f, '->', csv_file_path+'\n')
        report2csv_rpn_helper(f, csv_file_path)


def report2csv_rpn_helper(report_path, csv_path):
    f = open(report_path, 'r')
    dict = {}
    for line in f:
        if(re.search('Original chars:', line)):
            dict['original_chars'] = line.split(' ')[-1].rstrip()
        elif(re.search('Number of rules: ', line)):
            dict['number_of_rules'] = line.split(' ')[-1].rstrip()
        elif(re.search('Final sequence length: ', line)):
            dict['final_sequence_length'] = line.split(' ')[-1].rstrip()
        elif(re.search('Compression ratio: ', line)):
            dict['compression_ratio'] = line.split(' ')[-1].rstrip('% \n')
    print('Found dictionary = ', dict)

    f.close()
    with open(csv_path, "w") as csv_file:
        field_names = ['original_chars', 'number_of_rules',
                       'final_sequence_length', 'compression_ratio']
        writer = csv.DictWriter(
            csv_file, fieldnames=field_names, delimiter=',')
        writer.writeheader()
        dict_rows = [dict]
        writer.writerow(dict_rows[0])


def report2csv(report_path, csv_folder_path):
    if(not os.path.isdir(report_path)):
        print(report_path, 'is not a directory')
        exit(0)
    if(not os.path.isdir(csv_folder_path)):
        print('Creating CSV folder at', csv_folder_path)
        os.makedirs(csv_folder_path, exist_ok=True)
    if(rpn_report2csv):
        print("Translating Navarro's Repair reports into csv files")
        report2csv_rpn(report_path, csv_folder_path)
    if(gcis_report2csv):
        print("Translating GCIS' reports into csv files")
        report2csv_gcis(report_path, csv_folder_path)


def report2csv_gcis(report_path, csv_folder_path):
    files = [os.path.join(report_path, f) for f in os.listdir(
        report_path) if os.path.isfile(os.path.join(report_path, f))]
    files = [f for f in files if f.endswith('gcis.log')]
    for f in files:
        csv_file_path = os.path.join(csv_folder_path, os.path.basename(f))[
            :-9]+'-gcis-statistics.csv'
        print("Processing", os.path.basename(f))
        print(f, '->', csv_file_path+'\n')
        report2csv_gcis_helper(f, csv_file_path)


def report2csv_gcis_helper(report_path, csv_path):
    level = []
    alphabet_size = []
    string_size = []
    number_of_rules = []
    rules_total_length = []
    avg_rule_length = []
    avg_lcp = []
    avg_rule_suffix_length = []
    dict_size = []
    lcp_size = []
    total_rule_suffix_length = []
    rule_suffix_width = []
    tail_length = []
    tail_width = []
    rle_potential = []
    avg_rle = []
    reduced_string_length = []
    reduced_string_width = []
    f = open(report_path, 'r')
    for line in f:
        if(re.search("^Level.*", line)):
            level.append(line.split(" ")[-2])
        elif(re.search("^String Size.*", line)):
            string_size.append(line.split(" ")[-2])
        elif(re.search("^Alphabet Size.*", line)):
            alphabet_size.append(line.split(" ")[-2])
        elif(re.search("^Number of Rules.*", line)):
            number_of_rules.append(line.split(" ")[-2])
        elif(re.search("^Total rules length*", line)):
            rules_total_length.append(line.split(" ")[-2])
        elif(re.search("^Average Rule Length.*", line)):
            avg_rule_length.append(line.split(" ")[-2])
        elif(re.search("^Average LCP.*", line)):
            avg_lcp.append(line.split(" ")[-2])
        elif(re.search("^Average Rule Suffix Length.*", line)):
            avg_rule_suffix_length.append(line.split(" ")[-2])
        elif(re.search("^Dictionary Level Size.*", line)):
            dict_size.append(line.split(" ")[-2])
        elif(re.search("^LCP Size.*", line)):
            lcp_size.append(line.split(" ")[-2])
        elif(re.search("^Rule Suffix Length.*", line)):
            total_rule_suffix_length.append(line.split(" ")[-2])
        elif(re.search("^Rule Suffix Width.*", line)):
            rule_suffix_width.append(line.split(" ")[-2])
        elif(re.search("^Tail Length.*", line)):
            tail_length.append(line.split(" ")[-2])
        elif(re.search("^Tail Width.*", line)):
            tail_width.append(line.split(" ")[-2])
        elif(re.search("^Run Length Potential.*", line)):
            rle_potential.append(line.split(" ")[-2])
        elif(re.search("^Avg Run Length per Rule Suffix.*", line)):
            avg_rle.append(line.split(" ")[-2])
        elif(re.search("^Reduced String Length.*", line)):
            reduced_string_length.append(line.split(" ")[-2])
        elif(re.search("^Reduced String Width.*", line)):
            reduced_string_width.append(line.split(" ")[-2])

    f.close()
    # print(level,string_size,alphabet_size,number_of_rules,avg_rule_length,
    #     avg_lcp,avg_rule_suffix_length,dict_size,lcp_size,total_rule_suffix_length,
    #     rule_suffix_width,tail_length,tail_width,rle_potential,avg_rle,
    #     reduced_string_length,reduced_string_width,sep='\n')

    csvfile = open(csv_path, "w")
    writer = csv.writer(csvfile, delimiter=',')
    experiment_name = os.path.split(report_path)[1][:-4]
    writer.writerow([experiment_name])
    header = ['Level', 'String Size (symbols)', 'Alphabet Size', 'Number of Rules', 'Rules Total Length',
              'Avg Rule Length',
              'Avg LCP (symbols)', 'LCP Size (bits)', 'Avg Rule Suffix Length (symbols)',
              'Rule Suffix Length (total symbols)', 'Rule Suffix Width (bits)',
              'RLE potential (symbols)', 'Avg RLE per Rule Suffix (symbols)',
              'Tail Length (symbols)', 'Tail Width (bits)', 'Dictionary Size (bytes)']

    writer.writerow(header)
    for i, v in enumerate(level):
        writer.writerow([level[i], string_size[i],
                         alphabet_size[i], number_of_rules[i], rules_total_length[i], avg_rule_length[i],
                         avg_lcp[i], lcp_size[i],
                         avg_rule_suffix_length[i], total_rule_suffix_length[i],
                         rule_suffix_width[i], rle_potential[i],
                         avg_rle[i], tail_length[i],
                         tail_width[i], dict_size[i]])
    writer.writerow(['Reduced String Length', 'Reduced String Width'])
    writer.writerow([reduced_string_length[0], reduced_string_width[0]])
    csvfile.close()


''' argv[1] = path to the folder where report files are stored '''
''' argv[2] = path to the folder where the csv files must be stored '''

if __name__ == "__main__":
    if(len(sys.argv) < 3):
        print('Usage: report2csv report_folder csv_statistics_folder')
        exit(0)
    report2csv(sys.argv[1], sys.argv[2])
