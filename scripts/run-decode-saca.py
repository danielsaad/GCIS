import os
import sys
import subprocess
import time
import csv

class experiment_data:
    def __init__(self):
        self.experiment_data = ""
        self.experiment_time = []
        self.experiment_input = []

def run_sais_nong(input_folder_path,output_folder_path):
    e = experiment_data()
    e.experiment_name = "SAIS-NONG"
    input_filename = os.listdir(input_folder_path)
    for f in input_filename:
        input_file = os.path.join(input_folder_path,f)
        output_file = os.path.join(output_folder_path,f+'.sanong')

        print("SACA Nong: ",f)
        start_time = time.perf_counter() # Start time
        process = subprocess.run(['../bin/sais-nong', '-c', input_file,output_file],
        stdout=subprocess.PIPE,stderr=subprocess.PIPE)
        end_time = time.perf_counter() # End time
        print("Time Nong =","{:.2f}".format(end_time-start_time))
        e.experiment_input.append(f)
        e.experiment_time.append("{:.2f}".format(end_time-start_time))

    return e


def run_sais_yuta(input_folder_path,output_folder_path):
    e = experiment_data()
    e.experiment_name = "SAIS-YUTA"
    input_filename = os.listdir(input_folder_path)
    for f in input_filename:
        input_file = os.path.join(input_folder_path,f)

        output_file = os.path.join(output_folder_path,f+'.sayuta')
        print("SACA Yuta: ",f)
        start_time = time.perf_counter() # Start time
        process = subprocess.run(['../bin/sais-yuta', '-c', input_file,output_file],
        stdout=subprocess.PIPE,stderr=subprocess.PIPE)
        end_time = time.perf_counter() # End time
        print("Time Yuta =","{:.2f}".format(end_time-start_time))
        e.experiment_input.append(f)
        e.experiment_time.append("{:.2f}".format(end_time-start_time))

    return e



def run_decode_saca(input_folder_path,output_folder_path):

    e = experiment_data()
    e.experiment_name = "GCIS-SACA"
    input_filename = os.listdir(input_folder_path)
    for f in input_filename:
        input_file = os.path.join(input_folder_path,f)
        compressed_file = os.path.join(output_folder_path,f+'.gcis')
        output_file = os.path.join(output_folder_path,f+'.sa')

        print("GCIS Compressing: ",f)
        process = subprocess.run(['../bin/gc-is-codec','-c',input_file,compressed_file],
        stdout=subprocess.PIPE,stderr=subprocess.PIPE)

        print("GCIS Decoding + SAIS: ",f)
        start_time = time.perf_counter() # Start time
        process = subprocess.run(['../bin/gc-is-codec', '-s', compressed_file,output_file],
        stdout=subprocess.PIPE,stderr=subprocess.PIPE)
        end_time = time.perf_counter() # End time
        print("Time GCIS =","{:.2f}".format(end_time-start_time))

        e.experiment_input.append(f)
        e.experiment_time.append("{:.2f}".format(end_time-start_time))

    return e

#argv[1] contains the path to the experiments folder
#argv[2] contains the path to the output folder

def run_saca(input_folder_path,output_folder_path):
    
    decode_saca_data = run_decode_saca(input_folder_path,output_folder_path)
    saca_nong_data = run_sais_nong(input_folder_path,output_folder_path)
    saca_yuta_data = run_sais_yuta(input_folder_path,output_folder_path)
 
    with open(os.path.join(output_folder_path,'results-saca.csv'),'w') as csvfile:
        writer = csv.writer(csvfile, delimiter=',')
        header = ['Input File',decode_saca_data.experiment_name,saca_nong_data.experiment_name,saca_yuta_data.experiment_name]
        writer.writerow(header)
        time = list(zip(decode_saca_data.experiment_time,saca_nong_data.experiment_time,saca_yuta_data.experiment_time))
        for i,t in enumerate(time):
            row = [decode_saca_data.experiment_input[i],t[0],t[1],t[2]]
            writer.writerow(row)



if __name__ == "__main__":
    run_saca(sys.argv[1],sys.argv[2])
