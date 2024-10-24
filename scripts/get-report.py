import os
import subprocess
import repair
import gcis
import sys

gcis_report = True
repair_navarro_report = True

'''
    Run GCIS in order to get report.log statistics.
    Receives an experiment folder and stores the resulting reports
    in an output folder.
'''
def run_gcis(experiments_folder,output_folder):
    files = [os.path.join(experiments_folder,f) for f in os.listdir(experiments_folder)]
    for f in files:
        tmp_file = os.path.basename(f)+'.gcis'
        gcis.compress_gc_is_statistics(f,tmp_file)
        print('Removing index file')
        os.remove(tmp_file)
        print('Moving log file to output folder')
        report_file = os.path.join(output_folder,os.path.basename(f)+'.gcis.log')
        os.rename('report.log',report_file)


def run_repair_navarro(experiments_folder,output_folder):
    files = [os.path.join(experiments_folder,f) for f in os.listdir(experiments_folder)]
    for f in files:
        _,err = repair.compress_repair_navarro(f)
        print('Removing index file')
        os.remove(f+'.R')
        os.remove(f+'.C')
        print('Moving log file to output folder')
        report_file = os.path.join(output_folder,os.path.basename(f)+'.rpn.log') 
        with open(report_file,'w') as ouf:
            print(err.decode('utf-8'),file=ouf)


def get_report(experiments_folder,report_folder):
    if(not os.path.isdir(experiments_folder)):
        print(experiments_folder,'does not exists')
        return
    if(not os.path.isdir(report_folder)):
        print('Creating',report_folder)
        os.makedirs(report_folder,exist_ok=False)

    if(gcis_report):
        print('Running GCIS on',experiments_folder)
        run_gcis(experiments_folder,report_folder)
    if(repair_navarro_report):
        print('Running Repair-Navarro on',experiments_folder)
        run_repair_navarro(experiments_folder,report_folder)

'''
    argv[1] = experiments_folder
    argv[2] = report folder
'''
if __name__ == "__main__":
    if(len(sys.argv)<3):
        print('Usage: get-report <experiments_folder> <report_folder>')
    get_report(sys.argv[1],sys.argv[2])
    pass