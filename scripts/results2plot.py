import csv
import matplotlib.pyplot as plt
import numpy as np
import sys
import re
import os
import collections


def dtcr(expname,expsize,dt,cr,output_folder_path):
    row = 0
    for x in expname:

        fig = plt.figure()
        ax = fig.add_subplot(111)

        xgcis = dt[0][row]
        xrepair = dt[1][row]
        x7z = dt[2][row]
        xgz = dt[3][row]

        ygcis = cr[0][row]
        yrepair = cr[1][row]
        y7z = cr[2][row]
        ygz = cr[3][row]

        ax.plot(xgcis,ygcis,label='GCIS',color='green',marker='o')
        ax.plot(xrepair,yrepair,label='RePair',color='red',marker='^')
        ax.plot(x7z,y7z,label='7z',color='magenta',marker='*')
        ax.plot(xgz,ygz,label='gz',color='black',marker='s')

        size = round(float(expsize[row]), 2 )
        ax.set_title(expname[row] + ' ' + str(size) + ' (MB)')
        ax.set_xlabel('Decompression Time (s)')
        ax.set_ylabel('Compression Ratio (%)')
        ax.legend(loc='lower center', bbox_to_anchor=(0.5, -0.4))
        ax.grid(True)

#        fig.tight_layout()
#        plt.show()
        output_path = os.path.join(output_folder_path,'dtcr-' + expname[row]+'.png')
        plt.savefig(output_path,format='png',bbox_inches='tight',dpi=300)
        plt.close('all')
        row = row+1

def ctcr(expname,expsize,ct,cr,output_folder_path):
    row = 0
    for x in expname:

        fig = plt.figure()
        ax = fig.add_subplot(111)

        xgcis = ct[0][row]
        xrepair = ct[1][row]
        x7z = ct[2][row]
        xgz = ct[3][row]

        ygcis = cr[0][row]
        yrepair = cr[1][row]
        y7z = cr[2][row]
        ygz = cr[3][row]

        ax.plot(xgcis,ygcis,label='GCIS',color='green',marker='o')
        ax.plot(xrepair,yrepair,label='RePair',color='red',marker='^')
        ax.plot(x7z,y7z,label='7z',color='magenta',marker='*')
        ax.plot(xgz,ygz,label='gz',color='black',marker='s')

        size = round(float(expsize[row]), 2 )
        ax.set_title(expname[row] + ' ' + str(size) + ' (MB)')
#        ax.set_subtitle('Compression Time x Compression Ratio')
        ax.set_xlabel('Compression Time (s)')
        ax.set_ylabel('Compression Ratio (%)')
        ax.legend(loc='lower center', bbox_to_anchor=(0.5, -0.4))
        ax.grid(True)

#        fig.tight_layout()
#        plt.show()
        output_path = os.path.join(output_folder_path,'ctcr-' + expname[row]+'.png')
        plt.savefig(output_path,format='png',bbox_inches='tight',dpi=300)
        plt.close('all')
        row = row+1

def results2plot(results_folder_path, output_folder_path):

    os.makedirs(output_folder_path,exist_ok=True)

    ''' Lists for compression time
        Compression ratio and
        Decompression time '''
    exp_name = []
    exp_size = []
    ct_gcis = []
    cr_gcis = []
    dt_gcis = []
    ct_repair = []
    cr_repair = []
    dt_repair = []
    ct_7z = []
    cr_7z = []
    dt_7z = []
    ct_gz = []
    cr_gz = []
    dt_gz = []

    csv_input_path = os.path.join(results_folder_path,'results.csv')
    if(not os.path.isfile(csv_input_path)):
        print('results.csv do not exist')
        pass
    with open(csv_input_path) as csvfile:

        reader = csv.reader(csvfile)
        header = next(reader,None)

        ''' Get compression time, compression ratio and
            decompress time for every compressor in each experiment '''
        for row in reader:
#            print (row)
            exp_name.append(row[0])
            exp_size.append(row[1])

            cr_gcis.append(float(row[3]))
            ct_gcis.append(float(row[4]))
            dt_gcis.append(float(row[5]))

            cr_repair.append(float(row[8]))
            ct_repair.append(float(row[9]))
            dt_repair.append(float(row[10]))

            cr_7z.append(float(row[13]))
            ct_7z.append(float(row[14]))
            dt_7z.append(float(row[15]))

            cr_gz.append(float(row[18]))
            ct_gz.append(float(row[19]))
            dt_gz.append(float(row[20]))

    ''' Produce graphs '''
    ct = [ct_gcis,ct_repair,ct_7z,ct_gz]
    cr = [cr_gcis,cr_repair,cr_7z,cr_gz]
    dt = [dt_gcis,dt_repair,dt_7z,dt_gz]
    ctcr(exp_name,exp_size,ct,cr,output_folder_path)
    dtcr(exp_name,exp_size,dt,cr,output_folder_path)

''' argv[1] Contains the path to the results folder '''
''' argv[2] Contains the path to the results plots folder '''
if __name__ == "__main__":
    results2plot(sys.argv[1],os.sys.argv[2])
