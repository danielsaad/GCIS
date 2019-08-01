import os
import sys
import matplotlib.pyplot as plt
import re
import numpy as np
import csv
import collections


def plot_mem_peak_compress(expnames,gcisVmPeakMB,repairVmPeakMB,output_folder_path):
    fig = plt.figure()
    ax = fig.add_subplot(111)
    gap =  .3
    width = .1
    locs = [width + x*gap for x in range(len(gcisVmPeakMB))]
    locs2 = [x+width for x in locs]
    ax.bar(locs, gcisVmPeakMB, width=width, label='GCIS',color='green')
    ax.bar(locs2,repairVmPeakMB,width=width,label='RePair',color='red')
    max_y = max([max(gcisVmPeakMB),max(repairVmPeakMB)])
    ax.set_yticks(np.arange(0,max_y),max_y/10.0)
    ax.set_title('Memory Peak for Compression')
    ax.set_ylabel("Memory Peak (MB)")
    ax.grid(True)
    gcisLabels = ['GCIS' for x in gcisVmPeakMB]
    repairLabels = ['RePair' for x in repairVmPeakMB]
    plt.xticks([x+width/2.0 for x in locs], expnames,rotation='vertical')
    ax.legend(loc='best')
    output_path = os.path.join(output_folder_path,'memory-peak-compression.pdf')
    plt.savefig(output_path,format='pdf',bbox_inches='tight',dpi=600)

    print('GCIS COMPRESS PEAK',gcisVmPeakMB)
    print('REPAIR COMPRESS PEAK',repairVmPeakMB)


def plot_mem_peak_decompress(expnames,gcisVmPeakMB,repairVmPeakMB,output_folder_path):
    fig = plt.figure()
    ax = fig.add_subplot(111)
    gap =  .3
    width = .1
    locs = [width + x*gap for x in range(len(gcisVmPeakMB))]
    locs2 = [x+width for x in locs]
    ax.bar(locs, gcisVmPeakMB, width=width, label='GCIS',color='green')
    ax.bar(locs2,repairVmPeakMB,width=width,label='RePair',color='red')
    max_y = max([max(gcisVmPeakMB),max(repairVmPeakMB)])
    ax.set_yticks(np.arange(0,max_y),max_y/10.0)
    ax.set_title('Memory Peak for Decompression')
    ax.set_ylabel("Memory Peak (MB)")
    ax.grid(True)
    gcisLabels = ['GCIS' for x in gcisVmPeakMB]
    repairLabels = ['RePair' for x in repairVmPeakMB]
    plt.xticks([x+width/2.0 for x in locs], expnames,rotation='vertical')
    ax.legend(loc='best')
    output_path = os.path.join(output_folder_path,'memory-peak-decompression.pdf')
    plt.savefig(output_path,format='pdf',bbox_inches='tight',dpi=600)

    print('GCIS DECOMPRESS PEAK',gcisVmPeakMB)
    print('REPAIR DECOMPRESS PEAK',repairVmPeakMB)



''' Plot the peak memory for GCIS and RePair Regarding all experiments '''
def memory2plot(input_folder_path,output_folder_path):
    os.makedirs(output_folder_path,exist_ok=True)
    
    files = os.listdir(input_folder_path)
#    print(files)
    regexcompress = re.compile('.*-compress.csv$')
    regexdecompress = re.compile('decompress.csv$')
    filescompress = [f for f in files if re.search(regexcompress,f)]
    filesdecompress = [f for f in files if re.search(regexdecompress,f)]

    filescompress = sorted(filescompress,key=lambda s: s.lower())
    filesdecompress = sorted(filesdecompress,key=lambda s: s.lower())

    vmPeakMB = []
    for f in filescompress:
        with open(os.path.join(input_folder_path,f),'r') as csvfile:
            reader = csv.reader(csvfile,delimiter=';')
            vmPeakMB.append(round(float(collections.deque(reader, 1)[0][2])/1000000,2))

    gcisVmPeakMB = [vpeak for i,vpeak in enumerate(vmPeakMB) if i%2==0]
    repairVmPeakMB = [vpeak for i,vpeak in enumerate(vmPeakMB) if i%2!=0]
    expnames = list(set([f.split('-')[0] for f in filescompress]))
    expnames = sorted(expnames,key=lambda s: s.lower())


    plot_mem_peak_compress(expnames,gcisVmPeakMB,repairVmPeakMB,output_folder_path)

    vmPeakMB = []
    for f in filesdecompress:
        with open(os.path.join(input_folder_path,f),'r') as csvfile:
            reader = csv.reader(csvfile,delimiter=';')
            vmPeakMB.append(round(float(collections.deque(reader, 1)[0][2])/1000000,2))

    gcisVmPeakMB = [vpeak for i,vpeak in enumerate(vmPeakMB) if i%2==0]
    repairVmPeakMB = [vpeak for i,vpeak in enumerate(vmPeakMB) if i%2!=0]

    plot_mem_peak_decompress(expnames,gcisVmPeakMB,repairVmPeakMB,output_folder_path)



''' argv[1] = csv folder path containing the memory usage of the experiments '''
''' argv[2] = folder path containing the memory usage plots '''
if __name__ == '__main__':
    memory2plot(sys.argv[1],sys.argv[2])
