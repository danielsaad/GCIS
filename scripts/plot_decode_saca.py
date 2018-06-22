from matplotlib import pyplot as plt
import numpy
import os
import subprocess
import sys
import csv

class experiment_data:
    def __init__(self):
        self.experiment_name = []
        self.experiment_value = []
        self.input_name = []

def plot(e,output_folder):
    gap = .2
    width = .1
    colors = ['green','red','blue','orange','yellow']
    number_of_bars = len(e.experiment_name)
    labels = e.experiment_name
    print(labels)
    locs = [x for x in [[gap*loc] for loc in range(number_of_bars)]]

    for i in range(len(e.experiment_value)):
        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.set_title('Suffix Array Construction'+'\n' + e.input_name[i])
        ax.set_ylabel("Time(s)")
        ax.set_xticks([])
        ax.grid()
        values = [x for x in [[exp] for exp in e.experiment_value[i]]]
        print(values)
        for j in range(len(locs)):
            ax.bar(locs[j],values[j], width=width, label=labels[j],color=colors[j])
        ax.legend(bbox_to_anchor=(1.04,1), loc="upper left")
        plt.savefig(os.path.join(output_folder,e.input_name[i]+'.png'),format='png',bbox_inches='tight',dpi=600)


def setup(results_folder):
    
    decode_saca_file = os.path.join(results_folder,'results-saca.csv')
    
    if(not os.path.isfile(decode_saca_file)):
        print('Decode SACA file not available')
        exit(1)

    e = experiment_data()
    with open(decode_saca_file,'r') as csv_file:
        reader = csv.reader(csv_file)
        l = next(reader)
        e.experiment_name = l[1:] 
        # print(e.experiment_name)
        for l in reader:
            e.input_name.append(l[0])
            # print(e.input_name)
            e.experiment_value.append([float(x) for x in l[1:]])
            # print(e.experiment_value)
    return e
        
def plot_decode_saca(results_folder):
    experiment_data = setup(results_folder)
    plot(experiment_data,results_folder)

if __name__ == '__main__':
    plot_decode_saca(sys.argv[1])
