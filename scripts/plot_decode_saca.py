from matplotlib import pyplot as plt
import matplotlib.patches as mpatches

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


def plot(e, output_folder):
    gap = .2
    width = .08
    colors = ['black', 'green', 'red', 'blue',
              'orange', 'yellow', 'blue', 'blue', 'blue']
    labels = e.experiment_name
    for i,l in enumerate(labels):
        if(l.startswith('SAIS-YUTA')):
            labels[i] = 'SAIS'
        if(l.startswith('SAIS-DIVSUFSORT')):
            labels[i] = 'divsufsort'
    print(labels)
    locs = [x for x in [[gap*loc] for loc in range(len(colors))]]
    print('Locs = ', locs)
    for i in range(len(e.experiment_value)):
        print(i)
        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.set_title('Suffix Array Construction'+'\n' + e.input_name[i])
        ax.set_ylabel("Time(s)")
        ax.set_xticks([])
        ax.grid()
        values = [x for x in [[exp] for exp in e.experiment_value[i]]]

        # GCIS Decode
        ax.bar(locs[0], values[0], width=width, label=labels[0], color='black')

        # # GCIS SAIS
        # ax.bar(locs[1],values[1],width=width,label=labels[1],color='green')

        # # SAIS NONG
        # ax.bar(locs[2],values[0],width=width,color='black')
        # ax.bar(locs[2],values[3],width=width,label=labels[3],bottom=values[0],color='red')

        # # SAIS YUTA
        # ax.bar(locs[3],values[0],width=width,color='black')
        # ax.bar(locs[3],values[4],width=width,label=labels[4],bottom=values[0],color='yellow')

        # # SAIS DIVSUFSORT
        # ax.bar(locs[4],values[0],width=width,color='black')
        # ax.bar(locs[4],values[5],width=width,label=labels[5],bottom=values[0],color='orange')

        # GCIS SAIS + LCP
        ax.bar(locs[1], values[1], width=width, label=labels[1], color='green')
        value = [values[2][0] - values[1][0]]
        ax.bar(locs[1], value, width=width,
               bottom=values[1], color='green', hatch='////',edgecolor='black',linewidth=0)

        # SAIS Yuta + LCP
        ax.bar(locs[2], values[0], width=width, color='black')
        ax.bar(locs[2], [values[13][0]-values[0][0]], width=width,
               bottom=values[0], label=labels[4], color='yellow')
        value = [values[19][0] - values[13][0]]
        ax.bar(locs[2], value, width=width, bottom=[
               values[13][0]], color='yellow', hatch='////',edgecolor='black',linewidth=0)

        # SAIS DIVSUFSORT + LCP

        ax.bar(locs[3], values[0], width=width, color='black')
        ax.bar(locs[3], [values[16][0]-values[0][0]], width=width,
               bottom=values[0], label=labels[5], color='orange')
        if(values[22][0] > 0):
            value = [values[22][0] - values[16][0]]
            ax.bar(locs[3], value, width=width, bottom=[
                   values[16][0]], color='orange', hatch='////',edgecolor='black',linewidth=0)

        h, l = ax.get_legend_handles_labels()
        white_hatch = mpatches.Patch(
            facecolor='white', alpha=0.6, hatch=r'////', label='LCP-CALC')

        h += [white_hatch]
        ax.legend(frameon=False, handles=h, loc="upper left")
        plt.savefig(os.path.join(
            output_folder, e.input_name[i]+'.pdf'), format='pdf', bbox_inches='tight', dpi=600)


def setup(results_folder):

    decode_saca_file = os.path.join(results_folder, 'results-saca.csv')

    if(not os.path.isfile(decode_saca_file)):
        print('Decode SACA file not available')
        exit(1)

    e = experiment_data()
    with open(decode_saca_file, 'r') as csv_file:
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
    plot(experiment_data, results_folder)


'''
    argv[1]: results folder
'''
if __name__ == '__main__':
    plot_decode_saca(sys.argv[1])
