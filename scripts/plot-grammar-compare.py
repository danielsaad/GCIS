'''
    Compares GCIS ans Navarro's Repair using the csv obtained
    from the report data and plot the statistics
'''

import os
import sys
import csv
from matplotlib import pyplot as plt


def plot(csv_folder, plot_folder):
    gcis_files = [os.path.join(csv_folder, f) for f in os.listdir(csv_folder) if os.path.isfile(
        os.path.join(csv_folder, f)) and f.endswith('gcis-statistics.csv')]
    rpn_files = [os.path.join(csv_folder, f) for f in os.listdir(csv_folder) if os.path.isfile(
        os.path.join(csv_folder, f)) and f.endswith('rpn-statistics.csv')]

    gcis_files = sorted(gcis_files, key=str.casefold)
    rpn_files = sorted(rpn_files, key=str.casefold)

    gcis_total_rules = []
    gcis_first_level_rules = []
    repair_total_rules = []
    gcis_first_level_expansion = []
    gcis_grammar_size = []
    repair_grammar_size = []
    exp_name = []
    for g, r in zip(gcis_files, rpn_files):
        print('Gathering data from', os.path.basename(g).split('-')[0])
        with open(g, 'r') as gcis_f, open(r, 'r') as rpn_f:
            gcis_f_reader = csv.reader(gcis_f)
            rpn_f_reader = csv.DictReader(rpn_f)

            next(gcis_f_reader)
            next(gcis_f_reader)
            gcis_rules_per_level = []
            gcis_grammar_size.append(0)
            first = True
            for row in gcis_f_reader:
                if(row[0].startswith('Reduced')):
                    break
                if(first):
                    gcis_first_level_expansion.append(int(row[4]))
                    first = False
                gcis_grammar_size[-1] += int(row[4])
                gcis_rules_per_level.append(int(row[3]))
                gcis_grammar_size
            gcis_first_level_rules.append(int(gcis_rules_per_level[0]))
            gcis_total_rules.append(sum(gcis_rules_per_level))
            for row in rpn_f_reader:
                repair_total_rules.append(int(row['number_of_rules']))
                repair_grammar_size.append(repair_total_rules[-1]*2)

    for i, f in enumerate(gcis_files):
        exp_name.append(os.path.basename(f).split('-')[0])
        print(exp_name)
        print(gcis_total_rules[i], repair_total_rules[i])
        print(gcis_first_level_rules[i])
        print(gcis_first_level_expansion[i])

    plot_number_of_rules_comparison(exp_name, gcis_total_rules,
                                    repair_total_rules, plot_folder)
    plot_grammar_size_comparison(exp_name, gcis_grammar_size, gcis_total_rules,
                                 repair_grammar_size, repair_total_rules, plot_folder)
    plot_first_level_rules(exp_name, gcis_first_level_rules,
                           gcis_total_rules, plot_folder)
    plot_first_level_expansion(
        exp_name, gcis_first_level_expansion, plot_folder)


def plot_first_level_rules(exp_name, a, b, plot_folder):
    fig = plt.figure()
    ax = fig.add_subplot(111)
    gap = .3
    width = .1
    locs = [width + x*gap for x in range(len(a))]
    locs2 = [l+width for l in locs]
    ax.bar(locs, a, width=width, label='GCIS First Level', color='blue')
    ax.bar(locs2, b, width=width, label='GCIS All Levels', color='green')
    ax.set_title('Number of Rules vs Number of Rules in First Level')
    ax.set_ylabel('Number of Rules')
    ax.grid(True)

    ax.legend()
    x_labels = [x for x in exp_name]
    rects = ax.patches
    for rect, label in zip(rects, a):
        height = rect.get_height()
        ax.text(rect.get_x() + rect.get_width()/2, height, str(label), color='blue',
                rotation=90, fontweight='bold', ha='center', va='bottom', fontsize=6)

    for rect, label in zip(rects, b):
        height = rect.get_height()
        ax.text(rect.get_x() + 1.5*width, label, str(label), color='green',
                rotation=90, fontweight='bold', ha='center', va='bottom', fontsize=6)

    plt.xticks([x for x in locs], x_labels, rotation='vertical')
    # ax.legend(loc='best')

    output_path = os.path.join(
        plot_folder, 'grammar-comparison-number-of-rules-first-level-gcis.png')
    print('Generating', output_path+'\n')
    plt.savefig(output_path, format='png', bbox_inches='tight', dpi=600)
    plt.close('all')


def plot_first_level_expansion(exp_name, a, plot_folder):
    fig = plt.figure()
    ax = fig.add_subplot(111)
    gap = .2
    width = .1
    locs = [width + x*gap for x in range(len(a))]
    ax.bar(locs, a, width=width, label='GCIS First Level Expansion', color='green')
    ax.set_title('Total Length First Level Expansion (GCIS) ')
    ax.set_ylabel('Total Length')
    ax.grid(True)

    ax.legend()
    x_labels = [x for x in exp_name]
    rects = ax.patches
    for rect, label in zip(rects, a):
        height = rect.get_height()
        ax.text(rect.get_x() + rect.get_width()/2, height, str(label), color='green',
                rotation=90, fontweight='bold', ha='center', va='bottom', fontsize=6)

    plt.xticks([x for x in locs], x_labels, rotation='vertical')
    # ax.legend(loc='best')

    output_path = os.path.join(
        plot_folder, 'grammar-comparison-total-length-first-level-expansion-gcis.png')
    print('Generating', output_path+'\n')
    plt.savefig(output_path, format='png', bbox_inches='tight', dpi=600)
    plt.close('all')


def plot_number_of_rules_comparison(exp_name, a, b, plot_folder):
    fig = plt.figure()
    ax = fig.add_subplot(111)
    gap = .3
    width = .1
    locs = [width + x*gap for x in range(len(a))]
    locs2 = [x+width for x in locs]
    ax.bar(locs, a, width=width, label='GCIS', color='green')
    ax.bar(locs2, b, width=width, label='RePair', color='red')
    ax.set_title('GCIS vs RePair')
    ax.set_ylabel('Number of Rules')
    ax.grid(True)
    ax.legend()
    x_labels = [x for x in exp_name]
    rects = ax.patches
    for rect, label in zip(rects, a):
        height = rect.get_height()
        ax.text(rect.get_x() + rect.get_width()/2, height, str(label), color='green',
                rotation=90, fontweight='bold', ha='center', va='bottom', fontsize=6)

    for rect, label in zip(rects, b):
        height = rect.get_height()
        ax.text(rect.get_x() + 1.5*width, label, str(label), color='red',
                rotation=90, fontweight='bold', ha='center', va='bottom', fontsize=6)

    plt.xticks([x for x in locs], x_labels, rotation='vertical')
    # ax.legend(loc='best')

    output_path = os.path.join(
        plot_folder, 'grammar-comparison-number-of-rules.png')
    print('Generating', output_path+'\n')
    plt.savefig(output_path, format='png', bbox_inches='tight', dpi=600)
    plt.close('all')


def plot_grammar_size_comparison(exp_name, grammar_size_gcis, number_of_rules_gcis,
                                 grammar_size_repair, number_of_rules_repair, plot_folder):
    fig = plt.figure()
    ax = fig.add_subplot(111)
    gap = .3
    width = .1
    a = [x-y for (x, y) in zip(grammar_size_gcis, number_of_rules_gcis)]
    b = [x-y for (x, y) in zip(grammar_size_repair, number_of_rules_repair)]
    locs = [width + x*gap for x in range(len(a))]
    locs2 = [x+width for x in locs]
    ax.bar(locs, a, width=width, label='GCIS', color='green')
    ax.bar(locs2, b, width=width, label='RePair', color='red')
    ax.set_title('GCIS vs RePair')
    ax.set_ylabel('$\mathcal{G}-\mathcal{g}$')
    ax.grid(True)
    ax.legend()
    x_labels = [x for x in exp_name]
    rects = ax.patches
    for rect, label in zip(rects, a):
        height = rect.get_height()
        ax.text(rect.get_x() + rect.get_width()/2, height, str(label), color='green',
                rotation=90, fontweight='bold', ha='center', va='bottom', fontsize=6)

    for rect, label in zip(rects, b):
        height = rect.get_height()
        ax.text(rect.get_x() + 1.5*width, label, str(label), color='red',
                rotation=90, fontweight='bold', ha='center', va='bottom', fontsize=6)

    plt.xticks([x for x in locs], x_labels, rotation='vertical')
    # ax.legend(loc='best')

    output_path = os.path.join(
        plot_folder, 'grammar-comparison-size.png')
    print('Generating', output_path+'\n')
    plt.savefig(output_path, format='png', bbox_inches='tight', dpi=600)
    plt.close('all')


def plot_grammar_compare(csv_folder, plot_folder):
    if(not os.path.isdir(csv_folder)):
        print(csv_folder, 'is not a directory')
        exit(0)
    if(not os.path.isdir(plot_folder)):
        print('Creating plot folder', plot_folder)
        os.makedirs(plot_folder, exist_ok=True)

    plot(csv_folder, plot_folder)


if __name__ == "__main__":
    if(len(sys.argv) < 3):
        print('Usage: plot-grammar-compare <csv folder> <plot folder>')
        exit(0)
    plot_grammar_compare(sys.argv[1], sys.argv[2])
