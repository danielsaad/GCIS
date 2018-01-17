import os
import sys
from results2plot import results2plot
from report2csv import report2csv
from statistics2plot import statistics2plot
from memory2plot import memory2plot


def generate_graphs(output_folder_path):
    codec_plot_folder_path= os.path.join(output_folder_path,'Codec-Plots')
    report_folder_path = os.path.join(output_folder_path,'GCIS-Reports')
    statistics_csv_folder_path = os.path.join(output_folder_path,*['Statistics','csv'])
    statistics_plot_folder_path = os.path.join(output_folder_path,*['Statistics','plots'])
    memory_csv_folder_path = os.path.join(output_folder_path,*['Memory-Monitor','csv'])
    memory_plot_folder_path = os.path.join(output_folder_path,*['Memory-Monitor','plots'])
    results2plot(output_folder_path,codec_plot_folder_path)
    memory2plot(memory_csv_folder_path,memory_plot_folder_path)
    report2csv(report_folder_path,statistics_csv_folder_path)
    statistics2plot(statistics_csv_folder_path,statistics_plot_folder_path)

''' argv[1] = results folder path '''
if __name__ == "__main__":
    generate_graphs(sys.argv[1])
