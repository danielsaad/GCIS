import sys
import os
import matplotlib.pyplot as plt
import csv


def barplot(title,y_title,file_suffix,level,y,output_folder_path):
    fig = plt.figure()
    ax = fig.add_subplot(111)
    gap =  .3
    width = .1
    locs = [width + x*gap for x in range(len(y))]
    ax.bar(locs, y, width=width)
    ax.set_title(title)
    ax.set_ylabel(y_title)
    ax.grid(True)
    x_labels = ['Level '+ str(l) for l in level]

    rects = ax.patches
    for rect, label in zip(rects, y):
        height = rect.get_height()
        ax.text(rect.get_x() + rect.get_width()/2, height , str(label),color='blue', fontweight='bold', ha='center', va='bottom',fontsize=6)

    plt.xticks([x for x in locs], x_labels,rotation='vertical')
    #ax.legend(loc='best')

    output_path = os.path.join(output_folder_path,title+'-'+file_suffix+'.png')
    print('Generating',output_path+'\n')
    plt.savefig(output_path,format='png',bbox_inches='tight',dpi=300)
    plt.close('all')

def statistics2plot_helper(statistics_file_path,output_folder_path):
    os.makedirs(output_folder_path,exist_ok=True)
    expname = os.path.split(statistics_file_path)[1].split('-')[0]
    level = []
    string_size = []
    alphabet_size = []
    number_of_rules = []
    avg_rule_len = []
    avg_lcp_len = []
    avg_rule_suffix_len = []
    rle_potential = []
    avg_rle  = []
    dict_size = []

    with open(statistics_file_path,'r') as csvfile:
        csvreader = csv.reader(csvfile)
        # Skip the first two lines (header)
        next(csvreader,None)
        next(csvreader,None)
        for row in csvreader:
            if(row[0] == 'Reduced String Length'):
                break
            level.append(int(row[0]))
            string_size.append(int(row[1]))
            alphabet_size.append(int(row[2]))
            number_of_rules.append(int(row[3]))
            avg_rule_len.append(float(row[4]))
            avg_lcp_len.append(float(row[5]))
            avg_rule_suffix_len.append(float(row[7]))
            rle_potential.append(int(row[10]))
            avg_rle.append(float(row[11]))
            dict_size.append(round(float(row[14])/1000000.0,2))


    barplot(expname,'String Size','0-strsz',level,string_size,output_folder_path)
    barplot(expname,'Alphabet Size','1-alphsize',level,alphabet_size,output_folder_path)
    barplot(expname,'Number of Rules','2-number-of-rules',level,number_of_rules,output_folder_path)
    barplot(expname,'Average Rule Length','3-avg-rule-len',level,avg_rule_len,output_folder_path)
    barplot(expname,'Average LCP Length','4-avg-lcp-len',level,avg_lcp_len,output_folder_path)
    barplot(expname,'Average Rule Suffix Length','5-avg-rule-suffix-len',level,avg_rule_suffix_len,output_folder_path)
    barplot(expname,'Total RLE Length','6-rle-potential',level,rle_potential,output_folder_path)
    barplot(expname,'Average RLE Length per Rule Suffix','7-avg-rle',level,avg_rle,output_folder_path)
    barplot(expname,'Dictionary Size (MB)','8-dict-size',level,dict_size,output_folder_path)

def statistics2plot(statistics_folder_path,output_folder_path):
    files = [os.path.join(statistics_folder_path,f) for f in os.listdir(statistics_folder_path)]
    files.sort()
    for f in files:
        print(f)
        print('Processing',os.path.basename(f))
        statistics2plot_helper(f,output_folder_path)

''' argv[1] folder path to the csv files containing the GC-IS compression statistics
    argv[2] folder path to the csv files containing the GC-IS compression statistics graphs '''
if __name__ == '__main__':
    scan(sys.argv[1],sys.argv[2])
