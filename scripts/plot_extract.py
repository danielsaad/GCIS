import matplotlib.pyplot as plt
import csv
import numpy as np

def dominant_points(points):
    points = sorted(points, key=lambda x: (x[0], -x[1]))
    print(points)
    new_list = []
    for p in points:
        # print(p)
        if(new_list==[]):
            new_list.append(p)
        else:
            if p[1] < new_list[-1][1]:
                print('Including p',p)
                new_list.append(p)
    return new_list
with open('results_extract.csv') as infile:
    next(infile)
    reader = csv.reader(infile)
    for line in reader:
        repairsc_bps = []
        repairsn_bps = []
        access_gcis_us = []
        access_repairsc_us = []
        access_repairsn_us = []
        filename = line[0]
        filesize_mb = float(line[1])
        gcis_bps = [float(line[2])]
        access_gcis_us = [float(line[3])]
        for i in range(4,153,2):
            repairsc_bps.append(float(line[i]))
            access_repairsc_us.append(float(line[i+1]))

        points = zip(repairsc_bps, access_repairsc_us)
        points = dominant_points(points)
        repairsc_bps = [x for x,_ in points]
        access_repairsc_us = [x for _,x in points]

        print(access_repairsc_us[-1],flush=True)
        for i in range(152,len(line),2):
            repairsn_bps.append(float(line[i]))
            access_repairsn_us.append(float(line[i+1]))

        points = zip(repairsn_bps,access_repairsn_us)
        points = dominant_points(points)
        repairsn_bps = [x for x,_ in points]
        access_repairsn_us = [x for _,x in points]

        print(access_repairsn_us[-1],flush=True)
        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.grid(True)
        # print("bps = ",gcis_bps[i],repairsc_bps[i],repairsn_bps[i])
        # print("access = ",access_gcis_us[i],access_repairsc_us[i],access_repairsn_us[i])
        ax.set_title(filename + ' (' + "{:.2f}".format(float(filesize_mb))+'MB)')
        t_x = ax.set_xlabel("bps")
        t_x = ax.set_ylabel("Time ($\mu$s)")
        plt.semilogy(gcis_bps,access_gcis_us,label='GCIS',color='green',marker='o',markersize=3,linestyle='None')
        plt.semilogy(repairsc_bps,access_repairsc_us,label='GCC-C',color='red',marker='s',markersize=3,linestyle='--')
        plt.semilogy(repairsn_bps,access_repairsn_us,label='GCC-N',color='blue',marker='d',markersize=3,linestyle='-')
        max_bps = max([max(repairsc_bps), max(repairsn_bps)])
        print("max bps = ",max_bps) 
        upper = float(max_bps) * 1.1
        step = 0.1*float(max_bps)
        # ax.set_xticks(np.arange(0,upper,step))
        # ax.set_yticks(np.arange(1,300,50))
        ax.set_ylim(1,300)
        ax.set_xlim(0,max_bps+.1)
        # max_bps = float(max(gcis_bps[i],repair_bps[i]))
        # ax.set_xlim(0,max_bps+0.1)
        ax.legend(loc='best')
        # plt.show()
        fig.savefig(filename+'.pdf',dpi=600,format='pdf')