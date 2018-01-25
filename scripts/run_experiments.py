import sys
import os
import time
import filecmp
import csv
import gz
import sevenz
import gcis
import repair



enable_benchmark = True

enable_gz = True
enable_7z = True
enable_gc_is = True
enable_repair = True
enable_testing = True

enable_statistics = True


def calc_compression_ratio(path_f1,path_f2):
    st1 = os.stat(path_f1)
    st2 = os.stat(path_f2)
    compression_ratio = (st2.st_size/st1.st_size)
    return compression_ratio

def calc_file_size(path):
    st = os.stat(path)
    return st.st_size


def get_codec(name):
    if(name=="gc-is"):
        return(gcis.compress_gc_is,gcis.decompress_gc_is,".gcis")
    if(name=='gc-is-statistics'):
        return(gcis.compress_gc_is_statistics,gcis.decompress_gc_is_statistics,".gcis")
    if(name=="7z"):
        return(sevenz.compress_7z,sevenz.decompress_7z,".7z")
    if(name=='repair'):
        return(repair.compress_repair,repair.decompress_repair,'')
    else:
        return(gz.compress_gz,gz.decompress_gz,".gz")


class Experiment_data:

    """ Class that defines the experiment data """

    def __init__(self,experiment_name,input_size,dictionary_size,compression_ratio,compression_time,decompression_time,
                 codec_working):
        self.experiment_name = experiment_name
        self.compression_ratio = compression_ratio
        self.dictionary_size = dictionary_size
        self.input_size = input_size
        self.compression_time = compression_time
        self.decompression_time = decompression_time
        self.codec_working = codec_working



def run_experiment(file_path,exp_name):
    ''' Repair generates a different output (2 files instead of one)
        Also the compression ratio is calculated by analysing stdout, and not the compressed
        file'''

    if(exp_name =='repair'):
        return run_repair(file_path)
    # GC-IS, 7Z and GZ have the same experimental framework
    base,filename = os.path.split(file_path)

    #Create new directory to store temporary files
    compressed_files_dir = os.path.join(base,"compressedFiles")
    os.makedirs(compressed_files_dir,exist_ok=True)
    compress,decompress,extension = get_codec(exp_name)

    input_size = calc_file_size(file_path)

    compressed_file_path = os.path.join(compressed_files_dir,filename+extension)
    start = time.time()
    print("Compressing ("+exp_name+')',file_path,'->',compressed_file_path)
    compress(file_path,os.path.join(compressed_files_dir,compressed_file_path))
    end = time.time()
    compression_time = end-start
    dictionary_size = calc_file_size(compressed_file_path)

    decompressed_file_path = os.path.join(compressed_files_dir,os.path.splitext(compressed_file_path)[0])
    print("Decompressing ("+exp_name+')',compressed_file_path,'->',decompressed_file_path)
    start = time.time()
    decompress(compressed_file_path,decompressed_file_path)
    end = time.time()
    decompression_time = end-start

    compression_ratio = calc_compression_ratio(file_path,compressed_file_path)

    codec_working = True
    if(enable_testing):
        if(not (filecmp.cmp(file_path,decompressed_file_path,shallow=False))):
            print("Codec Failed: decompression of ",compressed_file_path,'differs from',file_path)
            codec_working = False
        else:
            print("Codec Succeeded: decompression of",compressed_file_path,'matches',file_path)


    os.remove(compressed_file_path)
    os.remove(decompressed_file_path)

    return Experiment_data(exp_name,input_size,dictionary_size,compression_ratio,compression_time,
                           decompression_time,codec_working)

def run_repair(file_path):
    base, filename = os.path.split(file_path)
    compress,decompress,extension = get_codec('repair')

    input_size = calc_file_size(file_path)

    print("Compressing (repair)",file_path,'->',file_path+'.[prel|seq]')
    start = time.time()
    compress(file_path)
    end = time.time()
    compression_time = end-start

    print("Decompressing (repair)",file_path+'.[prel|seq]','->',file_path)
    start = time.time()
    decompress(file_path)
    end = time.time()
    decompression_time = end-start

    dictionary_size = calc_file_size(file_path+'.prel') + calc_file_size(file_path+'.seq')
    compression_ratio = dictionary_size/calc_file_size(file_path)
    codec_working = True;
    os.remove(file_path+'.prel')
    os.remove(file_path+'.seq')
    os.remove(file_path+'.u')
    return Experiment_data('repair',input_size,dictionary_size,
        compression_ratio,compression_time,decompression_time,codec_working)

def setup_csv(csv_writer):
    header = ["Experiment"]
    header.append("Input Size (MB)")
    header_template = ["Dictionary Size (MB)","Compression Ratio","Compress Time","Decompress Time","Codec Ok"]
    if(enable_gc_is):
        for h in header_template:
            header.append(h+" (GC-IS)")
    if(enable_repair):
        for h in header_template:
            header.append(h+' (RePair)')
    if(enable_7z):
        for h in header_template:
            header.append(h+" (7z)")
    if(enable_gz):
        for h in header_template:
            header.append(h+" (gz)")

    csv_writer.writerow(header)

def get_row(ed):
    return [str(ed.dictionary_size/1000000),str(100 * ed.compression_ratio),
            str(ed.compression_time),str(ed.decompression_time),
            str(ed.codec_working)]

def codec(experiment_path,output_path):

    os.makedirs(output_path,exist_ok=True)

    csvfilename = os.path.join(output_path,'results.csv');
    csvfile = open(csvfilename,'w')
    writer = csv.writer(csvfile, delimiter=',')

    # get all regular files in experiments folder.
    files = [os.path.join(experiment_path,f) for f in os.listdir(experiment_path) if os.path.isfile(os.path.join(experiment_path,f))]
    files.sort()
    print(files)

    setup_csv(writer)

    # run experiments for each file
    for f in files:
        filename = os.path.split(f)[1]
        metadata_row = [filename,str(calc_file_size(f)/1000000)]
        if(enable_gc_is):
            experiment_data_gc_is = run_experiment(f,"gc-is")
            gc_is_row = get_row(experiment_data_gc_is)
        if(enable_7z):
            experiment_data_7z = run_experiment(f, "7z")
            sevenz_row = get_row(experiment_data_7z)
        if(enable_gz):
            experiment_data_gz = run_experiment(f, "gz")
            gz_row = get_row(experiment_data_gz)
        if(enable_repair):
            experiment_data_repair = run_experiment(f,'repair')
            repair_row = get_row(experiment_data_repair)

        line = metadata_row
        line += gc_is_row if enable_gc_is  else  []
        line += repair_row if enable_repair else []
        line += sevenz_row if enable_7z else  []
        line += gz_row if enable_gz else  []
        writer.writerow(line)
        csvfile.flush()

    csvfile.close()


''' Move .log files to the correct destination '''
def move_report_files(filename,report_folder_path):
    # Generate dictionary statistics from rerport 2 csv file
    # Check report.log exists
    report_file_path = 'report.log'
    if(os.path.isfile(report_file_path)):
        # rename the file to /output/path/report_file
        report_filename = filename+'.log'
        report_file_path = os.path.join(report_folder_path,report_filename)
        os.rename('report.log',report_file_path)
        print(report_filename, "generated")
    else:
        print("No report files were found")

''' Move memory_monitor files to the correct destination '''
def move_memory_monitor_files(filename,memory_folder_path,expname):
    # Generate Memory Footprint Graph
    mem_ftp_path = 'mem-mon-out.csv'
    # check if memory footprint csv exists
    if(os.path.isfile(mem_ftp_path)):
        print(mem_ftp_path,'generated')
        # Move the csv file and the plot to the ouput folder
        os.rename('mem-mon-out.csv',os.path.join(
            memory_folder_path,filename + '-' + expname + '.csv'))
    else:
        print("Memory footprint file not found")

''' Runs GC-IS and repair to collect GC-IS
    encode statistics and GC-IS/RePair Memory Consumption
    over time.
    experiment_path = path to the folder containing the texts
    output_path = path to the root of the reusults folder '''
def codec_statistics(experiment_path,output_path):
    os.makedirs(output_path,exist_ok=True)
    files = [os.path.join(experiment_path,f) for f in os.listdir(experiment_path) if os.path.isfile(os.path.join(experiment_path,f))]
    files.sort()
    print(files)

    ''' Create the GCIS-Reports Folder, where the enconding reports of
        GC-IS will be stored '''
    report_folder_path = os.path.join(output_path,'GCIS-Reports')
    memory_folder_path = os.path.join(output_path,*['Memory-Monitor','csv'])

    os.makedirs(report_folder_path,exist_ok=True)
    os.makedirs(memory_folder_path,exist_ok=True)

    for f in files:
        # GC-IS compressed file path
        filename = os.path.basename(f)
        gcis_output_path = os.path.join(output_path,'output.gc-is')
        # GC-IS decoded file path
        gcis_decode_path = os.path.join(output_path, os.path.split(f)[1])
        # Run GC-IS
        gcis.compress_gc_is_statistics(f,gcis_output_path)
        move_report_files(filename,report_folder_path)
        move_memory_monitor_files(filename,memory_folder_path,'gcis-compress')

        gcis.decompress_gc_is_statistics(gcis_output_path,gcis_decode_path)
        move_memory_monitor_files(filename,memory_folder_path,'gcis-decompress')

        # Remove both compressed and decompressed temporary files from GCIS
        os.remove(gcis_output_path)
        os.remove(gcis_decode_path)


        ''' Run Repair '''
        repair.compress_repair_statistics(f,output_path)
        move_memory_monitor_files(filename,memory_folder_path,'repair-compress')

        repair.decompress_repair_statistics(f,output_path)
        move_memory_monitor_files(filename,memory_folder_path,'repair-decompress')

        # Remove compressed temporary RePair files
        os.remove(f+'.prel')
        os.remove(f+'.seq')
        os.remove(f+'.u')



def run_experiments(input_folder_path,output_folder_path):
    if(enable_benchmark):
        codec(input_folder_path,output_folder_path)
    if(enable_statistics):
        codec_statistics(input_folder_path,output_folder_path)

#argv[1] contains the path to the experiments folder
#argv[2] contains the path to the output folder

if __name__ == "__main__":
    run_experiments(sys.argv[1],sys.argv[2])
