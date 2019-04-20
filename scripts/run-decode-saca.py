import os
import sys
import time
import csv
import gcis
import sais
import statistics

number_of_samples = 1


class experiment_data:
    def __init__(self):
        self.experiment_name = ""
        self.experiment_time = []
        self.experiment_input = []


def run_decode_nong(input_folder_path, output_folder_path):
    e = experiment_data()
    e.experiment_name = "DECODE-SAIS-NONG"
    input_filename = os.listdir(input_folder_path)
    for f in input_filename:
        input_file = os.path.join(input_folder_path, f)
        output_file = os.path.join(output_folder_path, f+'.decode.nong.sa')
        compressed_file = os.path.join(output_folder_path, f+'.gcis')

        print("GCIS Compressing: ", f)
        if(os.path.isfile(compressed_file)):
            print(compressed_file, 'already exists, skipping')
        else:
            gcis.compress_gc_is(input_file, compressed_file)

        print("GCIS Decoding and SAIS Nong")
        dt = st = tt = []
        for i in range(number_of_samples):
            (decode_time, saca_time) = sais.decode_sais_nong(
                compressed_file, output_file)
            dt.append(float(decode_time))
            st.append(float(saca_time))
            tt.append(float(decode_time)+float(saca_time))
            print("Decode Time =", decode_time, 'seconds')
            print("Saca Nong Time =", saca_time, 'seconds')

        e.experiment_input.append(f)
        e.experiment_time.append(('{:.2f}'.format(statistics.mean(dt)), '{:.2f}'.format(
            statistics.mean(st)), '{:.2f}'.format(statistics.mean(tt))))

    return e


def run_decode_yuta(input_folder_path, output_folder_path):
    e = experiment_data()
    e.experiment_name = "DECODE-SAIS-YUTA"
    input_filename = os.listdir(input_folder_path)
    for f in input_filename:
        input_file = os.path.join(input_folder_path, f)
        output_file = os.path.join(output_folder_path, f+'.decode.yuta.sa')
        compressed_file = os.path.join(output_folder_path, f+'.gcis')

        print("GCIS Compressing: ", f)
        if(os.path.isfile(compressed_file)):
            print(compressed_file, 'already exists, skipping')
        else:
            gcis.compress_gc_is(input_file, compressed_file)

        print("GCIS Decoding and SAIS Yuta")
        dt = st = tt = []
        for i in range(number_of_samples):
            (decode_time, saca_time) = sais.decode_sais_yuta(
                compressed_file, output_file)
            dt.append(float(decode_time))
            st.append(float(saca_time))
            tt.append(float(decode_time)+float(saca_time))
            print("Decode Time =", decode_time, 'seconds')
            print("Saca Yuta Time =", saca_time, 'seconds')

        e.experiment_input.append(f)
        e.experiment_time.append(('{:.2f}'.format(statistics.mean(dt)), '{:.2f}'.format(
            statistics.mean(st)), '{:.2f}'.format(statistics.mean(tt))))

    return e


def run_decode_lcp_yuta(input_folder_path, output_folder_path):
    e = experiment_data()
    e.experiment_name = "DECODE-SAIS-LCP-YUTA"
    input_filename = os.listdir(input_folder_path)
    for f in input_filename:
        input_file = os.path.join(input_folder_path, f)
        output_file = os.path.join(output_folder_path, f+'.decode.yuta.salcp')
        compressed_file = os.path.join(output_folder_path, f+'.gcis')

        print("GCIS Compressing: ", f)
        if(os.path.isfile(compressed_file)):
            print(compressed_file, 'already exists, skipping')
        else:
            gcis.compress_gc_is(input_file, compressed_file)

        print("GCIS Decoding and SAIS Yuta + LCP")
        dt = st = tt = []
        for i in range(number_of_samples):
            (decode_time, saca_time) = sais.decode_sais_lcp_yuta(
                compressed_file, output_file)
            dt.append(float(decode_time))
            st.append(float(saca_time))
            tt.append(float(decode_time)+float(saca_time))
            print("Decode Time =", decode_time, 'seconds')
            print("Saca LCP Yuta Time =", saca_time, 'seconds')

        e.experiment_input.append(f)
        e.experiment_time.append(('{:.2f}'.format(statistics.mean(dt)), '{:.2f}'.format(
            statistics.mean(st)), '{:.2f}'.format(statistics.mean(tt))))

    return e


def run_decode_divsufsort(input_folder_path, output_folder_path):
    e = experiment_data()
    e.experiment_name = "DECODE-SAIS-DIVSUFSORT"
    input_filename = os.listdir(input_folder_path)
    for f in input_filename:
        input_file = os.path.join(input_folder_path, f)
        output_file = os.path.join(
            output_folder_path, f+'.decode.divsufsort.sa')
        compressed_file = os.path.join(output_folder_path, f+'.gcis')

        print("GCIS Compressing: ", f)
        if(os.path.isfile(compressed_file)):
            print(compressed_file, 'already exists, skipping')
        else:
            gcis.compress_gc_is(input_file, compressed_file)

        print("GCIS Decoding and Divsufsort")
        dt = st = tt = []
        for i in range(number_of_samples):
            (decode_time, saca_time) = sais.decode_sais_divsufsort(
                compressed_file, output_file)
            dt.append(float(decode_time))
            st.append(float(saca_time))
            tt.append(float(decode_time)+float(saca_time))
            print("Decode Time =", decode_time, 'seconds')
            print("Saca Divsufsort Time =", saca_time, 'seconds')

        e.experiment_input.append(f)
        e.experiment_time.append(('{:.2f}'.format(statistics.mean(dt)), '{:.2f}'.format(
            statistics.mean(st)), '{:.2f}'.format(statistics.mean(tt))))

    return e


def run_decode_divsufsort_lcp(input_folder_path, output_folder_path):
    e = experiment_data()
    e.experiment_name = "DECODE-SAIS-DIVSUFSORT-LCP"
    input_filename = os.listdir(input_folder_path)
    for f in input_filename:
        input_file = os.path.join(input_folder_path, f)
        output_file = os.path.join(
            output_folder_path, f+'.decode.divsufsort.salcp')
        compressed_file = os.path.join(output_folder_path, f+'.gcis')

        print("GCIS Compressing: ", f)
        if(os.path.isfile(compressed_file)):
            print(compressed_file, 'already exists, skipping')
        else:
            gcis.compress_gc_is(input_file, compressed_file)

        print("GCIS Decoding and Divsufsort + LCP")
        dt = st = tt = []
        for i in range(number_of_samples):
            (decode_time, saca_time) = sais.decode_sais_divsufsort_lcp(
                compressed_file, output_file)
            dt.append(float(decode_time))
            st.append(float(saca_time))
            tt.append(float(decode_time)+float(saca_time))
            print("Decode Time =", decode_time, 'seconds')
            print("Saca Divsufsort Time =", saca_time, 'seconds')

        e.experiment_input.append(f)
        e.experiment_time.append(('{:.2f}'.format(statistics.mean(dt)), '{:.2f}'.format(
            statistics.mean(st)), '{:.2f}'.format(statistics.mean(tt))))

    return e


def run_sais_nong(input_folder_path, output_folder_path):
    e = experiment_data()
    e.experiment_name = "SAIS-NONG"
    input_filename = os.listdir(input_folder_path)
    for f in input_filename:
        input_file = os.path.join(input_folder_path, f)
        output_file = os.path.join(output_folder_path, f+'.nong.sa')

        print("SACA Nong: ", f)
        tt = []
        for i in range(number_of_samples):
            start_time = time.perf_counter()  # Start time
            sais.sais_nong(input_file, output_file)
            end_time = time.perf_counter()  # End time
            tt.append(end_time - start_time)
            print("Time Nong =", "{:.2f}".format(end_time-start_time))
        e.experiment_input.append(f)
        e.experiment_time.append("{:.2f}".format(statistics.mean(tt)))

    return e


def run_sais_yuta(input_folder_path, output_folder_path):
    e = experiment_data()
    e.experiment_name = "SAIS-YUTA"
    input_filename = os.listdir(input_folder_path)
    for f in input_filename:
        input_file = os.path.join(input_folder_path, f)
        output_file = os.path.join(output_folder_path, f+'.yuta.sa')

        print("SACA Yuta: ", f)
        tt = []
        for i in range(number_of_samples):
            start_time = time.perf_counter()  # Start time
            sais.sais_yuta(input_file, output_file)
            end_time = time.perf_counter()  # End time
            print("Time Yuta =", "{:.2f}".format(end_time-start_time))
            tt.append(end_time-start_time)

        e.experiment_input.append(f)
        e.experiment_time.append("{:.2f}".format(statistics.mean(tt)))

    return e


def run_sais_lcp_yuta(input_folder_path, output_folder_path):
    e = experiment_data()
    e.experiment_name = "SAIS-LCP-YUTA"
    input_filename = os.listdir(input_folder_path)
    for f in input_filename:
        input_file = os.path.join(input_folder_path, f)
        output_file = os.path.join(output_folder_path, f+'.yuta.salcp')

        print("SACA Yuta + LCP: ", f)
        tt = []
        for i in range(number_of_samples):
            start_time = time.perf_counter()  # Start time
            sais.sais_lcp_yuta(input_file, output_file)
            end_time = time.perf_counter()  # End time
            print("Time Yuta =", "{:.2f}".format(end_time-start_time))
            tt.append(end_time-start_time)

        e.experiment_input.append(f)
        e.experiment_time.append("{:.2f}".format(statistics.mean(tt)))

    return e


def run_sais_divsufsort(input_folder_path, output_folder_path):
    e = experiment_data()
    e.experiment_name = "SAIS-DIVSUFSORT"
    input_filename = os.listdir(input_folder_path)
    for f in input_filename:
        input_file = os.path.join(input_folder_path, f)
        output_file = os.path.join(output_folder_path, f+'.divsufsort.sa')

        print("SACA Divsufsort: ", f)
        tt = []
        for i in range(number_of_samples):
            start_time = time.perf_counter()  # Start time
            sais.sais_divsufsort(input_file, output_file)
            end_time = time.perf_counter()  # End time
            print("Time Divsufsort =", "{:.2f}".format(end_time-start_time))
            tt.append(end_time-start_time)

        e.experiment_input.append(f)
        e.experiment_time.append("{:.2f}".format(statistics.mean(tt)))

    return e


def run_sais_divsufsort_lcp(input_folder_path, output_folder_path):
    e = experiment_data()
    e.experiment_name = "SAIS-DIVSUFSORT-LCP"
    input_filename = os.listdir(input_folder_path)
    for f in input_filename:
        input_file = os.path.join(input_folder_path, f)
        output_file = os.path.join(output_folder_path, f+'.divsufsort.salcp')

        print("SACA Divsufsort: ", f)
        tt = []
        for i in range(number_of_samples):
            start_time = time.perf_counter()  # Start time
            sais.sais_divsufsort_lcp(input_file, output_file)
            end_time = time.perf_counter()  # End time
            print("Time Divsufsort + LCP =",
                  "{:.2f}".format(end_time-start_time))
            tt.append(end_time-start_time)

        e.experiment_input.append(f)
        e.experiment_time.append("{:.2f}".format(statistics.mean(tt)))

    return e


def run_decode(input_folder_path, output_folder_path):
    e = experiment_data()
    e.experiment_name = "GCIS-DECODE"
    input_filename = os.listdir(input_folder_path)
    for f in input_filename:
        input_file = os.path.join(input_folder_path, f)
        compressed_file = os.path.join(output_folder_path, f+'.gcis')
        decompressed_file = os.path.join(output_folder_path, f+'.tmp.txt')

        print("GCIS Compressing: ", f)
        if(os.path.isfile(compressed_file)):
            print(compressed_file, 'already exists, skipping')
        else:
            gcis.compress_gc_is(input_file, compressed_file)

        print("GCIS Decoding")
        tt = []
        for i in range(number_of_samples):
            total_time = gcis.decompress_gc_is(
                compressed_file, decompressed_file)
            tt.append(total_time)
            print("Time Decode =", "{:.2f}".format(total_time))
            print('Removing tmp file', decompressed_file)
            os.remove(decompressed_file)
        e.experiment_input.append(f)
        e.experiment_time.append("{:.2f}".format(statistics.mean(tt)))

    return e


def run_decode_saca(input_folder_path, output_folder_path):

    e = experiment_data()
    e.experiment_name = "GCIS-SACA"
    input_filename = os.listdir(input_folder_path)
    for f in input_filename:
        input_file = os.path.join(input_folder_path, f)
        compressed_file = os.path.join(output_folder_path, f+'.gcis')
        output_file = os.path.join(output_folder_path, f+'.gics.sa')

        print("GCIS Compressing: ", f)
        if(os.path.isfile(compressed_file)):
            print(compressed_file, 'already exists, skipping')
        else:
            gcis.compress_gc_is(input_file, compressed_file)

        print("GCIS Decoding + SAIS: ", f)
        tt = []
        for i in range(number_of_samples):
            start_time = time.perf_counter()  # Start time
            gcis.decompress_saca(compressed_file, output_file)
            end_time = time.perf_counter()  # End time
            total_time = end_time - start_time
            print("Time GCIS =", "{:.2f}".format(total_time))
            tt.append(total_time)

        e.experiment_input.append(f)
        e.experiment_time.append("{:.2f}".format(statistics.mean(tt)))

    return e


def run_decode_saca_lcp(input_folder_path, output_folder_path):

    e = experiment_data()
    e.experiment_name = "GCIS-SACA-LCP"
    input_filename = os.listdir(input_folder_path)
    for f in input_filename:
        input_file = os.path.join(input_folder_path, f)
        compressed_file = os.path.join(output_folder_path, f+'.gcis')
        output_file = os.path.join(output_folder_path, f+'.gics.salcp')

        print("GCIS Compressing: ", f)
        if(os.path.isfile(compressed_file)):
            print(compressed_file, 'already exists, skipping')
        else:
            gcis.compress_gc_is(input_file, compressed_file)

        print("GCIS Decoding + SAIS + LCP: ", f)
        tt = []
        for i in range(number_of_samples):
            start_time = time.perf_counter()  # Start time
            gcis.decompress_saca(compressed_file, output_file)
            end_time = time.perf_counter()  # End time
            total_time = end_time - start_time
            print("Time GCIS + LCP =", "{:.2f}".format(total_time))
            tt.append(total_time)

        e.experiment_input.append(f)
        e.experiment_time.append("{:.2f}".format(statistics.mean(tt)))

    return e


def run_saca(input_folder_path, output_folder_path):

    # Create experiment folder
    print("Creating Experiment Folder if it not exists")
    os.makedirs(output_folder_path, exist_ok=True)

    # Run each variation
    print("Running Decode")
    decode_data = run_decode(input_folder_path, output_folder_path)
    print("Running Decode-Saca")
    decode_saca_data = run_decode_saca(input_folder_path, output_folder_path)
    print("Running Decode-Saca + LCP")
    decode_saca_lcp_data = run_decode_saca_lcp(
        input_folder_path, output_folder_path)
    print("Running Sais-Nong")
    saca_nong_data = run_sais_nong(input_folder_path, output_folder_path)
    print("Running Yuta Saca")
    saca_yuta_data = run_sais_yuta(input_folder_path, output_folder_path)
    print("Running Yuta Saca + LCP")
    saca_yuta_lcp_data = run_sais_lcp_yuta(
        input_folder_path, output_folder_path)
    print('Running Divsufsort Saca')
    saca_divsufsort_data = run_sais_divsufsort(
        input_folder_path, output_folder_path)
    print('Running Divsufsort Saca + LCP')
    saca_divsufsort_lcp_data = run_sais_divsufsort_lcp(
        input_folder_path, output_folder_path)
    print("Running decode + Nong Saca")
    decode_saca_nong_data = run_decode_nong(
        input_folder_path, output_folder_path)
    print("Running decode + Yuta Saca")
    decode_saca_yuta_data = run_decode_yuta(
        input_folder_path, output_folder_path)
    print("Running decode + Yuta Saca + LCP")
    decode_saca_yuta_lcp_data = run_decode_lcp_yuta(
        input_folder_path, output_folder_path)
    print('Running decode + Divsufsort Saca')
    decode_saca_divsufsort_data = run_decode_divsufsort(
        input_folder_path, output_folder_path)
    print('Running decode + Divsufsort Saca + LCP')
    decode_saca_divsufsort_lcp_data = run_decode_divsufsort(
        input_folder_path, output_folder_path)

    # Save data into csv file
    with open(os.path.join(output_folder_path, 'results-saca.csv'), 'w') as csvfile:
        writer = csv.writer(csvfile, delimiter=',')

        header = ['Input File',
                  decode_data.experiment_name,
                  decode_saca_data.experiment_name,
                  decode_saca_lcp_data.experiment_name,
                  saca_nong_data.experiment_name,
                  saca_yuta_data.experiment_name,
                  saca_divsufsort_data.experiment_name,
                  saca_yuta_lcp_data.experiment_name,
                  saca_divsufsort_lcp_data.experiment_name,
                  decode_saca_nong_data.experiment_name + ' Decompress',
                  decode_saca_nong_data.experiment_name + ' SACA',
                  decode_saca_nong_data.experiment_name,
                  decode_saca_yuta_data.experiment_name + ' Decompress',
                  decode_saca_yuta_data.experiment_name + ' SACA',
                  decode_saca_yuta_data.experiment_name,
                  decode_saca_divsufsort_data.experiment_name + ' Decompress',
                  decode_saca_divsufsort_data.experiment_name + ' SACA',
                  decode_saca_divsufsort_data.experiment_name,
                  decode_saca_yuta_lcp_data.experiment_name,
                  decode_saca_yuta_lcp_data.experiment_name + ' Decompress',
                  decode_saca_yuta_lcp_data.experiment_name + ' SACA',
                  decode_saca_divsufsort_lcp_data.experiment_name,
                  decode_saca_divsufsort_lcp_data.experiment_name + ' Decompress',
                  decode_saca_divsufsort_lcp_data.experiment_name + ' SACA']

        writer.writerow(header)

        time = list(zip(decode_data.experiment_input,
                        decode_data.experiment_time,
                        decode_saca_data.experiment_time,
                        decode_saca_lcp_data.experiment_time,
                        saca_nong_data.experiment_time,
                        saca_yuta_data.experiment_time,
                        saca_divsufsort_data.experiment_time,
                        saca_yuta_lcp_data.experiment_time,
                        saca_divsufsort_lcp_data,
                        decode_saca_nong_data.experiment_time,
                        decode_saca_yuta_data.experiment_time,
                        decode_saca_divsufsort_data.experiment_time),
                    decode_saca_yuta_lcp_data.experiment_time,
                    decode_saca_divsufsort_lcp_data.experiment_time)

        for i, t in enumerate(time):

            row = [decode_saca_data.experiment_input[i], t[1], t[2],
                   t[3], t[4], t[5], t[6], t[7], t[8], t[9][0], t[9][1], t[9][2], t[10][0], t[10][1], t[10][2],
                   t[11][0], t[11][1], t[11][2], t[12][0], t[12][1], t[12][2], t[13][0], t[13][1], t[13][2]]
            writer.writerow(row)


# argv[1] contains the path to the experiments folder
# argv[2] contains the path to the output folder
if __name__ == "__main__":
    if(len(sys.argv) < 3):
        print("Error")
        print("Usage: run-decode-saca.py <input_folder> <results_folder>")
    run_saca(sys.argv[1], sys.argv[2])
