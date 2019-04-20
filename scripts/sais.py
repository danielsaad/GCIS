import sys
import os
import subprocess
import re


# Run SAIS YUTA
def sais_yuta(input_file, output_file):
    p = subprocess.Popen(['../bin/sais-yuta', input_file, output_file],
                         stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    return p.communicate()

# Run SAIS YUTA + LCP
def sais_lcp_yuta(input_file, output_file):
    p = subprocess.Popen(['../bin/sais-lcp-yuta', input_file, output_file],
                         stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    return p.communicate()

# Run SAIS NONG


def sais_nong(input_file, output_file):
    p = subprocess.Popen(['../bin/sais-nong', input_file, output_file],
                         stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    return p.communicate()


# Run Divsufsort
def sais_divsufsort(input_file, output_file):
    p = subprocess.Popen(['../bin/sais-divsufsort', input_file, output_file],
                         stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    return p.communicate()


# Run Divsufsort + LCP
def sais_divsufsort_lcp(input_file, output_file):
    p = subprocess.Popen(['../bin/sais-divsufsort-lcp', input_file, output_file],
                         stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    return p.communicate()


'''
    Run GCIS to decode compressed file and run SAIS YUTA 
    Returns a pair (decode_time,sais_yuta_time)
'''
def decode_sais_yuta(compressed_file, output_file):
    p = subprocess.Popen(['../bin/decode-sais-yuta', compressed_file,
                          output_file], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    out, _ = p.communicate()
    lines = out.decode("utf-8").split('\n')
    for l in lines:
        if(l.startswith('Decompression Time:')):
            decompress_time = l.split()[2]
        elif(l.startswith('Suffix Array Construction Time:')):
            saca_time = l.split()[4]

    return (decompress_time, saca_time)


'''
    Run GCIS to decode compressed file and run SAIS YUTA + LCP 
    Returns a pair (decode_time,sais_yuta_time)
'''
def decode_sais_lcp_yuta(compressed_file, output_file):
    p = subprocess.Popen(['../bin/decode-sais-lcp-yuta', compressed_file,
                          output_file], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    out, _ = p.communicate()
    lines = out.decode("utf-8").split('\n')
    for l in lines:
        if(l.startswith('Decompression Time:')):
            decompress_time = l.split()[2]
        elif(l.startswith('Suffix Array Construction Time:')):
            saca_time = l.split()[4]

    return (decompress_time, saca_time)



''' 
    Run GCIS to decode compressed file and run SAIS NONG 
    Returns a pair (decode_time,sais_nong_time)
'''


def decode_sais_nong(compressed_file, output_file):
    p = subprocess.Popen(['../bin/decode-sais-nong', compressed_file,
                          output_file], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    out, _ = p.communicate()
    lines = out.decode("utf-8").split('\n')
    for l in lines:
        if(l.startswith('Decompression Time:')):
            decompress_time = l.split()[2]
        elif(l.startswith('Suffix Array Construction Time:')):
            saca_time = l.split()[4]

    return (decompress_time, saca_time)


'''
    Run GCIS to decode compressed file and run SAIS Divsufsort 
    Returns a pair (decode_time,sais_divsufsort_time)
'''
def decode_sais_divsufsort(compressed_file, output_file):
    p = subprocess.Popen(['../bin/decode-sais-divsufsort', compressed_file,
                          output_file], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    out, _ = p.communicate()
    lines = out.decode("utf-8").split('\n')
    for l in lines:
        if(l.startswith('Decompression Time:')):
            decompress_time = l.split()[2]
        elif(l.startswith('Suffix Array Construction Time:')):
            saca_time = l.split()[4]

    return (decompress_time, saca_time)


'''
    Run GCIS to decode compressed file and run SAIS Divsufsort 
    Returns a pair (decode_time,sais_divsufsort_time)
'''
def decode_sais_divsufsort_lcp(compressed_file, output_file):
    p = subprocess.Popen(['../bin/decode-sais-divsufsort-lcp', compressed_file,
                          output_file], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    out, _ = p.communicate()
    lines = out.decode("utf-8").split('\n')
    for l in lines:
        if(l.startswith('Decompression Time:')):
            decompress_time = l.split()[2]
        elif(l.startswith('Suffix Array Construction Time:')):
            saca_time = l.split()[4]

    return (decompress_time, saca_time)
