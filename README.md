# A Grammar Compression Algorithm by the Induced Suffix Sorting Framework

## Introduction

The GCIS algorithm is capable of computing a context-free grammar by using the suffix-sorting framework by Nong *et al*. Moreover, by using GCIS is possible to extract substrings from the compressed text without decompressing it and to compute the Suffix and LCP arrays as a by-product of the decompression. 

## Compilation

First download the GCIS repository

```shell
git clone https://github.com/danielsaad/GCIS.git
```

In order to compile the program, it is necessary to download, compile and install the **Succinct Data Structure Library (SDSL)**.

```shell
git clone https://github.com/simongog/sdsl-lite
cd sdsl-lite
./install-sh <path to GCIS repository>
```

Now, it is possible to compile the library and the binaries

```shell
  mkdir build
  cd build
  cmake .. -DCMAKE_BUILD_TYPE=RELEASE
  make -j
  make install
```

Alternatively, one can install through the build.sh file

```shell
    ./build.sh
```

The binaries are available in the **bin** folder and the static library in the  **lib** folder.

## Basic Usage

### Choosing the CODEC

There are two CODECs which can be used on the GCIS algorithm:
- Simple8b and;
- Elias-Fano.

They can be chosen specifying in the command line the flags `-s8b` and `-ef` respectively.

Simple8b has some limitations, despite being able do encode and decode, it  cannot neither compute the SA and LCP arrays from decompression nor extract substrings.


### Compression

To compress a single text, one must execute:

```bash
./gc-is-codec -c <input_file> <compressed_file> <CODEC Flag>
```

where `<input_file>` and `<compressed_file>` correspond, respectively, to the paths of the input and compressed files.

### Decompression

To decompress a previously GCIS compressed file, one must execute:


```bash
./gc-is-codec -d <compressed_file> <decompressed_file> <CODEC Flag>
```

where `<compressed_file>` and `<decompressed_file>` stands for the paths of the compressed GCIS file and the decompressed text file respectively.

### Extraction

To extract substrings directly from the compressed GCIS file, it is necessary to run:

```bash
./gc-is-codec -e <compressed_file> <queries_file> <CODEC Flag>
```

where `<compressed_file>` stands for the compressed GCIS file path and `<query_file>` stands for the queries file path. The `<queries_file>` is a simple text file containing on each line the positions `[l,r]` of the substring `T[l,r]` from the original text to be extracted.

### SA and LCP arrays construction

To compute the suffix array from the compressed text, it is necessary to run:

```bash
./gc-is-codec -s <compressed_file> <suffix_array_file> <CODEC Flag>
```
where `<compressed_file>` stands for the GCIS compressed file path and `<suffix_array_file>` corresponds to the binary file where the suffix array contents shall be stored. 

If it is desired to compute the LCP array as well, one can simply use the previous command and modify the `-s` flag to `-l`:

```bash
./gc-is-codec -l <compressed_file> <suffix_array_file> <CODEC Flag>
```


## API

The GCIS library has the following types and functions

```cpp
//Type of the gc-is dictionary
/***
 * lcp_coder can be gcis_s8b_codec (simple8b encoding/decoding)
 * or gcis_eliasfano_codec (Elias-Fano encoding/decoding)
 ***/
gc_is_dictionary<lcp_coder>;

// Grammar compress the str by using the induced suffix sorting framework
void gc_is_dictionary<lcp_coder>::encode(char* str);

// Grammar decompress the dictionary into the original string.
char* gc_is_dictionary<lcp_coder>::decode(char* str);

// Extract substrings T[l,r] described by a pairs [l,r]
void gc_is_dictionary<lcp_coder>::extract_batch(vector<pair<int,int>>);

// Compute SA and LCP arrays and returns decoded text
char* gc_is_dictionary<lcp_coder>::decode_saca(uint_t** SA, int_t **LCP);

// Compute SA and LCP arrays and returns decoded text
char* gc_is_dictionary<lcp_coder>::decode_saca_lcp(uint_t** SA, int_t **LCP);

// Outputs the total number of bytes to represent the grammar.
uint64_t gc_is_dictionary<lcp_coder>::size_in_bytes();

// Serialize compressor data into a binary file
void gc_is_dictionary<lcp_coder>::serialize(std::ostream&);

// Load the binary file containing the compressor into main memory
void gc_is_dictionary<lcp_coder>::load(std::istream &i) = 0;


```
