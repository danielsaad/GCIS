# A Grammar Compression Algorithm by the Induced Suffix Sorting Framework

# Introduction

This library is able to compress an input and generate a grammar based on
the induced suffix sorting of (Nong _et al_. 2009). We call this algorithm GCIS.

# Compilation

First download the GCIS repository

```shell
git clone https://github.com/danielsaad/GCIS.git
```

In order to compile the program, it is necessary to download, compile and install the **Succinct Data Structure Library**.

```shell
git clone https://github.com/simongog/sdsl-lite
cd sdsl-lite
./install-sh <path to GC-IS repository>
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

# API

The GC-IS library has the folowing types and functions

```cpp
//Type of the gc-is dictionary
gc_is_dictionary<lcp_coder>;
//lcp_coder can be gcis_unary_codec (unary coding) or gcis_s8b_codec (simple8b coding)
// Grammar compress the str by using the induced suffix sorting framework
gc_is_dictionary<lcp_coder>::encode(char* str);

// Grammar decompress the dictionary into the original string.
gc_is_dictionary<lcp_coder>::decode(char* str);

// Grammar decompress the dictionary into the original string without decompressing the dictionary (Very slow). */
gc_is_dictionary<lcp_coder>::inline_decode(char* str);

// Outputs the total number of bytes to represent the grammar.
gc_is_dictionary<lcp_coder>::size_in_bytes();
```
