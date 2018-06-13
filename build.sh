#!/usr/bin/env bash
mkdir -p external/sdsl-lite
git clone https://github.com/simongog/sdsl-lite.git external/sdsl-lite
./external/sdsl-lite/install.sh .
cd external/repair
 # compile repair
mkdir -p build
cd build
cmake ..
make install
cd ../../../ # go to root folder
#compile sais
cd external/sais-2.4.1
mkdir -p build
cmake .. -DCMAKE_BUILD_TYPE=RELEASE -DCMAKE_INSTALL_PREFIX=../../../
make -j
make install
cd ../../../
# compile gcis
mkdir -p build
cd build
cmake .. -DCMAKE_BUILD_TYPE=RELEASE
make -j
make install
cd ..

# copy binaries in other projects to bin folder
cp external/repair/bin/despair external/repair/bin/despair-memory external/repair/bin/repair external/repair/bin/repair-memory bin/
