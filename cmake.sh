#!/bin/bash
cmake -DCMAKE_INSTALL_PREFIX:PATH=./install \
      -DCMAKE_CXX_COMPILER:STRING=g++ \
      -DCMAKE_C_COMPILER:STRING=gcc \
      -DMPIEXEC_PREFLAGS:STRING=--oversubscribe \
      -DCMAKE_BUILD_TYPE:STRING=Release \
      -DEXAGOOP_ENABLE_MPI:BOOL=OFF \
      -DEXAGOOP_ENABLE_CUDA:BOOL=ON \
      -DAMReX_CUDA_ARCH=Volta \
      -DPYTHON_EXECUTABLE=$(which python3) \
      -DEXAGOOP_PRECISION:STRING=DOUBLE \
      -DAMReX_SPACEDIM=3 \
      ..
#make
cmake --build . --parallel $(sysctl -n hw.ncpu) #&> output.txt
#ctest -j $(sysctl -n hw.ncpu)
