# Build roadrunner and FEBio
Setup information for cluster/workstation to run FEBio simulations with roadrunner.

https://libroadrunner.readthedocs.io/en/latest/Installation/installation.html

## TODO
- [ ] install and use newer gcc version for compilation 11.2


## requirements
- cmake => 3.18
- gcc => 9.4
- Intel MKL (Math Kernel library), https://www.intel.com/content/www/us/en/developer/tools/oneapi/onemkl.html; intel-mkl:amd64 (2020.0.166-1)

```bash
# latest gcc version
# sudo add-apt-repository ppa:ubuntu-toolchain-r/ppa -y
# sudo apt-get update
# sudo apt-get install g++-10 gcc-10


# latest cmake version
wget -O - https://apt.kitware.com/keys/kitware-archive-latest.asc 2>/dev/null | gpg --dearmor - | sudo tee /etc/apt/trusted.gpg.d/kitware.gpg >/dev/null
sudo apt-add-repository "deb https://apt.kitware.com/ubuntu/ $(lsb_release -cs) main"
sudo apt update
sudo apt-get -y install cmake cmake-gui cmake-curses-gui

# Intel MKL
sudo apt-get -y install intel-mkl
```

## Setup build directory
```bash
mkdir buildroadrunner
cd buildroadrunner
```

## Build LLVM
```bash
git clone https://github.com/sys-bio/llvm-13.x.git
cd llvm-13.x/llvm
mkdir build
cd build

cmake -DCMAKE_INSTALL_PREFIX="../../../llvm-13.x-release" -DCMAKE_BUILD_TYPE="Release" ..
cmake --build . --target install --config Release
cd ../../../
```

## Roadrunner dependencies
```bash
git clone https://github.com/sys-bio/libroadrunner-deps.git --recurse-submodules
cd libroadrunner-deps
git switch release
mkdir build
cd build
rm -rf ../../libroadrunner-deps-release

cmake -DCMAKE_INSTALL_PREFIX="../../libroadrunner-deps-release" -DCMAKE_BUILD_TYPE="Release" ..
# cmake -DCMAKE_INSTALL_PREFIX="../../libroadrunner-deps-release" -DCMAKE_INSTALL_LIBDIR=lib -DCMAKE_BUILD_TYPE="Release" -DCMAKE_CXX_COMPILER="/opt/hlrs/non-spack/2022-02/compiler/gcc/11.2.0/bin/c++" -DCMAKE_C_COMPILER="/opt/hlrs/non-spack/2022-02/compiler/gcc/11.2.0/bin/gcc" ..

cmake --build . --parallel 10 --target install --config Release
cd ../../
``` 

## Roadrunner
```bash
git clone https://github.com/sys-bio/roadrunner.git

cd roadrunner
git switch release
mkdir build
cd build

cmake -DCMAKE_INSTALL_PREFIX="../../roadrunner-release" -DLLVM_INSTALL_PREFIX="../../llvm-13.x-release" -DRR_DEPENDENCIES_INSTALL_PREFIX="../../libroadrunner-deps-release" -DCMAKE_BUILD_TYPE="Release" ..

cmake --build . --parallel 10 --target install --config Release

cd ../../../
``` 

## FEBio
```bash
git clone https://github.com/febiosoftware/FEBio.git
cd FEBio
git checkout 6f1a2ebcd5dce7c5fc4c9f5bf60e3a9acccdd8d1

## In ./FEBioMix/FESBMPointSource.cpp ll.360 change the following:
double min_d = 0;//std::numeric_limits<double>::max();

## In ./FECore/DumpStream.h ll.68-71 change the following:
//#ifndef uchar
//#define uchar unsigned char
//#endif
typedef unsigned char uchar;

mkdir cbuild
cd cbuild

cmake ..

cmake -DUSE_MKL=ON -DMKLROOT="/usr/include/mkl" -DMKL_INC="/usr/include/mkl" -DMKL_LIB_DIR="/usr/lib/x86_64-linux-gnu" -DMKL_OMP_LIB="/usr/lib/x86_64-linux-gnu/libiomp5.so" ..

# cmake -DMKLROOT="/opt/hlrs/non-spack/2022-02/compiler/intel/2022.1.2_oneapi/mkl" -DMKL_INC="/opt/hlrs/non-spack/2022-02/compiler/intel/2022.1.2_oneapi/mkl/latest/include" -DMKL_LIB_DIR="/opt/hlrs/non-spack/2022-02/compiler/intel/2022.1.2_oneapi/mkl/latest/lib/intel64" -DMKL_OMP_LIB="/opt/hlrs/non-spack/2022-02/compiler/intel/2022.1.2_oneapi/compiler/latest/linux/compiler/lib/intel64_lin/libiomp5.so" -DCMAKE_CXX_COMPILER="/opt/hlrs/non-spack/2022-02/compiler/gcc/11.2.0/bin/c++" -DCMAKE_C_COMPILER="/opt/hlrs/non-spack/2022-02/compiler/gcc/11.2.0/bin/gcc" ..
make -j 10

cd ../../
```

## FEBio plugins
https://github.tik.uni-stuttgart.de/login
```bash

git clone https@github.tik.uni-stuttgart.de:isd/FEBioTPM.git
git clone git@github.tik.uni-stuttgart.de:isd/FEBioTPM_IRI.git

cd FEBioTPM

# CMakeLists.txt -> adapt paths
# set(ROADRUNNER_INSTALL_PREFIX "/home/mkoenig/git/porous_media/buildroadrunner/roadrunner-release")
# set(LLVM_INSTALL_PREFIX "/home/mkoenig/git/porous_media/buildroadrunner/llvm-13.x-release")

rm -rf build
rm CMakeCache.txt
mkdir build
cd build

cmake -DFEBio_SDK="../../FEBio" ..
# cmake -DFEBio_SDK="/zhome/academic/HLRS/isd/isdluis/FEBio" -DCMAKE_CXX_COMPILER="/opt/hlrs/non-spack/2022-02/compiler/gcc/11.2.0/bin/c++" -DCMAKE_C_COMPILER="/opt/hlrs/non-spack/2022-02/compiler/gcc/11.2.0/bin/gcc" ..
make -j 10

# setup plugin path

## Change FEBio XML: /FEBio/cbuild/bin/febio.xml
<febio_config version="3.0">
    <default_linear_solver type="pardiso"/>
    <import>/home/mkoenig/git/porous_media/buildroadrunner/FEBioTPM/build/lib/libFEBioTPM.so</import>
</febio_config>
```