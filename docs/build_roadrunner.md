# Build roadrunner and FEBio
Setup information for cluster/workstation to run FEBio simulations with roadrunner.

## TODO
- [ ] install and use newer gcc version for compilation 11.2


## requirements

```bash
sudo apt-get -y install cmake cmake-gui
```

## Setup build directory
```
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
cmake --build . --parallel 20 --target install --config Release
cd ../../../
```

## Roadrunner dependencies
```bash
git clone https://github.com/sys-bio/libroadrunner-deps.git --recurse-submodules
cd libroadrunner-deps
git switch release
mkdir build
cd build

cmake -DCMAKE_INSTALL_PREFIX="../../libroadrunner-deps-release" -DCMAKE_INSTALL_LIBDIR=lib -DCMAKE_BUILD_TYPE="Release" ..
# cmake -DCMAKE_INSTALL_PREFIX="../../libroadrunner-deps-release" -DCMAKE_INSTALL_LIBDIR=lib -DCMAKE_BUILD_TYPE="Release" -DCMAKE_CXX_COMPILER="/opt/hlrs/non-spack/2022-02/compiler/gcc/11.2.0/bin/c++" -DCMAKE_C_COMPILER="/opt/hlrs/non-spack/2022-02/compiler/gcc/11.2.0/bin/gcc" ..

cmake --build . --parallel 10 --target install --config Release 

cd ../../
``` 

## Roadrunner
```swith
git clone https://github.com/sys-bio/roadrunner.git

cd roadrunner
git switch release
rm -rf build
mkdir build
cd build

### Changes in ImportRoadrunnerAndDepedencies.cmake
### Auskommentieren 
# find_package(iconv CONFIG REQUIRED)
# find_package(LibXml2 CONFIG REQUIRED)
# SUNDIALS groß

## hinzufügen 
# find_package(expat CONFIG REQUIRED)
## sbml statt libsbml
# find_package(sbml-static CONFIG REQUIRED)


cmake -DCMAKE_INSTALL_PREFIX="../../roadrunner-release" -DCMAKE_INSTALL_LIBDIR=lib  -DLLVM_INSTALL_PREFIX="../../llvm-13.x-release" -DRR_DEPENDENCIES_INSTALL_PREFIX="../../libroadrunner-deps-release" -DCMAKE_BUILD_TYPE="Release" -DCMAKE_CXX_COMPILER="/opt/hlrs/non-spack/2022-02/compiler/gcc/11.2.0/bin/c++" -DCMAKE_C_COMPILER="/opt/hlrs/non-spack/2022-02/compiler/gcc/11.2.0/bin/gcc" ..
cmake --build . --parallel 10 --target install --config Release

cd ../../../

git clone https://github.com/febiosoftware/FEBio.git

cd FEBio

git checkout 6f1a2ebcd5dce7c5fc4c9f5bf60e3a9acccdd8d1

## In ./FEBioMix/FESBMPointSource.cpp ll.460 folgendes ändern:
double min_d = 0;//std::numeric_limits<double>::max();

## In ./FECore/DumpStream.h ll.68-71 folgendes ändern:
//#ifndef uchar
//#define uchar unsigned char
//#endif
typedef unsigned char uchar;

##

mkdir cbuild
cd cbuild

cmake -DMKLROOT="/opt/hlrs/non-spack/2022-02/compiler/intel/2022.1.2_oneapi/mkl" -DMKL_INC="/opt/hlrs/non-spack/2022-02/compiler/intel/2022.1.2_oneapi/mkl/latest/include" -DMKL_LIB_DIR="/opt/hlrs/non-spack/2022-02/compiler/intel/2022.1.2_oneapi/mkl/latest/lib/intel64" -DMKL_OMP_LIB="/opt/hlrs/non-spack/2022-02/compiler/intel/2022.1.2_oneapi/compiler/latest/linux/compiler/lib/intel64_lin/libiomp5.so" -DCMAKE_CXX_COMPILER="/opt/hlrs/non-spack/2022-02/compiler/gcc/11.2.0/bin/c++" -DCMAKE_C_COMPILER="/opt/hlrs/non-spack/2022-02/compiler/gcc/11.2.0/bin/gcc" ..
make

cd ../../

git clone git@github.tik.uni-stuttgart.de:isd/FEBioTPM_IRI.git

cd FEBioTPM_IRI

# CMakeLists.txt -> Pfade anpassen
# set(ROADRUNNER_INSTALL_PREFIX "/zhome/academic/HLRS/isd/isdluis/buildroadrunner/roadrunner-release")
# set(LLVM_INSTALL_PREFIX "/zhome/academic/HLRS/isd/isdluis/buildroadrunner/llvm-13.x-release")
# set(CMAKE_MODULE_PATH "${ROADRUNNER_INSTALL_PREFIX}/lib/cmake")


rm -rf build
mkdir build
cd build

cmake -DFEBio_SDK="/zhome/academic/HLRS/isd/isdluis/FEBio" -DCMAKE_CXX_COMPILER="/opt/hlrs/non-spack/2022-02/compiler/gcc/11.2.0/bin/c++" -DCMAKE_C_COMPILER="/opt/hlrs/non-spack/2022-02/compiler/gcc/11.2.0/bin/gcc" ..
make

## Change FEBio XML: /FEBio/cbuild/bin/febio.xml
<febio_config version="3.0">
    <default_linear_solver type="pardiso"/>
    <import>/zhome/academic/HLRS/isd/isdluis/FEBioTPM_IRI/build/lib/libFEBioTPM.so</import>
</febio_config>