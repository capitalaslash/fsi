#! /bin/bash

# default to dbg
_METHOD=dbg

# use env variable to select method
if [[ $METHOD != "" ]]; then
  _METHOD=$METHOD
fi

# set method specific flags
cmake_flags=""
if [[ $_METHOD == "opt" ]]; then
  cmake_flags+="-DCMAKE_BUILD_TYPE:STRING=Release"
elif [[ $_METHOD == "dbg" ]]; then
  cmake_flags+="-DCMAKE_BUILD_TYPE:STRING=Debug"
fi

SRC_DIR=$PWD/..
BUILD_DIR=$PWD/../../libmesh-app_build-$_METHOD
INSTALL_DIR=$PWD/../../libmesh-app_build-$_METHOD/install

mkdir -p $BUILD_DIR
cd $BUILD_DIR
rm -rf CMake*

cmake \
  $SRC_DIR \
  $cmake_flags \
  -DCMAKE_INSTALL_PREFIX:PATH=$INSTALL_DIR \
  -DLIBMESH_DIR=/data/aslash/software/libmesh-git/install \
  $*

