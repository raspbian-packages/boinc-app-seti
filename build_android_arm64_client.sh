#! /bin/sh

#
# See: http://boinc.berkeley.edu/trac/wiki/AndroidBuildApp#
#

# Script to compile various BOINC libraries for Android to be used
# by science applications

COMPILEBOINC="yes"
CONFIGURE="yes"
MAKECLEAN="yes"
export HOME=/mydisks/b/usr/

export BOINC="../boinc" #BOINC source code

export ANDROIDTC="$HOME/androidarm64-tc"
export TCBINARIES="$ANDROIDTC/bin"
export TCINCLUDES="$ANDROIDTC/aarch64-linux-android"
export TCSYSROOT="$ANDROIDTC/sysroot"
export STDCPPTC="$TCINCLUDES/lib/libstdc++.a"

export PATH="$TCBINARIES:$TCINCLUDES/bin:$PATH"
export CC=aarch64-linux-android-gcc
export CXX=aarch64-linux-android-g++
export LD=aarch64-linux-android-ld
export CFLAGS="--sysroot=$TCSYSROOT -DANDROID -DDECLARE_TIMEZONE -Wall -I$TCINCLUDES/include -O3 -fomit-frame-pointer -funsafe-math-optimizations -freciprocal-math -fPIE"
export CXXFLAGS="${CFLAGS}"
export LDFLAGS="-static-libstdc++ -static-libgcc -L$TCSYSROOT/usr/lib -L$TCINCLUDES/lib -llog -fPIE -pie"
export GDB_CFLAGS="--sysroot=$TCSYSROOT -Wall -g -I$TCINCLUDES/include"
export PKG_CONFIG_SYSROOT_DIR=$TCSYSROOT

# Prepare android toolchain and environment
#./build_androidtc_arm.sh


echo "==================building seti_boinc=========================="
if [ -n "$MAKECLEAN" ]; then
make clean
fi
if [ -n "$CONFIGURE" ]; then
./configure --host=aarch64-linux-android --prefix=$TCINCLUDES --disable-graphics --disable-server --with-boinc-platform=aarch64-android-linux-gnu --with-ssl=${TCINCLUDES}
fi
make 
make install

echo "=============================seti_boinc done============================="

