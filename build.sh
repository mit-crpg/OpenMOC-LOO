#!/bin/bash -ex

CFLAGS="-Wall -Werror -std=gnu++0x `pkg-config Magick++ --cflags` -fopenmp -O3 -march=native"
LDFLAGS="-lexpat -lsilo `pkg-config Magick++ --libs`"

if [ "$1" = '--release' ]
then
    CFLAGS="$CFLAGS"
elif [ "$1" = "--debug" ]
then
    CFLAGS="-Wall -Werror -std=gnu++0x `pkg-config Magick++ --cflags` -fopenmp -g"
elif [ "$1" = '--profile' ]
then
    CFLAGS="$CFLAGS -pg -g"
elif [ "$1" = "--gprofile" ]
then
    CFLAGS="$CFLAGS -lprofiler -g"
fi

if ! test -d bin
then
	mkdir bin
fi
gcc src/*.cpp -o bin/openmoc $CFLAGS $LDFLAGS

if [ "$1" == "--profile" ]
then
    rm gmon.out || true
    ./bin/openmoc > /dev/null
    gprof ./bin/openmoc
    rm gmon.out || true
elif [ "$1" == "--gprofile" ]
then
    export CPUPROFILE_FREQUENCY=1000
    export CPUPROFILE=gprofile.prof
    ./bin/openmoc > /dev/null
    pprof --gv ./bin/openmoc gprofile.prof
fi