#!/bin/sh -ex

rm -r \
	autom4te.cache aclocal.m4 configure config.h config.log config.status \
	Makefile.in src/Makefile.in src/.deps/ stamp-h1 install-sh missing depcomp \
	|| true >& /dev/null

if [[ "$1" == "--clean" ]]
then
	exit 0
fi

aclocal
automake --add-missing
autoconf

./configure
make

if [[ "$1" == "--distrib" ]]
then
	make distclean
	echo "OK for distribution"
fi
