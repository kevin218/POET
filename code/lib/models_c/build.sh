#!/bin/bash
#this will remove any partial builds, if there was a failure

if [ "$1" == "clean" ]
then
    rm ext_func/*.so
    exit
fi

if [ "$1" == "python2" ]
then
    python2.7 setup.py build_ext --inplace
else
    python setup.py build_ext --inplace
fi

rename_so_files.sh
mv -f *.so ext_func/
rm -rf build/
echo "~~You made it!~~"
