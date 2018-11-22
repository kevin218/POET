#!/bin/bash
#this will remove any partial builds, if there was a failure

if [ "$1" == "clean" ]
then
    rm ext_func/*.so
    exit
fi

# python2.7 setup.py build_ext --inplace
python setup.py build_ext --inplace
rename_so_files.sh
mv -f *.so ext_func/
rm -rf build/
echo "~~You made it!~~"
