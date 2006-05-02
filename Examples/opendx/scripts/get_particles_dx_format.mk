#!/bin/bash
# from young to old
#This script reads through ellipsis3D binary output particle files, and #converts them into something useful (ie. opendx format for plotting).
#You will need to change this (eg the numbering for a start). Requires #Binary_read.pl and getmaterialpoints.pl to work properly. Watch the naming #conventions in these to files too...


rm materialOnePoints_0* materials_0*
for i in `ls ../../../iso-test.*.particles | perl -p -e 's/\.\.\/\.\.\/\.\.\/iso-test\.//' | perl -p -e 's/\.particles//'`
do

perl Binary_read.pl -f ../../../iso-test.${i}.particles -F materials_${i}.dat
perl getmaterialpoints.pl ${i} 

sed -e 1s/materialOnePoints_NUMBER/materialOnePoints_${i}/ materialOnePoints.general > tmp

numberlines=`nl materialOnePoints_${i}.dat | tail -n1 |awk '{ print $1}'`
if [ "${numberlines}" ]
then
   echo ${numberlines}
else
   numberlines=0
fi
sed -e 2s/LINES/${numberlines}/ tmp > materialOnePoints_${i}.general
rm tmp


done
