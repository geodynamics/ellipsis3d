#!/bin/ksh
# from young to old
#This script reads through ellipsis3D binary output particle files, and #converts them into something useful (ie. opendx format for plotting).
#You will need to change this (eg the numbering for a start). Requires #Binary_read.pl and getmaterialpoints.pl to work properly. Watch the naming #conventions in these to files too...


rm materialOnePoints_0* materials_0*
for i in 00001 00051 00097 00151 00198 00252 00298 00352 00401 00449 00500 00551 00598 00649 00698 00751 00801 00849 00897 00949 00997 01051 01099 01151 01200 01247 01303 01348 01400 01450 01498 01547 01600 01651 01701 01751 01803 01851 01898 01949 01996   
do

perl Binary_read.pl -f ../iso-test.${i}.particles -F materials_${i}.dat
perl getmaterialpoints.pl ${i} 

sed -e 1s/materialOnePoints/materialOnePoints_${i}/ materialOnePoints.general > tmp

numberlines=`nl materialOnePoints_${i}.dat | tail -n1 |awk '{ print $1}'`
if [ "${numberlines}" ]
then
   echo ${numberlines}
else
   numberlines=0
fi
sed -e 2s/112/${numberlines}/ tmp > materialOnePoints_${i}.general
rm tmp


done
