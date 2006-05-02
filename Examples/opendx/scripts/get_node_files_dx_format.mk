#!/bin/bash
# from young to old
# Example script to convert Ellipsis3D node data files into DX format.
# Will need to modify to use.


rm nodes0*
for i in `ls ../../../iso-test.*.node_data | perl -p -e 's/\.\.\/\.\.\/\.\.\/iso-test\.//' | perl -p -e 's/\.node_data//'`
do


perl processing_data.pl ${i}

sed -e 1s/nodesNUMBER/nodes${i}/ nodes.general > tmp

numberlines=`nl nodes${i}.dat | tail -n1 |awk '{ print $1}'`
if [ "${numberlines}" ]
then
   echo ${numberlines}
else
   numberlines=0
fi
sed -e 2s/LINES/${numberlines}/ tmp > nodes${i}.general
rm tmp

done
