#!/bin/ksh
# from young to old
# Example script to convert Ellipsis3D node data files into DX format.
# Will need to modify to use.


rm nodes0*
for i in 00001 00051 00097 00151 00198 00252 00298 00352 00401 00449 00500 00551 00598 00649 00698 00751 00801 00849 00897 00949 00997 01051 01099 01151 01200 01247 01303 01348 01400 01450 01498 01547 01600 01651 01701 01751 01803 01851 01898 01949 01996 
do

# cp nodes.general nodes${i}.general

sed -e 1s/00001/${i}/ nodes.general > nodes${i}.general

perl processing_data.pl ${i}

done
