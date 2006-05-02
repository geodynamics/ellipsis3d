#!/usr/bin/perl -w

# Use appropriate maths libraries
use Math::Trig;

$time=0;
$time=$ARGV[0];

print "$time \n";

open(INPUT,"<materials_$time.dat") || die "Could not open input file\n";

# SKIP information from (4) header lines (line12)
    local(@header);
    local($i);

    for($i=1;$i<3;$i++) {
        $_=<INPUT>; s/\#// ; $header[$i] = $_; # header line
    }

#This reads each line of input file, and assigns fields to variables
local($index)=0;
while (<INPUT>)
  {
    local(@this_line)=split(' ');
    $x[$index]=$this_line[0]; #this is second field
    $z[$index]=$this_line[1];
    $y[$index]=$this_line[2];
    $t[$index]=$this_line[3];
     $index++;
  }


 open(OUTPUT,">>materialOnePoints_$time.dat") || die "Could not open output file\n";
 foreach $j (0..$#x) {
    if ($t[$j]>0.5) {
        print OUTPUT "$x[$j] $z[$j] $y[$j] $t[$j] \n"; 
    }
  }

    close(INPUT);
    close(OUTPUT);

# And done
