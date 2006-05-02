#!/usr/bin/perl -w

# Use appropriate maths libraries
use Math::Trig;

$time=0;

$time=$ARGV[0];




# Script to convert ellipsis to dx format. Node specified NAMES in script - you'll need to change these to use 'em.
#This is just meant to help.
# This is only meant to show you how it can be done, it is not in final working order
# and could change for different runs
 print "$time \n";
$name="../../../iso-test.$time.node_data";    #Note ";" at end of EVERY line in perl script
print "Trying to open $name \n";

#Now open input file

open(INPUT,"<$name") || die "Could not open input file\n";

# SKIP information from (4) header lines (line12)
    local(@header);
    local($i);

    for($i=1;$i<5;$i++) {
	$_=<INPUT>; s/\#// ; $header[$i] = $_; # header line
    }

   
#This reads each line of input file, and assigns fields to variables
local($index)=0;
while (<INPUT>) 
  {
    local(@this_line)=split(' ');
    $x[$index]=$this_line[1]; #this is second field
    $z[$index]=$this_line[2];
    $y[$index]=$this_line[3];
    $t[$index]=$this_line[7];
    $vx[$index]=$this_line[4]; #this is second field
    $vz[$index]=$this_line[5];
    $vy[$index]=$this_line[6];

     $index++;
  }
   
   
 open(OUTPUT,">>nodes$time.dat") || die "Could not open output file\n";  
 foreach $j (0..$#x) {
    print OUTPUT "$x[$j] $z[$j] $y[$j] $vx[$j] $vz[$j] $vy[$j] $t[$j] \n"; #Better to have an array and cycle through these -just a thought 
  }
# Note: the \n at the end of the print statement is the "end of line" symbol.

    close(INPUT);
    close(OUTPUT);

# And done
