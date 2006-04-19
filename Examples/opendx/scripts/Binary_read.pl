#! /usr/local/bin/perl

require "getopts.pl" ;

#
# Perl script to take particle data and
# plot using (in this case) GMT to 
# produce a postscript file of specified size.
# 
# Assumption is that this is a frame for a movie
# and hence that time information is meaningful
#


&Getopts('F:f:');

# Options:   -f: Filename for input data
# Options:   -F: Filename (root) for output data


# default values for parameters if not specified

if($opt_F eq "") {
  $opt_F = "ascii-conversion";
}

# Read the particle file ...   Assume 2D !! 

open(PAR,"< $opt_f") || die "File not found $opt_f\n";
open(OUT,"> $opt_F") || die "File not found $opt_F\n";
# open(OUT,">$name") || die "Cannot open file $name : $!\n";

# header information

# line 1
@oneline = split(' ',<PAR>);
@oneline[1] =~s/FORMAT=//;
if(@oneline[1] ne "binary") {
  die "Binary format input file required\n";
}

# line 2

@oneline = split(' ',<PAR>);
@oneline[1] =~s/PARTICLES=//;
@oneline[2] =~s/DIMENSIONS=//;
@oneline[4] =~s/DATA_SIZE=//;

$particles = @oneline[1];
$dim = @oneline[2];
$datasize = @oneline[4];

#line 3

@oneline = split(' ',<PAR>);
@oneline[1] =~s/DATA_ELEMENTS=//;
$datacolumns=@oneline[1];
$alldatacolumns=$datacolumns+5; #Assuming dims=3

#line 4

@oneline = split(' ',<PAR>);

#line 5-5+alldatacolumns !

for($i=1;$i<=$alldatacolumns;$i++) {

@oneline = split(' ',<PAR>);
@oneline[1] =~ s/=+\d// ;
$dataname[$i] = @oneline[1];

}

# Last line of ascii - the ^L separator - discard this line

<PAR>;

# 22/7  

read(PAR,$buffer,$datasize);
print "Datasize $datasize \n";
@test = unpack("f",$buffer);
$twentytwooverseven=@test[0];
print "test @test[0] \n";
$diff = $twentytwooverseven - 3.1;

if($diff < 0.0) {
  $diff *= -1.0;
}

if($diff > 0.1) {
  die "Binary file may be in the wrong format !!\n";
}


# Now we have blocks of binary data

# Material properties: signed integer

read(PAR,$buffer,$datasize*(1+$particles));
@material = unpack("i*",$buffer);

# X and Z coordinates: float

read(PAR,$buffer,$datasize*(1+$particles));
@Xcoord = unpack("f*",$buffer);

read(PAR,$buffer,$datasize*(1+$particles));
@Zcoord = unpack("f*",$buffer);

read(PAR,$buffer,$datasize*(1+$particles));
@Ycoord = unpack("f*",$buffer);

# Particle weight: float
read(PAR,$buffer,$datasize*(1+$particles));
@weight = unpack("f*",$buffer);

# Output values

  read(PAR,$buffer,$datasize*(1+$particles));
  @outputvalue1 = unpack("f*",$buffer);

  read(PAR,$buffer,$datasize*(1+$particles));
  @outputvalue2 = unpack("f*",$buffer);

  read(PAR,$buffer,$datasize*(1+$particles));
  @outputvalue3 = unpack("f*",$buffer);

  read(PAR,$buffer,$datasize*(1+$particles));
  @outputvalue4 = unpack("f*",$buffer);

# Printing out

# Header lines

print OUT "# Nb_of_particles=$particles Nb_of_data=$datacolumns \n" ;
print OUT "# Xloc  |  Z loc  |  Y loc   |  Visc     ";

for($i=$dim+3;$i<$dim+3+$datacolumns;$i++) {
  print OUT $dataname[$i],"      |     " ;
}
print OUT "\n" ;

for($i=1;$i<=$particles;$i++) {
#  printf OUT "$i @material[$i]  @weight[$i]  @Xcoord[$i]  @Zcoord[$i]  @outputvalue[$i] \n";
  printf OUT "%2.5f  %2.5f %2.5f %3.f %10.8e %10.8e %10.8e %10.8e \n",
  @Xcoord[$i],@Zcoord[$i],@Ycoord[$i],@material[$i],@outputvalue1[$i] ;
}



exit;


$particles=$oneline[0];
$polygons=$oneline[1];

#printf "Number of particles = %d\n",$particles;
#printf "Number of polygons =  %d\n",$polygons;

for($i=0;$i<$particles;$i++) {
  @oneline = split(' ',<PAR>);
  $X[$i] = $oneline[0];
  $Y[$i] = $oneline[1];
  $R[$i] = $oneline[2];
  $Phi[$i] = $oneline[3]; 
}

$nodeoffset[0]=0;
for($i=0;$i<=$polygons;$i++) {
  @oneline = split(' ',<PAR>);
  $polypoints[$i]=$oneline[0];
  $nodeoffset[$i+1] = $nodeoffset[$i] + $polypoints[$i];
  
  for($j=0;$j<$polypoints[$i];$j++) {
    @oneline = split(' ',<PAR>);
    $PX[$nodeoffset[$i] + $j] = $oneline[0];
    $PY[$nodeoffset[$i] + $j] = $oneline[1];
  }
}

close(PAR);


# Plotting command with scaling etc

$scale=1.0;
$psxy="psxy -P -Jx$scale -R$opt_R -A -N";
$psbasemap="psbasemap -P -Jx$scale -R$optR -B1.0";

# Output file

$outfile = sprintf("%s.%05d.ps",$opt_F,$opt_N);

# Start the plotting mechanism


# Background and/or scale bars etc

#open(PLOT,"|$psbasemap -F255/255/255 -G255/255/240 -K > $outfile");
#close(PLOT);

open(PLOT,"|$psbasemap -F255/255/255 -G255/255/255 -K > $outfile");
close(PLOT);

# Plot Particle positions

open(PLOT,"|$psxy -F0/0/0 -Sh -W0.5/150/30/0t0.3_0.2:0.2p -L $pa_colour -K -O >> $outfile");
for($i=0;$i<$particles;$i++) {
  if($opt_C) {
    printf(PLOT "%g %g %g \n",$X[$i],$Y[$i],$Phi[$i],$R[$i]*$scale*1.1); # optional shading argument
  }
  else {
    printf(PLOT "%g %g %g \n",$X[$i],$Y[$i],$R[$i]*$scale*1.1); 
  }
}
close(PLOT);


# Plot Line Boundary positions

for($i=0;$i<$polygons;$i++) {
  open(PLOT,"|$psxy -G0/255/0 -L -O -K >> $outfile");
  for($j=0;$j<$polypoints[$i];$j++) {
    printf(PLOT "%g %g \n",$PX[$nodeoffset[$i]+$j],$PY[$nodeoffset[$i]+$j]);
  }
  close(PLOT);
}

# Large roller
open(PLOT,"|$psxy -Sc -G0/200/100 -O -K -W5w0/0/255 -L >> $outfile");
printf(PLOT "%g %g %g\n",1.0,0.8,0.6);  
printf(PLOT "%g %g %g\n",4.0,0.8,0.6);  
close(PLOT);

# inner roller
open(PLOT,"|$psxy -Sc -G0/100/30 -O -K >> $outfile");
printf(PLOT "%g %g %g\n",1.0,0.8,0.4);  
printf(PLOT "%g %g %g\n",4.0,0.8,0.4);  
close(PLOT);


# OK. The postscript file is finished, do we need to 
# delete the bulky original data or compress the new file ?

if($opt_D) {
  unlink($opt_f);
}

if($opt_Z) {
  system("gzip -f $outfile");
}

exit(1);
