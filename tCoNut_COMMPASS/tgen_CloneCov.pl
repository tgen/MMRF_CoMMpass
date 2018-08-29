#!/usr/bin/perl -w
###########################################################################
##
##  Program:  tgen_CloneCov.pl
##  Author:   David W. Craig, Ph.D.
##  Created:  2012
##
##
##  $Id: tgen_CloneCov.pl,v 0.09 dcraig Exp $
##
############################################################################
$VERS=0.092;

#######################
# Variables & Main Loop
#######################
$|=1;
$bin=100;
$mqf=20;
$flags="-f 0x0002 -F 0x0400 -F0x0200 -q $mqf"; 
$good=0;$dup=0;$poor_map=0;$skip=0;
&command_arg;
$i_min=10; $i_max=2000;
($i_min,$i_max)=calc_insert($in,$i_min,$i_max);
print "Version: $VERS $skip $mqf\n";
&init_clone_array;
&main_loop;
&print_out;
#################
# Main Loop
#################
sub main_loop {
  $c=0;
  open (IN,"$st view $flags $in $lo |") || die "Can't open $in\n";
  print "-Start Processing\n";
  LOOP: while (<IN>) {
    $line=$_;
    @temp=split(/\t/,$line);
    ++$c;
    if ($c % 1000000 == 0) {print "Read: $c at $temp[2]\t$temp[3]\n";}
    $temp[2]=~s/chr//g;
    if ($temp[2] eq "X") { $temp[2] = 23;}
    if ($temp[2] eq "Y") { $temp[2] = 24;}
    if ($temp[2] eq "MT") { $temp[2] = 25;}
    if ($temp[2] =~/\D/) {last LOOP; next LOOP;}
    if ($temp[8]>$i_min && $temp[8]<$i_max) { 
       ++$good;
       $st=int($temp[3]/$bin);
       $en=int($temp[7]/$bin);
       $chr=$temp[2];
       for ($j=$st;$j<=$en;++$j) {
         ++$cln[$chr][$j];
       }
       next LOOP;
    }
    if ($temp[4] < $mqf) {++$poor_map;next LOOP;} # Poor Mapping
  }
  close (IN);
  print "Good_reads:\t$good\n";
  print "Duplicates:\t$dup\n";
  print "Poor_MapQ:\t$poor_map\n";
  print "Total_Reads:\t$c\n";
}
#################
# Print Out
#################
sub print_out{
  print "-Printing output\n";
  open (OUT,">$out.dat") || die "Can't open $out.dat\n";
  for ($i=1;$i<=25;++$i) { 
    for ($j=0;$j<=$ln_chr[$i];++$j) {
      $pos=$j*$bin;
      if ($cln[$i][$j] >0 || $nz==0) {
        if (exists($cln[$i][$j] )) {
          print OUT "$i\t$pos\t$cln[$i][$j]\n";
        }
      }
    }
  }
  close (OUT);
}


#################
# Initialize
#################

sub init_clone_array {
  print "-Initializing\n";
  open (IN,"$st view -H $in |") || die "Can't open $in\n";
  while (<IN>) {
    chomp;
    @buf=split(/\t/);
    if ($buf[1] =~/SN:/) {
      $chr=$buf[1];
      $chr=~s/SN://g;
      $chr=~s/chr//g;
      $len=$buf[2];
      $len=~s/LN://g;
      $len=int($len/$bin)+1;
      if ($chr eq "X") {$chr=23}
      if ($chr eq "Y") {$chr=24}
      if ($chr eq "MT") {$chr=25}
      if (!($chr =~/\D/) ) { 
        $ln_chr[$chr]=$len;
        print "Initial chr: $chr to $len\n";
        for ($i=0;$i<=$len;++$i) { 
          $cln[$chr][$i]=0;
        }
      } 
    }
  }
  close (IN);
  print "-Done Initializing\n";
}
sub command_arg {
  $td="/Home/dcraig";
  $out="out.cln";
  $lo="";
  $nz=0;
  $ma="PG:Z";
  if ($#ARGV < 1) {
    print "./tgen_clone_coverage.v006.pl I=aut.twin.4982.bam T=/Home/dcraig O=aut.twin.4982 > aut.twin.4982.c.out\n"; 
    die;
  } else {
    foreach $a (@ARGV) {
      if ($a=~/I=/) {
        $in=$a;
        $in=~s/I=//g;
        print "-Input_file:\t$in\n";
       } elsif ($a=~/T=/) {
        $td=$a;
        $td=~s/T=//g;
        print "-Temp_dir:\t$in\n";
       } elsif ($a=~/L=/) {
        $lo=$a;
        $lo=~s/L=//g;
        print "-Location:\t$lo\n";
       } elsif ($a=~/O=/) {
        $out="$a.cln";
        $out=~s/O=//g;
        print "-Out_head:\t$in\n";
       } elsif ($a=~/M=/) {
        $ma=$a;
        $ma=~s/M=//g;
        print "-Match: $ma\n";
       } elsif ($a=~/S=/) {
        $st=$a;
        $st=~s/S=//g;
        print "-samtools: $st\n";
       } elsif ($a eq "--nz") {
        $nz=1;
        print "-Only outputing non-zero coverage\n";
       }else {die "DIE: Can't parse $a\n";}
    }
  } 
}
sub calc_insert {
 print "-Calculating insert size for $_[0]\n";
 my $in=$_[0];
 my $i_min=$_[1];
 my $i_max=$_[2];
 open (IN,"$st view $flags -s 0.001 $in |") || die "Can't open $in\n";
 my @insert_vals=();
 print "-Before $in max_insert=$i_max\tmin_insert=$i_min\n";
 LOOP:while (<IN>) {
    my @temp=split(/\t/);
    if ($temp[6] eq "*") {next LOOP;}
    if (abs($temp[8]) > $i_min && abs($temp[8]) < $i_max) {
      push (@insert_vals,abs($temp[8]));
    } #else { print "temp8 $temp[8] - $temp[6]\n"}
    if ($#insert_vals >10000) {last LOOP}
 }
 close (IN);
 @sort_insert_vals=sort { $a <=> $b } @insert_vals;
 $i_min=$sort_insert_vals[10];
 $i_max=$sort_insert_vals[9990];
 $i_med=$sort_insert_vals[5000];
print "number of vals: $#sort_insert_vals $sort_insert_vals[9990] $sort_insert_vals[9900]\n";
 if ($i_min<10 || $i_min>500) {$i_min=10}
 if ($i_max>2000 || $i_max<500) {$i_max=2000}
 print "-After: $in max_insert=$i_max\tmin_insert=$i_min\tmed_insert=$i_med\n";
 return ($i_min,$i_max);
}


