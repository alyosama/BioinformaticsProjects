#!/usr/bin/perl

# author Aly Osama
# id 11p6006
# email alyosamah@gmail.com
# Reading the input sequences file
# o sequence_1.txt, then sequence_2.txt, and finally sequence_3.txt.
# Initializing the scores matrix.
# Computing iteratively the scores.
# Backtracking and deducing the alignment.
# Display the score and the alignment on the screen.
################################################################################

use strict;
use warnings;
use List::Util qw[min max];
######################### Reading the input sequences file #####################
# o sequence_1.txt, then sequence_2.txt, and finally sequence_3.txt.

my $testDir="sequence_test_data/";
my @filenames=("Sequence_Pair_1.txt","Sequence_Pair_2.txt","Sequence_Pair_3.txt");

#### Read Sequence Pairs ######

my @sequence_pairs=();
for(my $i=0;$i<scalar @filenames;++$i){
  open FILE , $testDir.$filenames[$i] or die "Invalid File Name";
  my $x=<FILE>; $x =~ s/\s+//g;
  unless(length $x == ($x =~ tr/AaCcGgTt//) ){
    print "Invalid First DNA Sequence in ",$filenames[$i];
    exit;
  }
  push @{ $sequence_pairs[$i] }, uc $x;

  my $y=<FILE>; $y =~ s/\s+//g;
  unless(length $y == ($y =~ tr/AaCcGgTt//) ){
    print "Invalid Second DNA Sequence in ",$filenames[$i];
    exit;
  }
  push @{ $sequence_pairs[$i] }, uc $y;
  close FILE;
}


#### Perform Alignment #####

for(my $i=0;$i<scalar @filenames;++$i){
  print "The Alignment between $sequence_pairs[$i][0] and $sequence_pairs[$i][1] is : \n\n";
  my $maxScore=smith_waterman($sequence_pairs[$i][0],$sequence_pairs[$i][1]);
  print "Alignment Score is $maxScore\n\n";

}

################################################################################
####################### Local Alignment Fuction ###############################
################################################################################
sub smith_waterman{
  # Match Score: m=5
  # Mismatch Score: s=-4
  # Linear Gap penalty (both internal and terminal gaps): d=-5
  my $m=5;
  my $s=-4;
  my $d=-5;
  my ($A,$B)=@_;
  my @A = split( '', $A );
  my @B = split( '', $B );
  ############################# Initialize the array ###########################
  ######################## Computing iteratively the scores ####################
  my $array = [ ];

  for(my $i=0; $i <= (scalar @A); $i++) {
	  for(my $j=0; $j <= (scalar @B) ; $j++) {
      if($i==0 or $j==0){
		     $array->[$i][$j]=0;
      }
      else{
        if($A[$i-1] eq $B[$j-1]){
          $array->[$i][$j]=max(0,$array->[$i-1][$j-1]+$m,$array->[$i][$j-1]+$d,$array->[$i-1][$j]+$d);
        }
        else{
          $array->[$i][$j]=max(0,$array->[$i-1][$j-1]+$s,$array->[$i][$j-1]+$d,$array->[$i-1][$j]+$d);
        }

      }
	  }
  }

  ########################## Deducing the alignment.############################
  my @sequence_A=();
  my @sequence_B=();
  my @matching_str=();

  my $a=scalar @A;
  my $b=scalar @B;

  ############### Find Max element at last col or row ##########################
  my $maxIndexA=$a;
  my $maxIndexB=$b;
  my $maxScore=$array->[$a][$b];

  for(my $i=0; $i <= (scalar @A); $i++) {
    for(my $j=0; $j <= (scalar @B) ; $j++){
      if($maxScore<$array->[$i][$j]){
        $maxIndexA=$i;
        $maxIndexB=$j;
        $maxScore=$array->[$i][$j];
      }
    }
  }


################ Backtracking and deducing the alignment.#####################
$a=$maxIndexA;
$b=$maxIndexB;
do{
  if(($array->[$a][$b] == $array->[$a][$b-1]+$d ) and ($b!=0)){
    unshift @sequence_A , '-';
    unshift @sequence_B , $B[$b-1];
    unshift @matching_str , ' ';
    $b--;
  }
  elsif(($array->[$a][$b] == $array->[$a-1][$b]+$d) and ($a!=0)){
    unshift @sequence_A , $A[$a-1];
    unshift @sequence_B , '-';
    unshift @matching_str , ' ';
    $a--;
  }
  elsif(($array->[$a][$b] == $array->[$a-1][$b-1] + $m) and ($a!=0 or $b!=0) ){
    unshift @sequence_A , $A[$a-1];
    unshift @sequence_B , $B[$b-1];
    unshift @matching_str , '|';
    $a--;
    $b--;
  }
  elsif (($array->[$a][$b] == $array->[$a-1][$b-1]+$s)and ($a!=0 or $b!=0)){
    unshift @sequence_A , $A[$a-1];
    unshift @sequence_B , $B[$b-1];
    unshift @matching_str , '.';
    $a--;
    $b--;
  }
  else{
    $a=0;
    $b=0;
  }

} until($a==0 or $b==0);

  ################## Display the alignment on the screen #######################
  my $sequence_A= join( '', @sequence_A);
  my $matching_str= join( '', @matching_str);
  my $sequence_B= join( '', @sequence_B);
  print "$sequence_A\n$matching_str\n$sequence_B\n";


  ################################ Print 2D Array ##############################

  # print "@$_\r\n" for @$array;

  # printf "%3s %3s ", " "," ";
  # printf "%3s ",$_ for @B;
  # print "\r";
  # for(my $i=0; $i <= (scalar @A); $i++) {
  #   if($i!=0){printf "%3s ", $A[$i-1]; }
  #   else{printf "%3s ", " ";}
  #
	#   for(my $j=0; $j <= (scalar @B); $j++) {
  #       printf "%3d ",$array->[$i][$j];
  #   }
  #   print "\n";
  # }

  #printf "%3d %3d \n", scalar @A,scalar @B;

  return $maxScore;

};
