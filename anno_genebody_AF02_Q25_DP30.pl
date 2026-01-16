#!/usr/bin/perl-w
use strict;

open(FH1, "$ARGV[0]") or die $!;
my $out_fh;
open $out_fh,'>',"$ARGV[1]" or die $!;

my $cds;
my $intron;
my $utr;
my $intergenic;
my $exon;
my $c;
while(<FH1>){
   chomp;
   my @a=split(/\t/);
   if($a[5] >=0.2 && $a[4] >=25 && $a[6] <=30 && $a[6] >=5){
     if(($a[14]=~ /CDS/) or ($a[14]=~ /splice site/) or ($a[14]=~ /exon/)){
        $cds=$cds+1;
        next; 
       } 
     elsif($a[14]=~ /intron/){
        $intron=$intron+1;
        next; 
       } 
     else{
        $intergenic=$intergenic+1;
        next; 
       } 
     }
     else{
        next;
     }  
}

my $sum = $cds+$intron+$intergenic;
my $x = $cds/$sum;
my $y = $intron/$sum;
my $z = $intergenic/$sum;


print $out_fh "sum", "\t", $sum, "\n", "CDS", "\t", $cds, "\t", $x, "\n", "Intron", "\t", $intron, "\t", $y, "\n", "Intergenic", "\t", $intergenic, "\t", $z, "\n";

close FH1;
close out_fh;


