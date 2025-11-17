#!/usr/bin/perl-w
use strict;

open(FH1, "$ARGV[0]") or die $!;
my $out_fh;
open $out_fh,'>',"$ARGV[1]" or die $!;

my %h;
while(<FH1>){
   chomp;
   my @a=split(/\t/);
   if($a[5] >=0.5 && $a[4] >=25 && $a[6]<=30 && $a[6]>=5){
     if(($a[14]=~ /CDS/) && ($a[14]=~ /,/)){
     	if ($a[15]=~/Missense/){
        	$h{Missense}=$h{Missense}+1;
        	next; 
        }
        if($a[15]=~/Nonsense/){
        	$h{Nonsense}=$h{Nonsense}+1;
        	next;
        }  
        if($a[15]=~/Frame shift/){
        	$h{"Frame shift"}=$h{"Frame shift"}+1;
        	next;
        }
        else{
             $h{Silent}=$h{Silent}+1;
             next;
        }
     }   
     if(($a[14]=~ /CDS/) && ($a[14]!~ /,/)){
     	  $h{$a[15]}=$h{$a[15]}+1;
          next; 
     }                  
     if($a[14]=~ /splice site/){
        $h{"splice site"}=$h{"splice site"}+1; 
        next;
     }
     else{
        next;
     }  
  }
}

foreach my $i (sort keys %h){
print $out_fh $i,"\t",$h{$i},"\n";
}

close FH1;
close out_fh;


