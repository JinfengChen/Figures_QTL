#!/usr/bin/perl
use Getopt::Long;
use File::Basename;
use Data::Dumper;
use FindBin qw($Bin);


GetOptions (\%opt,"mping:s","project:s","help");


my $help=<<USAGE;
perl $0 --mping 
Convert reloceTE confident to gff
input:
mping	TTA	HEG4_2	Chr1:588832..588834	+	5	5	GCATCTTCTTGCATTGGTAGCAAGAAAACGGCAACATATGGGCCTCCGATGGAAGCCACGTCCTGTCCAATTCCCAAACTAACCCA
output:
Chr1    HEG4_2.3.mPing.all      Non_reference   588832  588834  .       .       .       ID=mPing_1;
USAGE


if ($opt{help} or keys %opt < 1){
    print "$help\n";
    exit();
}

$opt{project} ||="RelocaTE";
readtable($opt{mping});

#mping	TTA	HEG4_2	Chr1:588832..588834	+	5	5	GCATCTTCTTGCATTGGTAGCAAGAAAACGGCAACATATGGGCCTCCGATGGAAGCCACGTCCTGTCCAATTCCCAAACTAACCCA
sub readtable
{
my ($file)=@_;
my %hash;
my $prefix=basename($file,".txt");
my $cout;
open OUT, ">$prefix.gff";
open IN, "$file" or die "$!";
while(<IN>){
    chomp $_;
    next if ($_=~/^$/ or $_=~/^TE/);
    my @unit=split("\t",$_);
    $cout++;
    my ($chr,$start,$end);
    if ($unit[3]=~/(\w+)\:(\d+)\.\.(\d+)/){
       $chr  =$1;
       $start=$2;
       $end  =$3;
    }
    print OUT "$chr\t$unit[2]\t$opt{project}\t$start\t$end\t.\t.\t.\tID=mPing_$cout;Strain=$unit[2];TSD=$unit[1];\n";
    #print OUT "$unit[3]\t$unit[2]\t$opt{project}\t$start\t$end\t.\t.\t.\tID=mPing_$cout;\n";
}
close IN;
return \%hash;
}
