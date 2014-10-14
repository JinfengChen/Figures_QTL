#!/usr/bin/perl

use strict;
use FindBin qw ($Bin);
use Getopt::Long;
use SVG;

#################################Usage ##############################################
my %opt;
GetOptions(\%opt,"table:s","help:s");

my $help=<<USAGE;
perl $0 --table ../input/QTL.regions.list

USAGE

if(defined $opt{help} or keys %opt < 1){
        die  $help ;
}



my $svg=SVG->new(width=>1500,height=>800);
my $startw=100; 
my $endw  =1400;
my $starth=100;
my $endh  =700;

my %chrs=(
   'Chr1' => 1,
   'Chr3' => 1,
   'Chr5' => 1,
   'Chr6' => 1,
   'Chr8' => 1,
   'Chr10'=> 1
);

my $refchr=readtable("../input/MSU7.chrlen", \%chrs);
my $refcent=readtable("../input/MSU7.cent", \ %chrs);
my $scale=getscale($refchr,"600");
my @qtl=readQTLtable($opt{table});

#### chromosome line
my %refpos;
my %qrypos;
my @chr=sort {$a <=> $b} keys %$refchr;
for (my $i=0;$i<@chr;$i++){
    ###################reference chromosome
    my $refx=$i*220+$startw;
    my $refy=$starth;
    my $refw=20;
    my $refh=$refchr->{$chr[$i]}/$scale;
    $refpos{$chr[$i]}=$refx;
    #print "$i\t$chr[$i]\t$refx\n";
    my $refline=$svg->rectangle(
              x=>$refx,y=>$refy,
              width=>$refw,height=>$refh,
              rx=>8,ry=>8,
              style=>{
                   fill=>'indianred',
                   stroke=>'none'
              }
              
    );
    my $refcx=$refx+$refw/2;
    my $refcy=$refcent->{$chr[$i]}/$scale+$refy;
    my $refrx=$refw/2;
    my $refry=$refrx*1.2;
    $svg=centromere($svg,$refcx,$refcy,$refrx,$refry,"lightblue");
    my $refnote  =$svg->text(
              x=>$refx+15, y=>$starth-20,
              style=>{
                   fontsize=>'2','text-anchor'=>'middle','stroke-width'=>1
              }
    )->cdata('Chr'.$chr[$i]);
}


#####draw QTL block
my @position;
foreach my $refqtl (@qtl){
    my $chrx= $refpos{$refqtl->[0]};
    my $chry= $starth;
    my $rank = overlap(\@position, $refqtl->[0], $refqtl->[1], $refqtl->[2]);
    my @temp = ($refqtl->[0], $refqtl->[1], $refqtl->[2]);
    push (@position, \@temp);
    #print "$refqtl->[0], $refqtl->[1], $refqtl->[2], $refqtl->[3], $rank\n"; 
    if ($refqtl->[3] =~ /Heading/){
       my $x = $chrx + 30 + $rank*20;
       my $y = $chry + $refqtl->[1]/$scale;
       my $w = 15;
       my $h = ($refqtl->[2] - $refqtl->[1])/$scale;
       my $color = 'orchid';
       drawrect($x, $y, $w, $h, $color);
    }else{
       my $x = $chrx + 30 + $rank*20;
       my $y = $chry + $refqtl->[1]/$scale;
       my $w = 15;
       my $h = ($refqtl->[2] - $refqtl->[1])/$scale;
       my $color = 'darkslateblue';
       drawrect($x, $y, $w, $h, $color);

    }
}



writesvg("QTL_DRAW.svg",$svg);

#####

sub overlap{
my ($pos, $chr, $start, $end) = @_;
my $o = 0;
if (@$pos == 0){
   return 0;
}else{
   for(my $i=0; $i< @$pos; $i++){
       my $tempc = $pos->[$i]->[0];
       my $temps = $pos->[$i]->[1];
       my $tempe = $pos->[$i]->[2];
       if ($tempc eq $chr){
          if ($start >= $temps and $start <= $tempe){
              $o++;   
          }elsif($end >= $temps and $end <= $tempe){
              $o++;
          }elsif($temps >= $start and $temps <= $end){
              $o++;
          }elsif($tempe >= $start and $tempe <= $end){
              $o++;
          }
       }
       #print "$tempc\t$temps\t$tempe\t$o\n";
   }
}
return $o;
}

sub drawrect{
my ($x, $y, $w, $h, $color)= @_;
my $qtl=$svg->rectangle(
   x=>$x, y=>$y,
   width=>$w, height=>$h,
   style=>{
       fill=>$color,
       stroke=>'none'
   }
);   


}



sub getscale
{
my ($chr,$height)=@_;
my @len=sort {$b <=> $a} values %$chr;
#print "$len[0]\n";
my $rate=$len[0]/$height;
return $rate;
}


sub centromere
{
my ($svg,$cx,$cy,$rx,$ry,$fillcolor)=@_;
my $tag = $svg->ellipse(
        cx=>$cx, cy=>$cy,
        rx=>$rx, ry=>$ry,
        style=>{
            'fill'=>$fillcolor,
        }
    );
return $svg;
}


sub readtable
{
my ($file, $chrs)=@_;
my %hash;
open IN, "$file" or die "$!";
while(<IN>){
    chomp $_;
    next if ($_=~/^$/);
    my @unit=split("\t",$_);
    next unless exists $chrs->{$unit[0]};
    $unit[0] = $1 if ($unit[0]=~/(\d+)/);
    $hash{$unit[0]}=$unit[1];
}
close IN;
return \%hash;
}

#Chr	Start	End	Name	Length	Trait
#Chr3	962389	1665631	HeadingDate1	703243	HeadingDays
#Chr6	1989540	2342014	HeadingDate2	352475	HeadingDays
sub readQTLtable
{
my ($file)=@_;
my @qtl;
open IN, "$file" or die "$!";
while(<IN>){
    chomp $_;
    next if ($_=~/^$/);
    my @unit=split("\t",$_);
    $unit[0] = $1 if ($unit[0]=~/(\d+)/);
    my @temp = ($unit[0], $unit[1], $unit[2], $unit[5]);
    push (@qtl, \@temp);
}
close IN;
return @qtl;
}




################################### sub for write svg to file
sub writesvg {
my ($file,$svg)=@_;
#print "$file\n";
open OUT, ">$file" or die "can not open my file";
       print OUT $svg->xmlify;
close OUT;
system "/rhome/cjinfeng/software/tools/draw/svg2xxx_release/svg2xxx $file -t pdf";
}


