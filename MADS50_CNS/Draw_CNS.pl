#!/usr/bin/perl
use Getopt::Long;
use SVG;
my %opt;

GetOptions (\%opt,"headers:s","project:s","help");


my $help=<<USAGE;
perl $0 --project 
USAGE


if ($opt{help} or keys %opt < 1){
    print "$help\n";
    exit();
} 

$opt{project} ||="Intron_CNS";


our $height=600;
our $width=800;
my $svg=SVG->new(width=>$width,height=>$height);

my $length = 30074;
my $igs = 6000;
my $ige = 25000;
my $xstart = 20;
my $xend   = $width-20;
my $ystart = 50;
my $yend   = $height-50;
my $yint1  = 40;
my $yint2  = 50;
my $rate = ($width-2*$xstart)/($length-($ige-$igs+1));

my @genomes= ('Brachyantha', 'Brachypodium', 'Setaria', 'Sorghum');
#ReadAlign('MADS50_align.fa', \@genomes); # delete line 552 in sorghum identity
$svg = Draw_CNS($svg, 'Brachyantha', $rate, 'MADS50_Brachyantha.identity', 5*$yint2+$ystart, $xstart, $xend, $igs, $ige);
$svg = Draw_CNS($svg, 'Brachypodium', $rate, 'MADS50_Brachypodium.identity', 6*$yint2+$ystart, $xstart, $xend, $igs, $ige);
$svg = Draw_CNS($svg, 'Setaria', $rate, 'MADS50_Setaria.identity', 7*$yint2+$ystart, $xstart, $xend, $igs, $ige);
$svg = Draw_CNS($svg, 'Sorghum', $rate, 'MADS50_Sorghum.identity', 8*$yint2+$ystart, $xstart, $xend, $igs, $ige);

$svg = Draw_gene($svg, 'Nipponbare', $rate, 'MADS50_NB.anno', 0*$yint1+$ystart, $xstart, $xend, $igs, $ige);
$svg = Draw_gene($svg, 'HEG4', $rate, 'MADS50_HEG4.anno', 1*$yint1+$ystart, $xstart, $xend, $igs, $ige);
$svg = Draw_gene($svg, 'EG4', $rate, 'MADS50_EG4.anno', 2*$yint1+$ystart, $xstart, $xend, $igs, $ige);
$svg = Draw_gene($svg, 'A119', $rate, 'MADS50_A119.anno', 3*$yint1+$ystart, $xstart, $xend, $igs, $ige);
$svg = Draw_gene($svg, 'A123', $rate, 'MADS50_A123.anno', 4*$yint1+$ystart, $xstart, $xend, $igs, $ige);

#my @genomes= ('Brachyantha', 'Brachypodium', 'Setaria', 'Sorghum');
#ReadAlign('MADS50_align.fa', \@genomes); # delete line 552 in sorghum identity
#$svg = Draw_CNS($svg, 'Brachyantha', $rate, 'MADS50_Brachyantha.identity', 5*$yint2+$ystart, $xstart, $xend, $igs, $ige);
#$svg = Draw_CNS($svg, 'Brachypodium', $rate, 'MADS50_Brachypodium.identity', 6*$yint2+$ystart, $xstart, $xend, $igs, $ige);
#$svg = Draw_CNS($svg, 'Setaria', $rate, 'MADS50_Setaria.identity', 7*$yint2+$ystart, $xstart, $xend, $igs, $ige);
#$svg = Draw_CNS($svg, 'Sorghum', $rate, 'MADS50_Sorghum.identity', 8*$yint2+$ystart, $xstart, $xend, $igs, $ige);

$svg= mPing_pos($svg, 735, 24, '-122');
$svg= mPing_pos($svg, 690, 24, '-784');
$svg= mPing_pos($svg, 665, 24, '-787');


my $outfile="$opt{project}.svg";
writesvg($outfile,$svg);

#############

sub mPing_pos{
my ($svg, $x, $y, $p) = @_;
$svg->text( 
    x=> $x, y=> $y,
    style=>{
           'font-size'=>'70%','text-anchor'=>'start','font-weight'=>'100'
    }
)->cdata($p);
return $svg;
}


sub Draw_gene{
my ($svg, $title, $rate, $file, $y, $xstart, $xend, $igs, $ige) = @_;
my @anno = readanno($file);
$svg = draw_break_line($svg, $rate, $y, $xstart, $xend, $igs, 0.2);
$svg->text(
            x=>$xstart+20,y=>$y-20,
            style=>{
               'font-size'=>'70%','text-anchor'=>'start','font-weight'=>'100'
            }
)->cdata($title);
foreach my $feature (@anno){
    next if ($feature->[0] > $igs and $feature->[0] < $ige);
    next if ($feature->[1] > $igs and $feature->[1] < $ige);
    my $s = $feature->[0] < $igs ? $xstart+$feature->[0]*$rate : $xstart + ($feature->[0]-($ige-$igs+1))*$rate;
    my $l = ($feature->[1]-$feature->[0]+1)*$rate;
    my $c = $feature->[2] eq 'exon' ? 'blue' : 'aquamarine';
    my $h = 8;
    if ($feature->[2] =~ /mPing/){
       my $pos = $svg->line(
            x1=>$s, y1=> 30,
            x2=>$s, y2=> 450,
            style=>{
                stroke=>'gray',
                'stroke-dasharray' => "2,2", 
                d => "M5 40 l215 0"
            }
       );
       $s = $s - 430*$rate + 1;
       $l = 430*$rate;
       $c = 'orange';
       if ($feature->[2] eq 'mPing'){
           my $mping = $svg->rectangle(
                  x=>$s, y=>$y-$h/2, 
                  width=>$l, height=>$h,
                  style=>{
                     fill=>$c
                  }
           );
       }else{
           my $mping0 = $svg->pattern(
                  id =>'transform',
                  x=>$s, y=>$y-$h/2,
                  width=>$l, height=>$h,
                  patternUnits => "userSpaceOnUse",
                  patternTransform => "rotate(30)"
           );
           
       }
    }else{
       my $exon = $svg->rectangle(
              x=>$s, y=>$y-$h/2,
              width=>$l, height=>$h,
              style=>{
                fill=>$c
              }
       );
    }
}
return $svg;
}

sub readanno{
my ($file) = @_;
my @array;
open IN,"$file" or die "$!";
while (<IN>){
    next if (length $_ < 2);
    my @unit=split(" ",$_);
    if ($unit[0] ne '<'){
        push (@array, \@unit);
    }
}
close IN;
return @array
}

sub draw_break_line{
my ($svg, $rate, $y, $xstart, $xend, $igs, $linewidth) = @_;

my $breakx1 = $igs*$rate;
my $breakx2 = $breakx1 + 5;

my $geneline1=$svg->line(
          x1=>$xstart, y1=>$y,
          x2=>$breakx1, y2=>$y,
          style=>{stroke=>'black', 'stroke-width'=> $linewidth}
       );

my $geneline2=$svg->line(
          x1=>$breakx2, y1=>$y,
          x2=>$xend+5, y2=>$y,
          style=>{stroke=>'black', 'stroke-width'=> $linewidth}
       );

my $break1 = $svg->line(
          x1=>$breakx1+2, y1=>$y-2,
          x2=>$breakx1-2, y2=>$y+2,
          style=>{stroke=>'black'}
       );
my $break2 = $svg->line(
          x1=>$breakx2+2, y1=>$y-2,
          x2=>$breakx2-2, y2=>$y+2,
          style=>{stroke=>'black'}
       );
return $svg;
}

sub ReadAlign{
my ($file, $genomes) = @_;
$/=">";
my %hash;
open IN,"$file" or die "$!";
while (<IN>){
    next if (length $_ < 2);
    my @unit=split("\n",$_);
    my $temp=shift @unit;
    my @temp1=split(" ",$temp);
    my $head=$temp1[0];
    my $seq=join("\n",@unit);
    $seq=~s/\>//g;
    $seq=~s/\n//g;
    $seq=~s/\s//g;
    $hash{$head} = $seq;    
}
$/="\n";

my $win = 50;
my $ref;
my $count = 0;
my $last = 0;
my @range;
my $totalbase = 0;
for(my $i=0; $i< length($hash{'Rice'}); $i++){
    my $base = substr($hash{'Rice'},$i,1);
    if ($base ne '-'){
        $count++;
        $totalbase++;
        #print "N: $count\n";
        if ($count == $win){
           my @pos = ($last, $i);
           print "$totalbase\t$last\t$i\n";
           push (@range, \@pos);
           $count = 0;
           $last = $i+1;
        }
    }
}
print "Total base of rice: $totalbase\n";

for(my $i=0; $i < @$genomes; $i++){
    open OUT, ">MADS50_$genomes->[$i].identity" or die "$!";
    for(my $j =0; $j < @range; $j++){
        my $l = $range[$j]->[1] - $range[$j]->[0] + 1;
        my $ref= substr($hash{'Rice'}, $range[$j]->[0], $l);
        my $qry= substr($hash{$genomes->[$i]}, $range[$j]->[0], $l);
        my $id = identity($ref, $qry);
        my $rs = $j*$win;
        my $re = $rs + $win;
        print OUT "$rs\t$re\t$id\n";
    }
}

}


sub identity{
my ($ref, $qry) = @_;
my $count=0;
for(my $i =0; $i< length($ref); $i++){
    my $refb = substr($ref, $i, 1);
    my $qryb = substr($qry, $i, 1);
    if ($refb ne '-' and $refb eq $qryb){
       $count++;
    }
}
#my $id = $count/length($ref);
my $id = $count/50;
return $id;
}


sub Draw_CNS{
my ($svg, $title, $rate, $file, $y, $xstart, $xend, $igs, $ige) = @_;
my @cns = readcns($file);
$svg = draw_break_line($svg, $rate, $y, $xstart, $xend, $igs, 1);
$svg->text(
            x=>$xstart+20,y=>$y-40,
            style=>{
               'font-size'=>'70%','text-anchor'=>'start','font-weight'=>'100'
            }
)->cdata($title);
my $plotheight = 30;
foreach my $feature (@cns){
    next if ($feature->[0] > $igs and $feature->[0] < $ige);
    next if ($feature->[1] > $igs and $feature->[1] < $ige);
    my $s = $feature->[0] < $igs ? $xstart+$feature->[0]*$rate : $xstart + ($feature->[0]-($ige-$igs+1))*$rate;
    my $l = ($feature->[1]-$feature->[0]+1)*$rate;
    my $c = $feature->[2] > 0.7 ? 'purple' : 'powderblue';
    my $h = $feature->[2] * $plotheight;
    my $exon = $svg->rectangle(
              x=>$s, y=>$y-5-$h,
              width=>$l, height=>$h,
              style=>{
                fill=>$c
              }
    );
}
return $svg;
}





sub readcns{
my ($file) = @_;
my @array;
open IN,"$file" or die "$!";
while (<IN>){
    next if (length $_ < 2);
    my @unit=split("\t",$_);
    if ($unit[0] ne '<'){
        push (@array, \@unit);
    }
}
return @array;
}

sub legend{
my ($svg,$x,$y,$name,$color)=@_;
my $xtext=$x+18; 
 $svg->rectangle(
            x=>$x,y=>$y,
            width=>10,height=>10,
            style=>{
                fill=>$color
            }
 );
 $svg->text(
            x=>$xtext,y=>$y+6,
            style=>{
               'font-size'=>'70%','text-anchor'=>'start','font-weight'=>'100'
            }
 )->cdata($name);
 
return $svg;
}



sub maxlen
{
my (@header)=@_;
my $max=0;
for(my $i=0;$i<@header;$i++){
   my $fa=$header[$i].".fasta";
   my $len=getfastalen($fa);
   $max= $len > $max ? $len : $max;
}
return $max;
}


#############
sub parseACT
{
my ($file)=@_;
my @com;
open IN, "$file" or die "$!";
while(<IN>){
   chomp $_;
   my @unit=split(" ",$_);
   push (@com,[$unit[2],$unit[3],$unit[5],$unit[6]]);
}
close IN;
return \@com;
}

#############
sub parsecoord
{
my ($file)=@_;
my @com;
open IN, "$file" or die "$!";
while(<IN>){
   chomp $_;
   next unless ($_=~/^\s*\d+/);
   my @unit=split(" ",$_);
   push (@com,[$unit[1],$unit[2],$unit[4],$unit[5]]);
}
close IN;
return \@com;
}


#############
sub parseblastm8
{
my ($file)=@_;
my @com;
open IN, "$file" or die "$!";
while(<IN>){
   chomp $_;
   my @unit=split("\t",$_);
   push (@com,[$unit[6],$unit[7],$unit[8],$unit[9]]);
}
close IN;
return \@com;
}

sub drawlinkACT
{
my ($act,$rate,$h1,$h2)=@_;
my $identity=50;
my $lencut=50;
my @links;
open IN, "$act" or die "$!";
while (<IN>){
    chomp $_;
    next if ($_ eq "");
    my @array=split(" ",$_);
    if ($array[1] >= $identity and $array[3]-$array[2] > $lencut ){
        push (@links,"$array[2]\t$array[3]\t$array[5]\t$array[6]");
    }
}
close IN;

foreach (@links){
    my @unit=split("\t",$_);
    my $qleft=$unit[0]/$rate+200;
    my $qright=$unit[1]/$rate+200;
    my $tleft=$unit[2]/$rate+200;
    my $tright=$unit[3]/$rate+200;
    my $color;
    if ($tright < $tleft){
         $color='red';
    }else{
         $color='#778899';
    }

    my $qheight=$h1;
    my $theight=$h2;
    my $xv=[$qleft,$qright,$tright,$tleft];
    my $yv=[$qheight,$qheight,$theight,$theight];
    my $points =$svg->get_path(
                     x=>$xv,y=>$yv,
                     -type=>'polyline',
                     -closed=>'true'
              );
    my $tag=$svg->polyline(
                     %$points,
                     style=>{
                        fill=>$color
                        #fill=>'#FFDAB9' 
                     }
              );
}
}

#######drawlink blastm8
sub drawlinkm8
{
my ($m8,$rate)=@_;
my $identity=50;
my $lencut=50;
my @links;
open IN, "$m8" or die "$!";
while (<IN>){
    chomp $_;
    next if ($_ eq "");
    my @array=split(" ",$_);
    if ($array[2] >= $identity and $array[7]-$array[6] > $lencut ){
        push (@links,"$array[6]\t$array[7]\t$array[8]\t$array[9]");
    }
}
close IN;

foreach (@links){
    my @unit=split("\t",$_);
    my $qleft=$unit[0]/$rate+50;
    my $qright=$unit[1]/$rate+50;
    my $tleft=$unit[2]/$rate+50;
    my $tright=$unit[3]/$rate+50;
    my $color;
    if ($tright < $tleft){
         $color='red';
    }else{
         $color='#778899';
    }

    my $qheight=110;
    my $theight=290;
    my $xv=[$qleft,$qright,$tright,$tleft];
    my $yv=[$qheight,$qheight,$theight,$theight];
    my $points =$svg->get_path(
                     x=>$xv,y=>$yv,
                     -type=>'polyline',
                     -closed=>'true'
              );
    my $tag=$svg->polyline(
                     %$points,
                     style=>{
                        fill=>$color
                        #fill=>'#FFDAB9'
                     }
              );
}
}



sub drawdotplot
{
my ($svg,$compare,$rate,$x,$y)=@_;
foreach my $match (@$compare){
   my $x1=$match->[0]/$rate+$x;
   my $y1=$y-$match->[2]/$rate;
   my $x2=$match->[1]/$rate+$x;
   my $y2=$y-$match->[3]/$rate;
   #print "$x1\t$x2\t$y1\t$y2\n";
   my $line=$svg->line(
     x1=>$x1,y1=>$y1,
     x2=>$x2,y2=>$y2,
     style=>{stroke=>'black'}
   );   
}
return $svg;
}


sub drawXgene
{
my ($svg,$refgenegff,$rate,$x1,$x2,$y,$head)=@_;
print "Start:$x1\tEnd:$x2\n";
my $strandline=$svg->line(
     x1=>$x1,y1=>$y,
     x2=>$x2,y2=>$y,
     style=>{stroke=>'black'}
);
my $title=$svg->text(
     x=>$x1-190,y=>$y
)->cdata("$head");

foreach my $g (keys %$refgenegff){
    my @line=split("\n",$refgenegff->{$g});
    my @pos;
    my $strand;
    foreach my $e (@line){
        my @unit=split("\t",$e);
        if ($unit[2] eq "mRNA"){
           $strand=$unit[6];
        }else{
           push (@pos,[$unit[3],$unit[4]]);
        }
    }
    @pos=sort {$a->[0] <=> $b->[1]} @pos;
    my $gstart=$pos[0][0]/$rate+$x1;
    my $gend  =$pos[$#pos][1]/$rate+$x1; 
    if ($strand eq "+"){
=pod
       my $geneid=$svg->text(
          x=>$gstart,y=>$y-50,
          style=>{
             'font-size'=>'50%','text-anchor'=>'start','font-weight'=>'100'
          }
       )->cdata("$g");
=cut
       my $geneline=$svg->line(
          x1=>$gstart,y1=>$y-6,
          x2=>$gend,y2=>$y-6,
          style=>{stroke=>'black'}
       );
       foreach my $e (sort {$a->[0] <=> $b->[1]} @pos){
           my $start=$e->[0]/$rate+$x1;
           my $elen =($e->[1]-$e->[0]+1)/$rate;
           my $exony=$y-10;
           my $exon=$svg->rectangle(
              x=>$start, y=>$exony,
              width=>$elen,height=>8,
              style=>{
                fill=>'black'
              }
           );
       }   
    }else{
=pod
       my $geneid=$svg->text(
          x=>$gstart,y=>$y-10,
          style=>{
             'font-size'=>'50%','text-anchor'=>'start','font-weight'=>'100'
          }
       )->cdata("$g");
=cut
       my $geneline=$svg->line(
          x1=>$gstart,y1=>$y+6,
          x2=>$gend,y2=>$y+6,
          style=>{stroke=>'black'}
       );
       foreach my $e (sort {$a->[0] <=> $b->[1]} @pos){
           my $start=$e->[0]/$rate+$x1;
           my $elen =($e->[1]-$e->[0]+1)/$rate;
           my $exony=$y+2;
           my $exon=$svg->rectangle(
              x=>$start, y=>$exony,
              width=>$elen,height=>8,
              style=>{
                fill=>'black'
              }
           );
       }
    }
}
return $svg;
}
#####
sub drawXTE
{
my ($svg,$reftegff,$rate,$x1,$y)=@_;
foreach my $te (keys %$reftegff){
    my @line=split("\t",$reftegff->{$te});
    my $gstart=$line[3]/$rate+$x1;
    my $gend  =$line[4]/$rate+$x1;
    my $strand=$line[6];
    my $type=$1 if ($line[8]=~/Class=(.*?);/);
        my $color="gray";
    if ($type=~/DNA/){
       $color="blue";
    }elsif($type=~/LTR/){
       $color="red";
    }
    $type=~s/DNA\///;
    $type=~s/LTR\///;
    #print "$te\t$y\t$gstart\t$gend\t$strand\t$type\n";
    if ($strand eq "+"){
=pod
       my $geneid=$svg->text(
          x=>$gstart,y=>$y-14,
          style=>{
             'font-size'=>'50%','text-anchor'=>'start','font-weight'=>'100'
          }
       )->cdata("$type");
=cut 
             my $qleft =$gstart;
              my $qright=$gend;
              my $tright=$gend;
              my $tleft =$gstart;
              my $qheight=$y-2;
              my $theight=$y-10;
              my $xv=[$qleft,$qright,$tright,$tleft];
              my $yv=[$qheight,$qheight,$theight,$theight];
              my $points =$svg->get_path(
                     x=>$xv,y=>$yv,
                     -type=>'polyline',
                     -closed=>'true'
              );
              my $tag=$svg->polyline(
                     %$points,
                     style=>{
                        fill=>$color
                     }
              );

    }else{
=pod
       my $geneid=$svg->text(
          x=>$gstart,y=>$y+14,
          style=>{
             'font-size'=>'50%','text-anchor'=>'start','font-weight'=>'100'
          }
       )->cdata("$type");
=cut
              my $qleft =$gstart;
              my $qright=$gend;
              my $tright=$gend;
              my $tleft =$gstart;
              my $qheight=$y+2;
              my $theight=$y+10;
              my $xv=[$qleft,$qright,$tright,$tleft];
              my $yv=[$qheight,$qheight,$theight,$theight];
              my $points =$svg->get_path(
                     x=>$xv,y=>$yv,
                     -type=>'polyline',
                     -closed=>'true'
              );
              my $tag=$svg->polyline(
                     %$points,
                     style=>{
                        fill=>$color
                     }
              );      
    }
}
return $svg;
}



####
sub drawbox
{
my ($svg,$x1,$x2,$y1,$y2)=@_;
my $hline=$svg->line(
     x1=>$x1,y1=>$y2,
     x2=>$x2,y2=>$y2,
     style=>{stroke=>'black'}
);
my $vline=$svg->line(
     x1=>$x1,y1=>$y1,
     x2=>$x1,y2=>$y2,
     style=>{stroke=>'black'}
);
return $svg;
}


#########
sub drawxaxis
{
my ($svg,$x1,$x2,$y,$min,$max,$step,$rate)=@_;
#print "$x1\t$x2\t$y\n";
my $xaxis=$svg->line(
     x1=>$x1,y1=>$y,
     x2=>$x2,y2=>$y,
     style=>{stroke=>'black'}
);

my $bottomline=$svg->line(
     x1=>$x2,y1=>$y,
     x2=>$x2,y2=>$y+5,
     style=>{stroke=>'black'}
);
my $tail =int ($max/1000);
=pod
my $bottomtext=$svg->text(
     x=>$x2,y=>$y+20,
     style=>{
         'font-size'=>'70%','text-anchor'=>'start','font-weight'=>'100'
     }
)->cdata("$tail kb");
=cut

for(my $i=$min;$i<$max;$i=$i+$step){
     my $tempx=$x1+($i-$min+1)/$rate;
     #print "$tempx\t$min\t$step\t$i\t$rate\n";
     #print "$tempx\t$y\n";
     my $line=$svg->line(
         x1=>$tempx,y1=>$y,
         x2=>$tempx,y2=>$y+5,
         style=>{stroke=>'black'}
     );
     my $tempi=int ($i/1000);
     my $text=$svg->text(
         x=>$tempx+3,y=>$y+20,
         style=>{
             'font-size'=>'70%','text-anchor'=>'end','font-weight'=>'100'
         }
     )->cdata($tempi);
}
return $svg;
}

#########
sub drawyaxis
{
my ($svg,$y1,$y2,$x,$min,$max,$step,$rate)=@_;
#print "$x1\t$x2\t$y\n";
my $yaxis=$svg->line(
     x1=>$x,y1=>$y1,
     x2=>$x,y2=>$y2,
     style=>{stroke=>'black'}
);

my $bottomline=$svg->line(
     x1=>$x,y1=>$y2,
     x2=>$x+5,y2=>$y2,
     style=>{stroke=>'black'}
);
my $tail =int ($max/1000);
=pod
my $bottomtext=$svg->text(
     x=>$x,y=>$y2+20,
     style=>{
         'font-size'=>'70%','text-anchor'=>'start','font-weight'=>'100'
     }
)->cdata("$tail kb");
=cut
for(my $i=$min;$i<=$max;$i=$i+$step){
     my $tempy=$y1-($i-$min+1)/$rate;
     #print "$tempy\t$min\t$step\t$i\t$rate\n";
     #print "$tempx\t$y\n";
     my $line=$svg->line(
         x1=>$x,y1=>$tempy,
         x2=>$x+5,y2=>$tempy,
         style=>{stroke=>'black'}
     );
     my $tempi=int ($i/1000);
     my $text=$svg->text(
         x=>$x+20,y=>$tempy+3,
         style=>{
             'font-size'=>'70%','text-anchor'=>'end','font-weight'=>'100'
         }
     )->cdata($tempi);
}
return $svg;
}
#####
sub parseGFF
{
my ($gff)=@_;
my %hash;  ##hash to store every record by key of Seq_id
my $seq;   ##Scaffold
my $id;    ##ID for element
my $record;##all line for this record
my $index; ##key, Seq_id
open IN, "$gff" or die "$!";
while (<IN>){
    chomp $_;
    next if ($_=~/^#/);
    my @unit=split("\t",$_);
    if ($unit[2]=~/mRNA/){
        $seq=$unit[0];
        if ($unit[8]=~/ID=(.*?);/ or $unit[8] =~/ID=(.*)/){
            $id=$1;
        }
        $record="$_\n";
        $hash{$id}=$record;
    }elsif($unit[0] eq $seq and $unit[8] =~ /Parent=$id/){
        $hash{$id}.="$_\n";
    }

}
close IN;
return \%hash;
}

#####
sub parseTEGFF
{
my ($gff)=@_;
my %hash;  ##hash to store every record by key of Seq_id
my $seq;   ##Scaffold
my $id;    ##ID for element
my $record;##all line for this record
my $index; ##key, Seq_id
open IN, "$gff" or die "$!";
while (<IN>){
    chomp $_;
    next if ($_=~/^#/);
    my @unit=split("\t",$_);
    if ($unit[8]=~/ID=(.*?);/ or $unit[8] =~/ID=(.*)/){
        $id=$1;
        $hash{$id}="$_";
    }

}
close IN;
return \%hash;
}


################################### sub for write svg to file
sub writesvg {
my ($file,$svg)=@_;
#print "$file\n";
open OUT, ">$file" or die "can not open my file";
       print OUT $svg->xmlify;
close OUT;
       system "/rhome/cjinfeng/software/tools/draw/svg2xxx_release/svg2xxx $file -t pdf -m 700";
}

