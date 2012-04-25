#!/usr/bin/perl
#
# COPYRIGHT (C)  2011,  University of Tennessee
# Author:  David Jenkins (david.d.jenkins@gmail.com)
#
# File:  get_energies.pl
# Date: 10 Jun 2011
# Description: 
# 

if($#ARGV+1 != 1)
{
    print STDERR "get_energies.pl STAAR.csv\n";
    exit;
}

my $staar_file = $ARGV[0];
open(fp, "<$staar_file");
my @staarout = <fp>;
close(fp);

my $num_lines=scalar(@staarout);

for(my $i=0;$i<$num_lines;$i++)
{
    if( substr($staarout[$i],0,1) eq "#"){ next; }
    my @line = split(',',$staarout[$i]);
    my $pdb = $line[9];
    my $res1 = $line[0];
    my $res2 = $line[1];
    my $loc1 = $line[6];
    my $loc2 = $line[7];
    my $chain1 = $line[13];
    my $chain2 = $line[14];
    my $filename = $line[12];
    my $model = $line[11];
    my $resolution = $line[10];
    my $dist = $line[27];
    my $angle = $line[30];

    # Pull out the file name number
    my $index_of_slash = rindex($filename,'/');
    my $out_directory = substr($filename,0,$index_of_slash);
    $out_directory =~ s/inp/out/g;

    $filename =~ m/^[\w\/]*gamessinp-(\d+).inp/;
    my $fnumber = $1;
    $pdb =~ m/^[\w\/]*([\d\w]{4}?).pdb.gz/;
    $pdb = $1;

    my $norm=`grep " FINAL RHF" $out_directory/gamessout-$fnumber-g1.out | awk '{print \$5}'`;
    my $benz=`grep " FINAL RHF" $out_directory/gamessout-$fnumber-g2.out | awk '{print \$5}'`;
    my $form=`grep " FINAL RHF" $out_directory/gamessout-$fnumber-g3.out | awk '{print \$5}'`;
    my $totalh = ($norm-$benz-$form);
    my $totalk = $totalh * 627.509;
    print $pdb.",".$resolution.",".$model.",".$res1.",".$loc1.",".$chain1.",".$res2.",".$loc2.",".$chain2.",".$dist.",".$angle.",".$totalh.",".$totalk."\n";
}

