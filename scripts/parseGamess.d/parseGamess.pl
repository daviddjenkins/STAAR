#!/usr/bin/perl
#
# COPYRIGHT (C)  2011,  University of Tennessee
# Author:  David Jenkins (david.d.jenkins@gmail.com)
#
# File:  parseGamess.pl
# Date: 03 May 2011
# Description: Parses the GAMESS output for the MK energies
#              and appends them to the STAAR.out table
#
# Assumptions: Assumes input files are gamessinp-*.inp and
#              output files are gamessout-*.out 
#              If named otherwise, edit getOutputFilename
#              function

# Check for arguments
if($#ARGV + 1 != 3)
{
    print STDERR "Usage: perl parseGamess.pl gamess_output_directory STAAR.out output_full_or_short\n";
    exit;
}

# Save the command line args
my $output_directory = $ARGV[0];
my $STAARout = $ARGV[1];
my $output_type = $ARGV[2];
if(!($output_type eq "full") && !($output_type eq "short"))
{
    print STDERR "Output type must be either \"full\" or \"short\"\n";
    exit;
}

# Open the STAAR.out table, read it, then close it
open(fp, "<$STAARout");
my @STAARtable = <fp>;
close(fp);

my $number_of_pairs = scalar @STAARtable;

# Prints out the header of the table
if($output_type eq "full")
{
    print @STAARtable[0];
}
else
{
    print "#PDB,Resolution,Model,Residue1,Loc1,Chain1,Residue2,Loc2,Chain2,Distance,Angle,Energy(Hartree),Energy(kcal/mol)\n";
}

# Remove all stuff
chomp($pair_line = @STAARtable[1]);
my @fields = split(/,/,$pair_line);
my $index_of_slash = rindex($fields[12],'/');
my $in_directory = substr($fields[12],0,$index_of_slash);

if ( -e "$in_directory/runs_redo.sge")
{
    unlink("$in_directory/runs_redo.sge")
}

if( -e "$STAARout.redo")
{
    unlink("$STAARout.redo");
}


# Go through all the potential pairs in the STAAR.out table
for( $i=0; $i<$number_of_pairs; $i++ )
{
    my $pair_line;
    chomp($pair_line = @STAARtable[$i]);

    if( substr($pair_line,0,1) eq "#"){ next; }
    # Get the line and split it into comma separated fields
    my @fields = split(/,/,$pair_line);

    # get the corresponding output filename
    my $out_filename = &getOutputFilename(@fields[12],$output_directory);

    # Check if this is a PO4 line, if it is, we are going to look for the one
    # with the most negative total energy
    if( @fields[1] eq "PO4" )
    {
        deal_with_PO4($output_directory, $output_filename, @STAARtable);
    }
    else
    {
        # Get the energies
        my $energies = &getEnergies($out_filename);

        # check the validity of the results
        my ($skipped,$total) = &checkForInvalidResults($energies, @fields[12]);

        # If the results are good, let's print them out
        if( not $skipped )
        {
            # Check for big results
            $skipped = &checkForBigResults($total, $fields[12], $STAARout.".redo", $pair_line, $output_directory);
            if( not $skipped )
            {
                my $pdb = $fields[9];
                my $res1 = $fields[0];
                my $res2 = $fields[1];
                my $loc1 = $fields[6];
                my $loc2 = $fields[7];
                my $chain1 = $fields[13];
                my $chain2 = $fields[14];
                my $filename = $fields[12];
                my $model = $fields[11];
                my $resolution = $fields[10];
                my @efields=split(/,/,$energies);
                my $totalh = $efields[11];
                my $totalk = $efields[12];
                my $dist = $fields[27];
                my $angle = $fields[30];

                $pdb =~ m/^[\w\/]*([\d\w]{4}?).pdb.gz/;
                $pdb = $1;
                
                if($output_type eq "full")
                {
                    print $pair_line.','.$out_filename.$energies."\n";
                }
                else
                {
                    print $pdb.",".$resolution.",".$model.",".$res1.",".$loc1.",".$chain1.",".$res2.",".$loc2.",".$chain2.",".$dist.",".$angle.",".$totalh.",".$totalk."\n";
                }
            }
        }
    }
}

if ( -e "$in_directory/runs_redo.sge")
{
    `qsub $in_directory/runs_redo.sge`;
}

# Create the output filename.  This is dependent on the way that I have 
# files named. Assumes input files are gamessinp-*.inp and output files
# are gamessout-*.out 
sub getOutputFilename
{
    my($in_filename)      = @_[0];
    my($output_directory) = @_[1];

    # Get the index of the last slash if it exists
    my $index_of_slash = rindex($in_filename,'/');

    # Create output file
    my $out_filename = $output_directory.substr($in_filename, $index_of_slash);

    # replace all inp with out 
    $out_filename =~ s/inp/out/g;
    return $out_filename;
}

# Parse the file with some magical bash one liner.  This seems magical, but
# what it does is greps for all the energies and places them all on the same
# line each separated by commas to be added to the STAAR.out table later.
# This should not ever need to change unless GAMESS changes their output. If
# they do, good luck! *insert evil laugh here*
sub getEnergies
{
    my $out_filename = @_[0];
    return `egrep "ELECTROSTATIC ENERGY             ES=|EXCHANGE REPULSION ENERGY        EX=|POLARIZATION ENERGY              PL=|CHARGE TRANSFER ENERGY           CT=|HIGH ORDER COUPLING ENERGY      MIX=|TOTAL INTERACTION ENERGY,   DELTA-E=" $out_filename | tr -s " " | cut -d'=' -f2 | tr ' ' ',' | tr -d '\n'`;
}

# Check the validity of the results
sub checkForInvalidResults
{
    my $results  = @_[0];
    my $filename = @_[1];
    my @fields;
    # There are no results.  happens with a handful of cases where a hydrogen is not
    # added to the formate in GLU or ASP residues when the distances are too far to 
    # determine bonds correctly.  This is a problem with Babel, not STAAR.  Although,
    # I could put a check for this in STAAR, but I have decided to output it anyway 
    # just in case we get curious about the details of the potential pair or when
    # Newton is just being crappy for no reason at all
    if( !$results )
    {
        print STDERR "Skipping ".$filename." : No results\n";
        return (1,0);
    }
    # This checks for GAMESS runs that didn't converge as indicated by ******** in 
    # the energy results
    elsif( $results =~ /\*{8}/ )
    {
        print STDERR "Skipping ".$filename." : Didn't converge within 200 iterations\n";
        return (1,0);
    }
    else
    {
        # Now, we know how have potentially good results, but we need to check if 
        # they are valid.  If the absolute_value(mix energy) > .25, we are skipping 
        # it.  Other wise, we will return from this function successfully!
        @fields = split(/,/,$results);
        if( abs(@fields[10]) > .25 )
        {
            print STDERR "Skipping ".$filename." : abs( MIX energy ) > .25\n";
            return (1,0);
        }
    }
    return (0,@fields[12]);
}

sub checkForBigResults
{
    my($total, $in_filename, $redo_file, $pair_line, $out_directory) = @_;
    if($total <= -24)
    {
        print STDERR "Skipping ".$in_filename." : No results\n";        
#         `bash /lustre/AQ/AllPDB/split_input_file $in_filename`;
#         open($fp,">>$redo_file");
#         print $fp $pair_line."\n";
#         close($fp);
#         my $index_of_slash = rindex($in_filename,'/');

#         # Create output file
#         my $in_directory = substr($in_filename,0,$index_of_slash);
#         $in_filename =~ m/^[\w\/]*gamessinp-(\d+).inp/;
#         my $fnumber = $1;

#         # print out job script information
#         if ( ! -e "$in_directory/runs_redo.sge") {
#             open($fp,">$in_directory/runs_redo.sge");
#             print $fp "#\$ -N Gam_redo$fnumber
# #\$ -j y
# #\$ -o /dev/null
# #\$ -q *
# #\$ -pe openmpi* 1
# #\$ -l dedicated=4
# #\$ -cwd
# module load gamess
# cd \$TMPDIR\n";
#         } else {
#             open($fp,">>$in_directory/runs_redo.sge");
#         }

#         print $fp "cp $in_directory/gamessinp-$fnumber-g*.inp .\n";
#         print $fp "mpirun -np 4 /data/AQ/bin/rungms gamessinp-$fnumber-g1.inp 01 1 > $out_directory/gamessout-$fnumber-g1.out
# mpirun -np 4 /data/AQ/bin/rungms gamessinp-$fnumber-g2.inp 01 1 > $out_directory/gamessout-$fnumber-g2.out
# mpirun -np 4 /data/AQ/bin/rungms gamessinp-$fnumber-g3.inp 01 1 > $out_directory/gamessout-$fnumber-g3.out
# ";
#         close($fp);
        return 1;
    }
    return 0;
}

sub dealWithPO4
{
    my ($output_directory, $output_filename, @STAARtable) = @_;
    my ($best_line,$best_energy);
    for(my $j=0; $j<3; $j++, $i++)
    {
        my $po4_line;
        chomp($po4_line = @STAARtable[$i]);
        my @po4fields = split(/,/,$po4_line);

        if(@po4fields[1] != "PO4")
        {
            print STDERR "Failed to handle the PO4 correctly\n";
            exit;
        }

        # get the corresponding output filename
        my $out_filename = &getOutputFilename(@po4fields[12],$output_directory);
        
        # Get the energies
        my $energies = &getEnergies($out_filename);
        
        # Check the validity of the results
        my ($skipped,$total) = &checkForInvalidResults($energies, @po4fields[12]);
        
        # If the results are good, let's print them out
        if( not $skipped )
        {
            if($j == 0)
            {
                $best_line   = $po4_line.','.$out_filename.$energies."\n";
                $best_energy = $total
            }
            elsif( $best_energy > $total )
            {
                print STDERR "Skipping ".($i-1)." : Another PO4 result had better energy\n";
                $best_line   = $po4_line.','.$out_filename.$energies."\n";
                $best_energy = $total
            }
            else
            {
                print STDERR "Skipping ".$i." : Another PO4 result had better energy\n";
            }

        }
    }
    if( $best_line )
    {
        print $best_line;
    }
    $i--;    
}
