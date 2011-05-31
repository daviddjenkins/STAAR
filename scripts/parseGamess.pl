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
if($#ARGV + 1 != 2)
{
    print STDERR "Usage: perl parseGamess.pl gamess_output_directory STAAR.out\n";
    exit;
}

# Save the command line args
my $output_directory = $ARGV[0];
my $STAARout = $ARGV[1];

# Open the STAAR.out table, read it, then close it
open(fp, "<$STAARout");
my @STAARtable = <fp>;
close(fp);

my $number_of_pairs = scalar @STAARtable;

# Prints out the header of the table
print @STAARtable[0];

# Go through all the potential pairs in the STAAR.out table
for( $i=1; $i<$number_of_pairs; $i++ )
{
    # Get the line and split it into comma separated fields
    my $pair_line;
    chomp($pair_line = @STAARtable[$i]);
    my @fields = split(/,/,$pair_line);

    # get the corresponding output filename
    my $out_filename = &getOutputFilename(@fields[12],$output_directory);

    # Check if this is a PO4 line, if it is, we are going to look for the one
    # with the most negative total energy
    if( @fields[1] eq "PO4" )
    {
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
    else
    {
        # Get the energies
        my $energies = &getEnergies($out_filename);

        # Check the validity of the results
        my ($skipped,$total) = &checkForInvalidResults($energies, @fields[12]);

        # If the results are good, let's print them out
        if( not $skipped )
        {
            print $pair_line.','.$out_filename.$energies."\n";
        }
    }
}

print STDERR "PO4 count: ".$count."\n";
print STDERR "PO4 not skipped: ".$count2."\n";
print STDERR "PO4 saved: ".$count3."\n";

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
    # just in case we get curious about the details of the potential pair
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
