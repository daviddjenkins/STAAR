#STAAR: __ST__atistical __A__nalysis of __A__romatic __R__ings
Version 2.0; October 2012
Copyright (C) 2012 University of Tennessee, Knoxville
-----------------------------------------------------------------------

STAAR searches a set of PDB files for possible anion-quadrupole 
interactions.

Other files to read in the top-level source directory:

    INSTALL           Installation instructions.
    LICENSE           The GNU General Public License (GPLv3).

STAAR uses the gzstream library for the decompression of PDB files that
are compressed.  The library source code can be found in its entirety
in the gzstream folder.  I did have to make a change to the code 
because the istream '!' operator was not working correctly to check
for a valid open. What I did was to add a failure flag to one of the 
classes and a get function named fail().  It returns true if the file
could not be opened. The gzstream library is freely available under
the LGPLv2.1 library that can be seen in gzstream/COPYING.lib.  For 
more information on the gzstream library, please visit their page at:
    http://www.cs.unc.edu/Research/compgeom/gzstream/

The STAAR code was not intended to be library of code that you could 
plug into to develop more code.  This is evident by the amount of 
hoops you would have to jump through just to compile it.  This is
mainly due to the OpenBabel dependency.  However, if you wish to use
it as a library, you should be able to as we have done it a couple
of times.  The only thing that is truly specifc to STAAR are the 
Options.{cpp,hpp} files, but from my recollection nothing else
depends on those files.  There are limited options in the PDB reader
as we were only concerned with the ATOM lines (well, for one of our
side projects, HETATM lines). If you need support for other PDB 
features, you will need to edit PDB.{cpp,hpp}.  Anyway, good luck.

The scripts folder contains some scripts and codes for parsing the 
GAMESS output, added the results to the STAAR output, dealing with
some of the database stuff, etc. Speaking on this folder, here are 
some explanations of the codes in there [Disclaimer: they weren't
created with general applicability in mind. They were simply created
for our own specific needs, but we are including them in case you
might need them for reference]:
     downloadPDB - This will download a supplied list of PDBs from 
                    the RCSB website
     getNewPDBFiles - Similar to downloadPDB, but it will only 
                      download whatever hasn't already been 
                      downloaded before
     make_scripts - Makes the submission scripts and submits the jobs
                    to the Newton cluster.  This is specific to our 
                    needs, so it may not help others much. Check the
                    comments inside the script for more explanations
     parseGamess.pl - Perl script to parse out the information from 
                      the GAMESS runs.  Again, this is somewhat 
                      specific to our needs, but you can use it as a 
                      guideline, if you wish
     splitInput - Splits a large list of PDBs into multiple smaller 
                  ones.  Again, specific to our needs
     stats - Gets stats about our runs.  This one is actually one you
             can use on your own, assuming you don't have results 
             split up. 

If you would rather not download and compile this code on your own 
system, please visit staar.bio.utk.edu where you can upload your zipped
PDB file and run through the analysis on our servers.  This website 
already includes the results for all PDB files from 9 May 2012.

If you choose to use this code, website, or results produced from this 
work, please cite the following paper:

> Jenkins, D. D., Harris, J. B., Howell, E. E., Hinde, R. J. and Baudry, J. 
> (2012), STAAR: Statistical analysis of aromatic rings. J. Comput. Chem.. 
> doi: [10.1002/jcc.23164](http://dx.doi.org/10.1002/jcc.23164)

-----------------------------------------------------------------------

