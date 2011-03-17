/****************************************************************************************************/
//  COPYRIGHT 2011, University of Tennessee
//  Author: David Jenkins (david.d.jenkins@gmail.com)
//  File: staar.cpp
//  Date: 12 Jan 2011
//  Version: 1.0
//  Description: This is the STAAR program. It will take a PDB file/directory parsing through it
//               looking for possible anion-quadrople interactions.  It will output into a
//               specified file all possible interactions that are within a specified distance 
//               threshold and the angles at which the two residues are interacting at.
//
//  Updates: Added the ability to input residues of interest through the cmd line (11 Feb 2011)
//
/***************************************************************************************************/
//
/***************************************************************************************************/
//  Redistribution and use in source and binary forms, with or without modification,
//  are permitted provided that the following conditions are met:
//  Redistributions of source code must retain the above copyright notice,
//  this list of conditions and the following disclaimer.
//  Redistributions in binary form must reproduce the above copyright notice,
//  this list of conditions and the following disclaimer in the documentation
//  and/or other materials provided with the distribution.
//  Neither the name of the University of Tennessee nor the names of its contributors
//  may be used to endorse or promote products derived from this software
//  without specific prior written permission.
//
//  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
//  "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
//  TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
//  PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR
//  CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
//  EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
//  PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS
//  OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
//  WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
//  OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
//  ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
/*************************************************************************************************/

#include <cstdlib>
#include <cstdio>
#include <iostream>
#include <fstream>
#include <float.h>

#include "Utils.hpp"
#include "Options.hpp"
#include "PDB.hpp"
#include "Seqres.hpp"
#include "Geometry.hpp"
#include "AminoAcid.hpp"
#include "Coordinates.hpp"
#include "CoutColors.hpp"

#define MAX_STR_LENGTH 1024

// Parses through a single PDB file checking for interactions
bool processSinglePDBFile(char* filename,
                          Options& opts,
                          ofstream& output_file);

// Traverses through a directory of PDB files processing each one
bool processPDBDirectory(Options& opts);

// Searches through all of the chains looking for interactions
void searchChainInformation(PDB & PDBfile,
                            unsigned int chain1,
                            unsigned int chain2,
                            string residue1,
                            string residue2,
                            Options & opts,
                            ofstream& output_file);

// Finds the closest distance among all of the centers
// associated with each amino acid
double findClosestDistance(AminoAcid& aa1,
			   AminoAcid& aa2,
			   float threshold,
			   unsigned int* closest_index1,
			   unsigned int* closest_index2);

// Finds the closest interaction among all of the possible 
// amino acid centers
void findBestInteraction( AminoAcid& aa1,
                          AminoAcid& aa2,
                          float threshold,
                          string filename,
			  char* gamessfolder,
                          ofstream& output_file);

// Writes the INP files
void outputINPfile(string input_filename,
		   char* filename, 
		   AminoAcid& aa1h, 
		   AminoAcid& aa2h);

void write_output_head(ofstream& out);

int main(int argc, char* argv[]){
  printHeader();
  int return_value;
  double start = getTime();

  // Parse and error check the command line arguments
  Options opts(argc, argv);

  // Looks like something is wrong with the cmd line args
  // or the user input the -h flag
  if(opts.fail())
    {
      return 1;
    }

  // Lets check to see what kind of input we received for the
  // PDB file.  If a directory, parse all files in the directory.
  // Otherwise just parse the single file
  if( isDirectory(opts.pdbfile) )
    {
      // go through each file in the directory
      return_value = processPDBDirectory(opts);
    }
  else
    {
      ofstream output_file(opts.outputfile);
      write_output_head(output_file);
      return_value = processSinglePDBFile(opts.pdbfile, opts, output_file);
      output_file.close();
    }

  cout << "Time taken: " << getTime() - start << "s" << endl;
  return !return_value;
}

bool processSinglePDBFile(char* filename,
                          Options& opts,
                          ofstream& output_file)
{

  // Read in the PDB file
  // This function actually reads and stores more information than
  // we really need, but the files are relatively small so it isn't
  // taking up much RAM (from what I saw, < 5MB each PDB file)
  // with a few exceptions
  PDB PDBfile(filename);

  if(PDBfile.fail())
    {
      cerr << red << "Error" << reset << ": Parsing PDB file, " << opts.pdbfile << ", failed!" << endl;
      return false;
    }

  PDBfile.populateChains(false);

  if( opts.numLigands )
    {
      PDBfile.findLigands( opts.ligands );
    }

  int numRes1 = opts.residue1.size();
  int numRes2 = opts.residue2.size();

  // Searching for interations within each chain
  for(unsigned int i = 0; i < PDBfile.chains.size(); i++)
    {
      // if we want to look for interactions between the ith chain
      // and each of the other chains, we set the indices to loop
      // though all the chains
      unsigned int start = 0;
      unsigned int end   = PDBfile.chains.size();

      // otherwise we just set the indices to go through the ith chain
      if(opts.sameChain)
        {
          start = i;
          end   = i+1;
        }
      for(unsigned int j = start; j < end; j++)
        {
          for(int ii = 0; ii < numRes1; ii++)
            {
              for(int jj = 0; jj < numRes2; jj++)
                {
                  // Search for different residue combinations
                   searchChainInformation(PDBfile,
                                          i,
                                          j,
                                          opts.residue1[ii],
                                          opts.residue2[jj],
                                          opts,
                                          output_file);
                }
            }
        }
      // for(unsigned int j = 0; j<opts.numLigands; j++)
      // 	{
	  
      // 	}
    }

  return true;
}

bool processPDBDirectory(Options& opts)
{
  DIR* directory;
  struct dirent* filename;
  int numberOfFiles;
  int count=0;
  ofstream output_file(opts.outputfile);
  write_output_head(output_file);

  // Count the number of files in the given PDB directory
  //     Used for checkpointing only
  numberOfFiles = countFilesInDirectory(opts.pdbfile);

  // Open the directory for traversal
  if( (directory = opendir( opts.pdbfile )) )
    {
      // Go through each file of the directory
      while( (filename = readdir( directory )) )
        {
          if( strcmp(filename->d_name, ".") != 0 && strcmp(filename->d_name, "..") != 0 )
            {
              count++;
              // Checkpoint
              cout << purple << "File " << count << "/" << numberOfFiles << " : " << filename->d_name << endl;

              //create the full path to the file
              char fullFilePath [MAX_STR_LENGTH];
              sprintf(fullFilePath, "%s/%s", opts.pdbfile, filename->d_name);

              // perform some work on the current file
              if ( !processSinglePDBFile(fullFilePath, opts, output_file) )
                {
                  closedir(directory);
                  return false;
                }
            }
        }
    }
  else
    {
      perror(opts.pdbfile);
      return false;
    }
  closedir(directory);
  output_file.close();
  return true;
}

void searchChainInformation(PDB & PDBfile,
                            unsigned int chain1,
                            unsigned int chain2,
                            string residue1,
                            string residue2,
                            Options & opts,
                            ofstream& output_file)
{

  Chain* c1 = &(PDBfile.chains[chain1]);
  Chain* c2 = &(PDBfile.chains[chain2]);
#ifdef DEBUG
  unsigned int length_chain1 = c1->seqres[0]->numberOfResidues;
  unsigned int length_chain2 = c2->seqres[0]->numberOfResidues;

  cout << purple << "Checking for res1= "<< residue1 << ", res2= "<< residue2 << endl;
       << " in chains " << c1->id << " and " << c2->id << "..." << endl;

  cout << purple << "length(Chain1)= " << length_chain1
       << " length(Chain2)= " << length_chain2 << endl;
#endif
  // Go through each AA in the first chain
  for(unsigned int i = 0; i < c1->aa.size(); i++)
    {
      // If this chain has a residue that we are looking for,
      // let's do some analysis!
      if( c1->aa[i].residue.compare(residue1) == 0 && !(c1->aa[i].skip))
        {
          for(unsigned int j = 0; j < c2->aa.size(); j++)
            {
              if(c2->aa[j].residue == residue2 && !(c2->aa[j].skip))
                {
                  // Find the best interaction out of all the centers
                  // for this AA pair
                  findBestInteraction( c1->aa[i],
                                       c2->aa[j],
                                       opts.threshold,
                                       PDBfile.filename,
				       opts.gamessfolder,
                                       output_file);
                }
            }
        }
    }
#ifdef DEBUG
  cout << purple << "\tDone" << endl;
#endif
}

// Finds the closest distance among all of the centers
// associated with each amino acid
double findClosestDistance(AminoAcid& aa1,
			   AminoAcid& aa2,
			   float threshold,
			   unsigned int* closest_index1,
			   unsigned int* closest_index2)
{
  double dist;
  double closest = FLT_MAX;

  // Go through all combination of distances looking 
  // for the closet pair
  for( unsigned int i = 0; i < aa1.center.size(); i++ )
    {
      for( unsigned int j = 0; j < aa2.center.size(); j++ )
        {
          dist = aa1.center[i].distance(aa2.center[j]);
          // flag it if is the closest and within the threshold
          if( dist < closest && dist < threshold )
            {
              *closest_index1 = i;
              *closest_index2 = j;
              closest = dist;
            }
        }
    }
  return closest;
}

void findBestInteraction( AminoAcid& aa1,
                          AminoAcid& aa2,
                          float threshold,
                          string input_filename,
			  char* gamessfolder,
                          ofstream& output_file)
{
  float closestDist = FLT_MAX;
  unsigned int closestDist_index1 = 0;
  unsigned int closestDist_index2 = 0;
  float dist;
  float angle;
  float angle1;
  float angleP;
  float angleh;
  float angleOxy;
  float angleOxy2;
  AminoAcid aa1h;
  AminoAcid aa2h;
  char output_filename[1024] = "N/A";
  static int numOutputted = 0;

  // Go through all combination of distances looking 
  // for the closet pair
  closestDist = findClosestDistance(aa1, 
				    aa2, 
				    threshold,
				    &closestDist_index1, 
				    &closestDist_index2);

  // AND WE HAVE A WINNER! (sorta)
  // if we found something output it to the file,
  // add hydrogens, and try the process again
  if( closestDist != FLT_MAX )
    {
      // Just some codes that were in the original STAAR
      char code1 = 'I';
      if( aa1.atom[0]->chainID != aa2.atom[0]->chainID )
        code1 = 'X';
      char code2 = 'S';
      //if( aa1.center.size() > 1 || aa2.center.size() > 1 )
      if( aa1.altLoc || aa2.altLoc )
        code2 = 'M';
      
      PDB pairWithHydrogen;
      pairWithHydrogen.addHydrogensToPair(aa1,aa2);
      
      // The following is just to put the benzene in aa1h
      // and the formate in aa2h just so that I can keep
      // this straight in my head.  This may not be necessary
      // but it helps me when I am looking through the code.

      // If these residues were in different chains
      if( pairWithHydrogen.chains.size() != 1 )
      	{
      	  // If the benzene was naturally first, 
      	  // it is in the first chain while the formate
      	  // is in the second chain
      	  if(aa1.atom[0]->resSeq < aa2.atom[0]->resSeq)
      	    {
      	      aa1h = pairWithHydrogen.chains[0].aa[0];
      	      aa2h = pairWithHydrogen.chains[1].aa[0];
      	    }
      	  // otherwise they are in the opposite order
      	  else
      	    {
      	      aa1h = pairWithHydrogen.chains[1].aa[0];
      	      aa2h = pairWithHydrogen.chains[0].aa[0];
      	    }
      	}
      // If these residues were in the same chain
      else
      	{
      	  // If the benzene was naturally first, 
      	  // it is in the first chain while the formate
      	  // is in the second chain
      	  if(aa1.atom[0]->resSeq < aa2.atom[0]->resSeq)
      	    {
      	      aa1h = pairWithHydrogen.chains[0].aa[0];
      	      aa2h = pairWithHydrogen.chains[0].aa[1];
      	    }
      	  // otherwise they are in the opposite order
      	  else
      	    {
      	      aa1h = pairWithHydrogen.chains[0].aa[1];
      	      aa2h = pairWithHydrogen.chains[0].aa[0];
      	    }
      	}

      // Now we are looking at the distances between the AA's
      // and center of charges of the GLU and ASP
      
      if( closestDist != FLT_MAX )
	{
	  // calculate the angles of this interaction
	  aa1.calculateAnglesPreHydrogens(aa2,
	  				  closestDist_index1,
	  				  closestDist_index2,
	  				  &angle,
	  				  &angle1,
	  				  &angleP);
	  
	  float dist;
	  float distOxy;
	  float distOxy2;
	  aa1h.calculateDistancesAndAnglesPostHydrogens(aa2h,
							aa2.center[closestDist_index2],
							&dist,
							&distOxy,
							&distOxy2,
							&angleh,
							&angleOxy,
							&angleOxy2);
	  
	  // Here, we are outputting the 
	  if( gamessfolder )
	    {
	      numOutputted++;
	      sprintf(output_filename, "%s/gamessinp-%d.inp", gamessfolder, numOutputted);
	      outputINPfile(input_filename, output_filename, aa1h, aa2h);
	    }
	  else 
	    {
	      sprintf(output_filename, "N/A");
	    }
	  
	  // and we finally output some results!
	  output_file << aa1.residue                        << "\t"
		      << aa2.residue                        << "\t"
		      << closestDist                        << "\t"
		      << angle                              << "\t"
		      << angleP                             << "\t"
		      << angle1                             << "\t"
		      << aa1.atom[0]->resSeq                << "\t"
		      << aa2.atom[0]->resSeq                << "\t"
		      << code1 << code2                     << "\t"
		      << input_filename                     << "\t"
		      << output_filename                    << "\t"
		      << aa1.atom[0]->chainID               << "\t"
		      << aa2.atom[0]->chainID               << "\t"
		      << aa1.center[closestDist_index1]     << "\t"
		      << aa2.center[closestDist_index2]     << "\t" 
		      << aa1h.center[0]                     << "\t" 
		      << aa2h.center[0]                     << "\t"
		      << dist                               << "\t"
		      << distOxy                            << "\t"
		      << distOxy2                           << "\t "
		      << angleh                             << "\t"
		      << angleOxy                           << "\t "
		      << angleOxy2                          << "\t" << endl;
	}      
    }
}

void outputINPfile(string input_filename, char* filename, AminoAcid& aa1h, AminoAcid& aa2h)
{
  ofstream inpout(filename);

  inpout << INPheader << endl;
  inpout << " $MOROKM IATM(1)=" << aa1h.atom.size() << "," << aa2h.atom.size() << " ICHM(1)=0,-1" << " $END" << endl;
  inpout << " $DATA" << endl;
  inpout << input_filename << endl;
  inpout << "C1" << endl;
  for(int i=0; i<aa1h.atom.size(); i++)
    {
      if( aa1h.atom[i]->element == " H" )
	{
	  inpout << "H      1.0     ";
	}
      else if( aa1h.atom[i]->element == " C")
	{
	  inpout << "C      6.0     ";
	}
      else if( aa1h.atom[i]->element == " O")
	{
	  inpout << "O      8.0     ";
	}
      inpout << aa1h.atom[i]->coord << endl;
    }
  
  for(int i=0; i<aa2h.atom.size(); i++)
    {
      if( aa2h.atom[i]->element == " H")
	{
	  inpout << "H      1.0     ";
	}
      else if( aa2h.atom[i]->element == " C")
	{
	  inpout << "C      6.0     ";
	}
      else if( aa2h.atom[i]->element == " O")
	{
	  inpout << "O      8.0     ";
	}
      inpout << aa2h.atom[i]->coord << endl;
    }
  
  inpout << " $END" << endl;
  inpout.close();
}

void write_output_head(ofstream& out)
{
  out <<"# res1  res2    dist    angle   angleP  angle1  loc1    loc2    code    pdbID                           gamessinput             center1                 center2                 center1h                center2h                dist    distOxy     distOxy2 angleh         angleOxy angleOxy2" << endl;
}
