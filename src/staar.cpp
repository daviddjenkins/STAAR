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

// Finds the closest interaction among all of the possible 
// amino acid centers
void findBestInteraction( AminoAcid& aa1,
                          AminoAcid& aa2,
                          float threshold,
                          string filename,
                          ofstream& output_file);

// Calculates all of the angles that are outputted to a file
void calculateAngles(AminoAcid& aa1,
                     AminoAcid& aa2,
                     int index1,
                     int index2,
                     float* angle,
                     float* angle1,
                     float* angleP);

void write_output_head(ofstream& out);

int main(int argc, char* argv[]){

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
  PDB PDBfile(filename, opts);

  if(PDBfile.fail())
    {
      cerr << "Parsing PDB file, " << opts.pdbfile << ", failed!" << endl;
      return false;
    }

  PDBfile.populateChains(opts);

  int numRes1 = opts.residue1.size();
  int numRes2 = opts.residue2.size();

  // Set the residue combinations to look for.
  // Possibly add this as a cmd line arg later for more a generic program
  //string residue1[] = {"TYR", "PHE", "TRP"};
  //string residue2[] = {"GLU", "ASP"};

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
              cout << "File " << count << "/" << numberOfFiles << " : " << filename->d_name << endl;

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
  // unsigned int length_chain1 = c1->seqres[0]->numberOfResidues;
  // unsigned int length_chain2 = c2->seqres[0]->numberOfResidues;

  // cout << "Checking for res1= "<< residue1 << ", res2= "<< residue2
  //      << " in chains " << c1->id << " and " << c2->id << "..." << endl;

  // cout << "length(Chain1)= " << length_chain1
  //      << " length(Chain2)= " << length_chain2 << endl;

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
                                       output_file);
                }
            }
        }
    }

  // cout << "\tDone" << endl;
}

void findBestInteraction( AminoAcid& aa1,
                          AminoAcid& aa2,
                          float threshold,
                          string input_filename,
                          ofstream& output_file)
{
  float closestDist = FLT_MAX;
  int closestDist_index1 = 0;
  int closestDist_index2 = 0;
  float dist;
  float angle;
  float angle1;
  float angleP;

  // Go through all combination of distances looking 
  // for the closet pair
  for( unsigned int i = 0; i < aa1.center.size(); i++ )
    {
      for( unsigned int j = 0; j < aa2.center.size(); j++ )
        {
          dist = aa1.center[i].distance(aa2.center[j]);
          // flag it if is the closest and within the threshold
          if( dist < closestDist && dist < threshold )
            {
              closestDist_index1 = i;
              closestDist_index2 = j;
              closestDist = dist;
            }
        }
    }

  // if we found something output it to the file
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

      // calculate the angles of this interaction
      calculateAngles(aa1,
                      aa2,
                      closestDist_index1,
                      closestDist_index2,
                      &angle,
                      &angle1,
                      &angleP);

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
                  << aa1.atom[0]->chainID               << "\t"
                  << aa2.atom[0]->chainID               << "\t"
                  << aa1.center[closestDist_index1]     << "\t"
                  << aa2.center[closestDist_index2]     << endl;
    }
}

void calculateAngles(AminoAcid& aa1,
                     AminoAcid& aa2,
                     int index1,
                     int index2,
                     float* angle,
                     float* angle1,
                     float* angleP )
{
  Coordinates planeP;
  Coordinates planeProject;
  
  // Calculate the equation of the plane
  float det = getPlaneEquation( *(aa1.center[index1].plane_info[ CG_PLANE_COORD_PTT]),
                                *(aa1.center[index1].plane_info[CD1_PLANE_COORD_PTT]),
                                *(aa1.center[index1].plane_info[CD2_PLANE_COORD_PTT]),
                                &planeP);

  // Calculate the angle between a plane and a line
  *angle = angleBetweenPlaneAndLine( planeP,
                                     aa1.center[index1],
                                     aa2.center[index2]);
  
  // Calculate the "plane project" coordinates
  planeProjectCoordinate( planeP, 
                          aa2.center[index2],
                          -1 * det, 
                          &planeProject);
  
  // find the second angle betwen the CG ATOM and the "plane project"
  *angle1 = findAngle(*aa1.center[index1].plane_info[ CG_PLANE_COORD_PTT],
                      aa1.center[index1],
                      planeProject);
  
  // finally, calculate the angle between the planes of AA1 and AA2
  *angleP = calculateAngleBetweenPlanes( planeP,
                                         aa2,
                                         index2 );

}

void write_output_head(ofstream& out)
{
  out <<"res1\tres2\tdist\tangle\tangleP\tangle1\tloc1\tloc2\tcode\tpdbID\tchain1\tchain2" << endl;
}
