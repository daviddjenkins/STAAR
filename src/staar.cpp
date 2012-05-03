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

#include <map>
#include <utility>
#include <sstream>

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
bool processSinglePDBFile(const char* filename,
                          Options& opts,
                          ofstream& output_file,
                          ofstream& pairlistfile,
                          const char* chains=NULL);

// Read PDB names from a list and parses them from the specified directory
bool processPDBList(Options& opts);

// Reads the PDB names and chains from a list and parses them from the
// specified directory
bool processPDBChainList(Options& opts);

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

void searchTripletInformation(PDB & PDBfile,
                            Options& opts,
                            unsigned int chain1,
                            unsigned int chain2,
                            unsigned int chain3,
                            ofstream& output_file,
                            ofstream& pairlistfile);

string buildKeyString(PDB & PDBfile,
                      AminoAcid &aa1,
                      AminoAcid &aa2);

void searchLigandsInformation(PDB & PDBfile,
                              Residue & ligand,
                              unsigned int chain1,
                              string residue1,
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
                          PDB& PDBfile,
                          char* gamessfolder,
                          bool ligand,
                          ofstream& output_file);

// Writes the INP files
void outputINPfile(string input_filename,
                   char* filename,
                   AminoAcid& aa1h,
                   AminoAcid& aa2h);

void write_output_head(ofstream& out,bool triplets);

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
  // PDB file.  If we have a list file of PDBs, just go through
  // list with the specified directory. If a directory, parse 
  // all files in the directory. Otherwise just parse the single 
  // file
  if( opts.pdblist )
    {
      return_value = processPDBList(opts);
    }
  else if( opts.chain_list )
    {
      return_value = processPDBChainList(opts);
    }
  else if( isDirectory(opts.pdbfile) )
    {
      // go through each file in the directory
      return_value = processPDBDirectory(opts);
    }
  else
    {
      ofstream output_file(opts.outputfile);
      if(!output_file)
        {
          cerr << red << "Error" << reset << ": Failed to open output file," << opts.outputfile << endl;
          return 1;
        }
      write_output_head(output_file,opts.triplets);
      //open the pair list file (if we're searching for triplets, and if the option to log pairs is specified
      ofstream pairlistfile;

      if(opts.triplets && opts.pairlistfile != NULL){
        pairlistfile.open(opts.pairlistfile);
        if( !pairlistfile)
          {
            cerr << red << "Error" << reset << ": Failed to open pairfile output file" << endl;
            perror("\t");
          }
      }

      return_value = processSinglePDBFile(opts.pdbfile, opts, output_file,pairlistfile);
      if(opts.triplets && opts.pairlistfile != NULL){
        pairlistfile.close();
      }
      output_file.close();
    }

  cout << "Time taken: " << getTime() - start << "s" << endl;
  return !return_value;
}

bool processSinglePDBFile(const char* filename,
                          Options& opts,
                          ofstream& output_file,
                          ofstream& pairlistfile,
                          const char* chains)
{
  int numRes1 = opts.residue1.size();
  int numRes2 = opts.residue2.size();

  // Read in the PDB file
  // This function actually reads and stores more information than
  // we really need, but the files are relatively small so it isn't
  // taking up much RAM (from what I saw, < 5MB each PDB file)
  // with a few exceptions
  PDB PDBfile_whole(filename, opts.resolution);

  if( PDBfile_whole.fail() )
    {
      //cerr << red << "Error" << reset << ": Parsing PDB file failed!" << endl;
      PDBfile_whole.printFailure();
      return false;
    }

  for(unsigned int model=0; model < PDBfile_whole.models.size(); model++)
    {
      PDB PDBfile = PDBfile_whole.models[model];
      PDBfile.filename = PDBfile_whole.filename;
      PDBfile.resolution = PDBfile_whole.resolution;

      PDBfile.setResiduesToFind(&opts.residue1, &opts.residue2);
      if(opts.numLigands)
        {
          PDBfile.setLigandsToFind(&opts.ligands);
        }
      PDBfile.populateChains(false);

      if( opts.numLigands )
        {
          PDBfile.findLigands( opts.ligands );
        }

      //If we're searching for triplets, then do that here
      if(opts.triplets == true){
      //this triple for loop is set up so that it searches through each
      //unique combination of chains. e.g. if there are 3 chains labeled 1, 2, and 3
      //it will search through the combinations: 111 112 113 122 123 133 222 223 233 333
        for(unsigned int i = 0; i < PDBfile.chains.size(); i++)
        {

	if(opts.sameChain){
          searchTripletInformation(PDBfile,opts,i,i,i,output_file,pairlistfile);
          return true;
	}

	// in the future we may have the below double for loop go through each chain combination.
          for(unsigned int j = i; j < PDBfile.chains.size(); j++)
          {
            for(unsigned int k = j; k < PDBfile.chains.size(); k++)
            {
		searchTripletInformation(PDBfile,opts,i,j,k,output_file,pairlistfile);
            }
          }
        }
        return true;  //end here if we're processing triplets
      }



      // Searching for interations within each chain
      for(unsigned int i = 0; i < PDBfile.chains.size(); i++)
        {
          // Check if we are supposed to look at certain chains
          if( chains )
            {
              if( !strchr(chains,PDBfile.chains[i].id) )
                {
                  continue;
                }
            }
          // if we want to look for interactions between the ith chain
          // and each of the other chains, we set the indices to loop
          // through all the chains
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
      
          // Go through the ligands, if there are any
          for(unsigned int j = 0; j<PDBfile.ligands.size(); j++)
            {
              for(int ii=0; ii<numRes1; ii++)
                {
                  searchLigandsInformation(PDBfile,
                                           *PDBfile.ligands[j],
                                           i,
                                           opts.residue1[ii],
                                           opts,
                                           output_file);
                }
            }
        }
    }
  return true;
}

bool processPDBList(Options& opts)
{
  // Open the list file
  ifstream listfp(opts.pdblist);
  if( !listfp )
    {
      cerr << red << "Error" << reset << ": Failed to open list file, " << opts.pdblist << endl;
      perror("\t");
      return false;
    }

  // Open the output file
  ofstream output_file(opts.outputfile);
  if( !output_file )
    {
      cerr << red << "Error" << reset << ": Failed to open output file" << endl;
      perror("\t");
    }

  // Write the header to output file
  write_output_head(output_file,opts.triplets);

  string line;

  //open the pair list file (if we're searching for triplets, and if the option to log pairs is specified
  ofstream pairlistfile;

  if(opts.triplets && opts.pairlistfile != NULL){
    pairlistfile.open(opts.pairlistfile);
    if( !pairlistfile)
      {
        cerr << red << "Error" << reset << ": Failed to open pairfile output file" << endl;
        perror("\t");
      }
  }

  // Go through each line of the PDB list file
  while(getline(listfp, line))
    {
      // Get and print the file name
      string filename(opts.pdbfile);
      filename += "/" + line + opts.extension;
      cout << purple << line << opts.extension << endl;

      // Process the PDB file
      processSinglePDBFile(filename.c_str(), opts, output_file, pairlistfile);
    }

  cout << endl;
  output_file.close();
  listfp.close();
  if(opts.triplets && opts.pairlistfile != NULL){
    pairlistfile.close();
  }

  return true;
}

bool processPDBChainList(Options& opts)
{
  // Open up the chain list file
  ifstream listfp(opts.chain_list);
  if( !listfp )
    {
      cerr << red << "Error" << reset << ": Failed to open list file, " << opts.chain_list << endl;
      perror("\t");
      return false;
    }

  // Open up the output file
  ofstream output_file(opts.outputfile);
  if( !output_file )
    {
      cerr << red << "Error" << reset << ": Failed to open output file" << endl;
      perror("\t");
    }

  // Write the header to the ouput file
  write_output_head(output_file,opts.triplets);

  string line;

  //open the pair list file (if we're searching for triplets, and if the option to log pairs is specified
  ofstream pairlistfile;

  if(opts.triplets && opts.pairlistfile != NULL){
    pairlistfile.open(opts.pairlistfile);
    if( !pairlistfile)
      {
        cerr << red << "Error" << reset << ": Failed to open pairfile output file" << endl;
        perror("\t");
      }
  }

  // Go through each line of the list file
  while(getline(listfp, line))
    {
      // Store the PDB directory path first
      string filename(opts.pdbfile);

      // Separate the line into files
      vector<string> fields = split(line,'\t');
      if( fields.size() < 2)
        {
          cerr << "Chain file is malformed!" << endl;
          return false;
        }

      // Create and output the file name
      filename += "/" + fields[0] + opts.extension;
      cout << purple << fields[0] << opts.extension << endl;

      // Process the file
      processSinglePDBFile(filename.c_str(), opts, output_file,pairlistfile, fields[1].c_str());
    }
  cout << endl;
  output_file.close();
  listfp.close();
  if(opts.triplets && opts.pairlistfile != NULL){
    pairlistfile.close();
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
  if( !output_file )
    {
      cerr << red << "Error" << reset << ": Failed to open output file" << endl;
      perror("\t");
    }
  write_output_head(output_file,opts.triplets);

  //open the pair list file (if we're searching for triplets, and if the option to log pairs is specified
  ofstream pairlistfile;

  if(opts.triplets && opts.pairlistfile != NULL){
    pairlistfile.open(opts.pairlistfile);
    if( !pairlistfile)
      {
        cerr << red << "Error" << reset << ": Failed to open pairfile output file" << endl;
        perror("\t");
      }
  }

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
              cout << purple << filename->d_name << endl;

              //create the full path to the file
              char fullFilePath [MAX_STR_LENGTH];
              sprintf(fullFilePath, "%s/%s", opts.pdbfile, filename->d_name);

              // perform some work on the current file
              processSinglePDBFile(fullFilePath, opts, output_file,pairlistfile);
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

  if(opts.triplets && opts.pairlistfile != NULL){
    pairlistfile.close();
  }
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
                                       PDBfile,
                                       opts.gamessfolder,
                                       false,
                                       output_file);
                }
            }
        }
    }
#ifdef DEBUG
  cout << purple << "\tDone" << endl;
#endif
}

void searchTripletInformation(PDB & PDBfile,
                            Options& opts,
                            unsigned int chain1,
                            unsigned int chain2,
                            unsigned int chain3,
                            ofstream& output_file,
                            ofstream& pairlistfile)
{
  Chain* c1 = &(PDBfile.chains[chain1]);
  Chain* c2 = &(PDBfile.chains[chain2]);
  Chain* c3 = &(PDBfile.chains[chain3]);

  map<string,int> pairmap;
  map<string,int>::iterator pairmapit;

  double threshold = 7.0;
  int residueType[3];

  // Go through each AA combination in the chain
  for(unsigned int i = 0; i < c1->aa.size(); i++)
  {
    residueType[0] = c1->aa[i].getType();
    if(residueType[0] == AATYPE_UNKNOWN)continue;

    //if we are searching the same chain then we do not want to go over the same amino acid combination twice.
    unsigned int j;
    if(chain1 == chain2)j = i+1;
    else j=0;

    for(; j < c2->aa.size(); j++)
    {
      residueType[1] = c2->aa[j].getType();
      if(residueType[1] == AATYPE_UNKNOWN)continue;

      unsigned int k;
      if(chain2 == chain3)k = j+1;
      else k = 0;

      for(; k < c3->aa.size(); k++)
      {
        float angle[2];
        float angle1[2];
        float angleP[2];
	//init these to -1, so -1 will be output for cases where this isn't calculated.
	angle[0] = angle[1] = angle1[0] = angle1[1] = angleP[0] = angleP[1] = -1.0;
        
        residueType[2] = c3->aa[k].getType();
	if(residueType[2] == AATYPE_UNKNOWN)continue;

        int aaCount[AATYPE_MAX]; //count each of the amino acid types, accumulate them here

	for(int kk=0;kk<AATYPE_MAX;kk++)aaCount[kk]=0;

	for(int kk=0;kk<3;kk++)aaCount[residueType[kk]] += 1;

        //if we have no PI's, then we don't care
        if(aaCount[AATYPE_PI] == 0){
          continue;
        }

        unsigned int closestDist_index[6] = {0};

	//c1->aa[i].calculateCenter(false);
	//c2->aa[j].calculateCenter(false);
	//c3->aa[k].calculateCenter(false);

        //compute all 3 distances
        double dist[3];
	dist[0] = findClosestDistance(c1->aa[i],c2->aa[j],FLT_MAX,&closestDist_index[0],&closestDist_index[1]);
	dist[1] = findClosestDistance(c1->aa[i],c3->aa[k],FLT_MAX,&closestDist_index[2],&closestDist_index[3]);
	dist[2] = findClosestDistance(c2->aa[j],c3->aa[k],FLT_MAX,&closestDist_index[4],&closestDist_index[5]);

	//cout << dist[0] << " " << dist[1] << " " << dist[2] << endl;
        //cout << green << "triplet found: " << c1->aa[i].residue << " " << c2->aa[j].residue << " " << c3->aa[k].residue
         //    << " distance1 " << dist[0] << " distance2 " << dist[1] << " distance3 " << dist[2] << endl;

	//There are two cases for computing distances. If we have 1 Pi, we want the min distance from that Pi to the Anion and Cation
	//If we have 2 or more Pi's we just want the minimum 2 out of 3 distances
	if(aaCount[AATYPE_PI] == 1){
          string keystr;

          //at this point we know these pairs may or may not be in a triplet. 
          if(dist[0] < 7.0){
            keystr = buildKeyString(PDBfile,c1->aa[i],c2->aa[j]);
            //insert a '0' if this key isn't in the map. If it is already in the map, then do nothing.
            pairmapit = pairmap.find(keystr);
            if(pairmapit == pairmap.end()){
              pairmap[keystr] = 0;
            }
          }

          if(dist[1] < 7.0){
            keystr = buildKeyString(PDBfile,c1->aa[i],c3->aa[k]);
            pairmapit = pairmap.find(keystr);
            if(pairmapit == pairmap.end()){
              pairmap[keystr] = 0;
            }
          }

          if(dist[2] < 7.0){
            keystr = buildKeyString(PDBfile,c2->aa[j],c3->aa[k]);
            pairmapit = pairmap.find(keystr);
            if(pairmapit == pairmap.end()){
              pairmap[keystr] = 0;
            }
          }

          //the below code checks to see if the distance from the Pi to the anion's/cation's is above the threshold. If it is
          //we end the loop
          if(residueType[0] == AATYPE_PI){
            if(dist[0] > threshold || dist[1] > threshold)continue;
            //position 1 is a pi, we want angles from 1 to 2 and 1 to 3
            c1->aa[i].calculateAnglesPreHydrogens(c2->aa[j],closestDist_index[0],closestDist_index[1],&angle[0],&angle1[0],&angleP[0]);
            c1->aa[i].calculateAnglesPreHydrogens(c3->aa[k],closestDist_index[2],closestDist_index[3],&angle[1],&angle1[1],&angleP[1]);
          }

          if(residueType[1] == AATYPE_PI){
            if(dist[0] > threshold || dist[2] > threshold)continue;
            //position 2 is a pi, we want angles from 2 to 1 and 2 to 3

            //carefully note these indices for closestDist_index[]. the distance from 1 TO 2 was calculated, not from 2 TO 1. The
            //the indices have to be in reverse order to account for this.
            c2->aa[j].calculateAnglesPreHydrogens(c1->aa[i],closestDist_index[1],closestDist_index[0],&angle[0],&angle1[0],&angleP[0]);
            c2->aa[j].calculateAnglesPreHydrogens(c3->aa[k],closestDist_index[4],closestDist_index[5],&angle[1],&angle1[1],&angleP[1]);
          }

          if(residueType[2] == AATYPE_PI){
            if(dist[1] > threshold || dist[2] > threshold)continue;
            //position 3 is a pi, we want angles from 3 to 1 and 3 to 2
            c3->aa[k].calculateAnglesPreHydrogens(c1->aa[i],closestDist_index[3],closestDist_index[2],&angle[0],&angle1[0],&angleP[0]);
            c3->aa[k].calculateAnglesPreHydrogens(c2->aa[j],closestDist_index[5],closestDist_index[4],&angle[1],&angle1[1],&angle[1]);
          }

          //At this point if we KNOW we have a triplet. We increment the count for all pairs involved in this triplet to indicate
          //that they are involved in another triplet.
          if(dist[0] < 7.0){
            keystr = buildKeyString(PDBfile,c1->aa[i],c2->aa[j]);
            pairmap[keystr] = (pairmap[keystr] + 1);
          }

          if(dist[1] < 7.0){
            keystr = buildKeyString(PDBfile,c1->aa[i],c3->aa[k]);
            pairmap[keystr] = (pairmap[keystr] + 1);
          }

          if(dist[2] < 7.0){
            keystr = buildKeyString(PDBfile,c2->aa[j],c3->aa[k]);
            pairmap[keystr] = (pairmap[keystr] + 1);
          }
        }


	if(aaCount[AATYPE_PI] > 1){
          //if the minimum 2 distances are < threshold we want to keep this as a potential interaction.
          
          //the below checks for the inverse - which is that at least 2 distances are above the threshold.
          if(((dist[0] > threshold) && (dist[1] > threshold)) ||
             ((dist[1] > threshold) && (dist[2] > threshold)) ||
             ((dist[0] > threshold) && (dist[2] > threshold))){
            continue;
          }
        }
#ifdef DEBUG
         cout << purple << "triplet within distance found " << c1->aa[i].residue << " " << c2->aa[j].residue << " " << c3->aa[k].residue << " "
             << " distance1 " << dist[0] << " distance2 " << dist[1] << " distance3 " << dist[2] << endl;
#endif

         output_file << c1->aa[i].residue                          << ","
                      << c2->aa[j].residue                  << ","
                      << c3->aa[k].residue                  << ","
                      << dist[0]                            << ","
                      << dist[1]                            << ","
                      << dist[2]                            << ","
                      << c1->aa[i].atom[0]->resSeq          << ","
                      << c2->aa[j].atom[0]->resSeq          << ","
                      << c3->aa[k].atom[0]->resSeq          << ","
                      << PDBfile.filename                   << ","
                      << PDBfile.resolution                 << ","
                      << PDBfile.model_number               << ","
                      << c1->aa[i].atom[0]->chainID         << ","
                      << c2->aa[j].atom[0]->chainID         << ","
                      << c3->aa[k].atom[0]->chainID         << ","
                      << angle[0]                           << ","
                      << angle1[0]                          << ","
                      << angleP[0]                          << ","
                      << angle[1]                           << ","
                      << angle1[1]                          << ","
                      << angleP[1]                          << ","
                      << endl;

         if(opts.pairlistfile != NULL){
           for(pairmapit=pairmap.begin(); pairmapit!=pairmap.end(); ++pairmapit)
           {
             pairlistfile << (*pairmapit).first << "," << (*pairmapit).second << endl;
           }
         }
          
      }
    }
  }
}

string buildKeyString(PDB & PDBfile,
                      AminoAcid &aa1,
                      AminoAcid &aa2){
  ostringstream keystr;
  keystr << aa1.residue << "," << aa1.atom[0]->resSeq << "," << aa1.atom[0]->chainID << ","
         << aa2.residue << "," << aa2.atom[0]->resSeq << "," << aa2.atom[0]->chainID
         << PDBfile.filename << "," << PDBfile.model_number;
  return keystr.str();
}


void searchLigandsInformation(PDB & PDBfile,
                              Residue & ligand,
                              unsigned int chain1,
                              string residue1,
                              Options & opts,
                              ofstream& output_file)
{
  Chain* c1 = &(PDBfile.chains[chain1]);
  for(int i=0; i<c1->aa.size(); i++)
    {
      if(c1->aa[i].residue == residue1 && !(c1->aa[i].skip))
        {
          findBestInteraction(c1->aa[i],
                              ligand,
                              opts.threshold,
                              PDBfile,
                              opts.gamessfolder,
                              true,
                              output_file);
        }
    }
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
	//cout << "here" << endl;
      if(aa1.center[i].skip) continue;
      for( unsigned int j = 0; j < aa2.center.size(); j++ )
        {
          if(aa2.center[j].skip) continue;
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
                          PDB & PDBfile,
                          char* gamessfolder,
                          bool ligand,
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

  // AND WE HAVE A WINNER! 
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
        {
          code2 = 'M';
        }

      PDB pairWithHydrogen;

      // Set the residues and ligands to find
      pairWithHydrogen.setResiduesToFind(PDBfile.residue1, PDBfile.residue2);
      pairWithHydrogen.setLigandsToFind(PDBfile.ligandsToFind);
      

      // Add the hydrogens
      pairWithHydrogen.addHydrogensToPair(aa1,aa2,closestDist_index1,closestDist_index2);

      // Set the filename
      pairWithHydrogen.filename = PDBfile.filename;

      // Separate the pair into 2 variables
      pairWithHydrogen.getPair(aa1.atom[0]->resSeq,
                               aa2.atom[0]->resSeq,
                               &aa1h,
                               &aa2h,
                               ligand);

      if( aa2h.skip == true || aa1h.skip == true )
        {
          if(code2 == 'M')
            {
              // Now we are going to undo the alt loc markings
              aa1.unmarkAltLocAtoms();
              aa2.unmarkAltLocAtoms();
            }

          return;
        }

      float dist;
      float distOxy;
      float distOxy2;

      if(!aa1h.calculateDistancesAndAnglesPostHydrogens(aa2h,
                                                        aa2.center[closestDist_index2],
                                                        threshold,
                                                        &dist,
                                                        &distOxy,
                                                        &distOxy2,
                                                        &angleh,
                                                        &angleOxy,
                                                        &angleOxy2))
        {
          return;
        }

      // calculate the angles of this interaction
      aa1.calculateAnglesPreHydrogens(aa2,
                                      closestDist_index1,
                                      closestDist_index2,
                                      &angle,
                                      &angle1,
                                      &angleP);

      // The following is pretty hackish.  For the time being since we don't
      // have an agreement on how to deal with the PO4 ligands completely, 
      // we will take 1 H out at a time and output the GAMESS input file.
      // Thus, for each PO4, we will have 3 GAMESS input files.
      // Good thing there aren't too many of these
      if(aa2h.residue == "PO4")
        {
          vector<Atom*>::iterator atom_iterator;
          int count = 0;
          for( atom_iterator  = aa2h.atom.begin();
               atom_iterator != aa2h.atom.end();
               ++atom_iterator )
            {
              if((*atom_iterator)->name == " H  ")
                {
                  (*atom_iterator)->skip = true;
                  if( gamessfolder )
                    {
                      numOutputted++;
                      sprintf(output_filename, "%s/gamessinp-%d.inp", gamessfolder, numOutputted);
                      outputINPfile(PDBfile.filename, output_filename, aa1h, aa2h);
                    }
                  else
                    {
                      sprintf(output_filename, "N/A");
                    }

                  // and we finally output some results!
                  output_file << aa1.residue                        << ","
                              << aa2.residue                        << ","
                              << closestDist                        << ","
                              << angle                              << ","
                              << angleP                             << ","
                              << angle1                             << ","
                              << aa1.atom[0]->resSeq                << ","
                              << aa2.atom[0]->resSeq                << ","
                              << code1 << code2                     << ","
                              << PDBfile.filename                   << ","
                              << PDBfile.resolution                 << ","
                              << output_filename                    << ","
                              << aa1.atom[0]->chainID               << ","
                              << aa2.atom[0]->chainID               << ","
                              << aa1.center[closestDist_index1]     << ","
                              << aa2.center[closestDist_index2]     << ","
                              << aa1h.center[0]                     << ","
                              << aa2h.center[0]                     << ","
                              << dist                               << ","
                              << distOxy                            << ","
                              << distOxy2                           << ","
                              << angleh                             << ","
                              << angleOxy                           << ","
                              << angleOxy2                          << endl;
                  (*atom_iterator)->skip = false;
                  if(count == 2) break;
                  ++count;
                }
            }
        }
      else 
        {
          // Here, we are outputting the
          if( gamessfolder )
            {
              numOutputted++;
              sprintf(output_filename, "%s/gamessinp-%d.inp", gamessfolder, numOutputted);
              outputINPfile(PDBfile.filename, output_filename, aa1h, aa2h);
            }
          else
            {
              sprintf(output_filename, "N/A");
            }

          // and we finally output some results!
          output_file << aa1.residue                        << ","
                      << aa2.residue                        << ","
                      << closestDist                        << ","
                      << angle                              << ","
                      << angleP                             << ","
                      << angle1                             << ","
                      << aa1.atom[0]->resSeq                << ","
                      << aa2.atom[0]->resSeq                << ","
                      << code1 << code2                     << ","
                      << PDBfile.filename                   << ","
                      << PDBfile.resolution                 << ","
                      << PDBfile.model_number               << ","
                      << output_filename                    << ","
                      << aa1.atom[0]->chainID               << ","
                      << aa2.atom[0]->chainID               << ","
                      << aa1.center[closestDist_index1]     << ","
                      << aa2.center[closestDist_index2]     << ","
                      << aa1h.center[0]                     << ","
                      << aa2h.center[0]                     << ","
                      << dist                               << ","
                      << distOxy                            << ","
                      << distOxy2                           << ","
                      << angleh                             << ","
                      << angleOxy                           << ","
                      << angleOxy2                          << endl;
        }

    }
}

void outputINPfile(string input_filename, char* filename, AminoAcid& aa1h, AminoAcid& aa2h)
{
  ofstream inpout(filename);

  inpout << INPheader << endl;
  inpout << " $MOROKM IATM(1)=" << aa1h.atom.size() << ",";
  if(aa2h.residue == "PO4")
    {
      inpout<< aa2h.atom.size() - 1 << " ICHM(1)=0,-1" << " $END" << endl;
    }
  else
    {
      inpout<< aa2h.atom.size() << " ICHM(1)=0,-1" << " $END" << endl;
    }
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
      if( !aa2h.atom[i]->skip )
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
          else if( aa2h.atom[i]->element == " P")
            {
              inpout << "P     15.0     ";
            }
          inpout << aa2h.atom[i]->coord << endl;
        }
    }

  inpout << " $END" << endl;
  inpout.close();
}

void write_output_head(ofstream& out,bool triplets)
{

if(triplets){
	out << "#res1,res2,res3,distance1to2,distance1to3,distance2to3,loc1,loc2,loc3,pdbID,resolution,model,chain1,chain2,chain3" << endl;

}
else{
  out <<"#res1,res2,dist,angle,angleP,angle1,loc1,loc2,code,pdbID,resolution,model,gamessinput,chain1,chain2,center1,,,center2,,,center1h,,,center2h,,,dist,distOxy,distOxy2,angleh,angleOxy,angleOxy2,gamessoutput,electrostatic(Hartree),electrostatic(kcal/mol),exchangerep(Hartree),exchangerep(kcal/mole),polarization(Hartree),polarization(kcal/mole),chargexfer(Hartree),chargexfer(kcal/mol),highordercoup(Hartree),highordercoup(kcal/mole),totalinter(Hartree),totalinter(kcal/mole)" << endl;
}
}
