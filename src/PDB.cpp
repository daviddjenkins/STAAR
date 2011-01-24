/****************************************************************************************************/
//  COPYRIGHT 2011, University of Tennessee
//  Author: David Jenkins (david.d.jenkins@gmail.com)
//  File: PDB.cpp
//  Date: 12 Jan 2011
//  Version: 1.0
//  Description: Actual does all PDB file reading. Also splits info into chains for further
//               analysis.
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
#include <algorithm>
#include "PDB.hpp"
#include "Utils.hpp"

// Constructor to initialize the PDB class object by 
// ensuring all the vectors are empty
PDB::PDB()
{
  filename = NULL;
  chains.clear();
  atoms.clear();
  hetatms.clear();
  seqres.clear();
}

// Constructor to parse the inputted PDB file
PDB::PDB(char* fn, Options& opt)
{
  failure = false;
  parsePDB(fn, opt);
}

// Destructor to empty the arrays
PDB::~PDB()
{
  filename = NULL;
  chains.clear();
  atoms.clear();
  hetatms.clear();
  seqres.clear();
}

bool PDB::fail()
{
  return failure;
}

// Parses the given PDB file
void PDB::parsePDB(char * fn, Options& opt)
{
  filename = fn;

  // Open the file
  ifstream PDBfile(fn, ifstream::in);
  string line; // this is a temp var to hold the current line from the file
  
  // Ensure the file opened correctly
  if( !PDBfile.is_open() || PDBfile.fail() ){
    cerr << "Failed to open PDB file " << fn << endl;
    perror("\t");
    failure=true;
    return;
  }

  // For each line in the file
  while( getline(PDBfile, line) && !failure)
    {
      // Parse the line if we are on a SEQRES line
      int found = line.find("SEQRES");
      if( found == 0 )
        {
          Seqres s(line);
          failure = s.fail();
          seqres.push_back(s);
        }

      // Parse the line if we are on an ATOM line
      found = line.find("ATOM");
      if( found == 0 )
        {
          Atom a(line);
          failure = a.fail();
          atoms.push_back(a);
        }

      // This is current commented out because we aren't
      // doing anything with the HETATM lines so there
      // is no reason to read them
      // // Parse the line if we are on a HETATM line
      // found = line.find("HETATM");
      // if( found == 0 )
      //   {
      //     Atom h(line);
      //     //failure = h.fail();
      //     hetatms.push_back(h);
      //   }
    }

  
  PDBfile.close();
  
}

// Search the chain by id
//    - Not used, right now
vector<Chain>::iterator PDB::findChainNumber(char id)
{
  Chain c(id);
  return find(chains.begin(), chains.end(), c);
}

// Organizes the the data by chains
void PDB::populateChains(Options& opt)
{
  // Flag to hold the current chain id
  char chainID =  '-';
  // Index for the chain
  int chainIndex = -1;
  
  // Go through each atom
  for(unsigned int i = 0; i < atoms.size(); i++)
    {
      // Check if this is a new chain
      if(atoms[i].chainID != chainID )
        {
          chainIndex++;
          chainID = atoms[i].chainID;
          chains.resize(chainIndex+1);
          chains[chainIndex].id = chainID;
        }

      AminoAcid aa;
      unsigned int residue_number = atoms[i].resSeq;
      while(residue_number == atoms[i].resSeq)
	{
	  aa.atom.push_back(&atoms[i]);
	  i++;
	  if( i == atoms.size() )
	    {
	      break;
	    }
	}
      i--;
      
      aa.residue = atoms[i].residueName;
      if(aa.residue == "TRP" || aa.residue == "PHE" ||
	 aa.residue == "TYR" || aa.residue == "ASP" ||
	 aa.residue == "GLU")
	{
	  aa.calculateCenter(opt.center);
	}
      else
	{
	  aa.skip = true;
	}


      // Store the reference of the atom in the corresponding chain info
      chains[chainIndex].addAminoAcid(aa);
    }
  
  //reset the flags
  chainID = '-';
  chainIndex = -1;

  // This is current commented out because we aren't
  // doing anything with the HETATM lines so there
  // is no reason to read them
  // // Go through each hetatm
  // for(int i = 0; i < hetatms.size(); i++)
  //   {
  //     // Check if this is a new chain
  //     if(hetatms[i].chainID != chainID )
  //       {
  //         chainIndex++;
  //         chainID = hetatms[i].chainID;
  //       }
  //     // Store the reference of the hetatm in the corresponding chain info
  //     chains[chainIndex].addHetatm(&hetatms[i]);
  //   }

  // // again, reset the flags
  // chainID = '-';
  // chainIndex = -1;

  // Go through each seqres
  for(unsigned int i = 0; i < seqres.size(); i++)
    {
      // Check if this is a new chain
      if(seqres[i].chainID != chainID )
        {
          chainIndex++;
          chainID = seqres[i].chainID;
        }

      //Store the reference of the seqres in the corresponding chain info
      chains[chainIndex].addSeqres(&seqres[i]);
    }

}
