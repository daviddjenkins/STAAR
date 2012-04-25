/****************************************************************************************************/
//  COPYRIGHT 2011, University of Tennessee
//  Author: David Jenkins (david.d.jenkins@gmail.com)
//  File: PDB.cpp
//  Date: 12 Jan 2011
//  Date Modified: 4 Feb 2011
//  Version: 1.0
//  Description: Actually does all PDB file reading. Also splits info into chains for further
//               analysis.
//
//  Updates: Added ability to read .gz files (4 Feb 2011)
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
#include <string>
#include <iterator>
#include "PDB.hpp"
#include "Utils.hpp"
#include "../gzstream/gzstream.h"
#include "CoutColors.hpp"

#ifndef NO_BABEL
using namespace OpenBabel;
#endif

// Constructor to initialize the PDB class object by 
// ensuring all the vectors are empty
PDB::PDB()
{
  filename      = NULL;
  chains.clear();
  atoms.clear();
  hetatms.clear();
  seqres.clear();
  ligands.clear();
  conect.clear();
  models.clear();
  ligandsToFind = NULL;
  residue1      = NULL;
  residue2      = NULL;
  resolution = -2;
  model_number=1;
}

// Constructor to parse the inputted PDB file
PDB::PDB(const char* fn, float res)
{
  failure       = false;
  ligandsToFind = NULL;
  residue1      = NULL;
  residue2      = NULL;
  resolution = -2;
  model_number=1;
  parsePDB(fn, res);
}

// Constructor to parse the inputted PDB file
PDB::PDB(istream& file, float res)
{
  failure       = false;
  ligandsToFind = NULL;
  residue1      = NULL;
  residue2      = NULL;
  resolution = -2;
  model_number=1;
  parsePDBstream(file, res);
}

// Destructor to empty the arrays
PDB::~PDB()
{
  filename      = NULL;
  ligandsToFind = NULL;
  residue1      = NULL;
  residue2      = NULL;
  chains.clear();
  atoms.clear();
  hetatms.clear();
  seqres.clear();
  ligands.clear();
  conect.clear();
  models.clear();
  resolution = -2;
  model_number=1;
}

void PDB::clear()
{
  filename      = NULL;
  ligandsToFind = NULL;
  residue1      = NULL;
  residue2      = NULL;
  chains.clear();
  atoms.clear();
  hetatms.clear();
  seqres.clear();
  ligands.clear();
  conect.clear();
  models.clear();
  resolution = -2;
  model_number=1;
}

bool PDB::fail()
{
  return failure;
}

void PDB::printFailure()
{
  if(failflag == FAILED_TO_OPEN_FILE)
    {
      cerr << red << "Error" << reset << ": Failed to open PDB file " << filename << endl;
      perror("\t");
    }
  else if(failflag == RESOLUTION_NOT_APPLICABLE || failflag == RESOLUTION_BLANK ||
          failflag == NO_RESOLUTION)
    {
      cout << cyan << "Skipping" << reset << " because no resolution was specified" << endl;
    }
  else if(failflag == RESOLUTION_TOO_HIGH)
    {
      cout << cyan << "Skipping" << reset << " because resolution was too high" << endl;
    }
  else if(failflag == RESOLUTION_TO_NUMBER_FAILED)
    {
      cout << cyan << "Skipping" << reset << " because resolution could not be converted to a number" << endl;
    }
  else if(failflag == MODEL_TO_NUMBER_FAILED)
    {
      cout << cyan << "Skipping" << reset << ": couldn't convert model number to an int" << endl;
    }
  else if(failflag == MULTIPLE_MODELS_SKIP)
    {
      cout << cyan << "Skipping" << reset << ": multiple models exist in PDB file" << endl;
    }
  else
    {
      cerr << red << "Unspecified failure" << endl;
    }
}

// Parses the given PDB file
void PDB::parsePDB(const char * fn, float resolution)
{
  filename = fn;

  // Open the file
  igzstream PDBfile(fn); 

  // Ensure the file opened correctly
  if( PDBfile.fail() )
    {
      failure = true;
      failflag = FAILED_TO_OPEN_FILE;
      return;
    }

  parsePDBstream(PDBfile, resolution);

  PDBfile.close();
}

void PDB::parsePDB(istream& file, float resolution)
{
  parsePDBstream(file, resolution);
}

void PDB::parsePDBstream(istream& PDBfile, float resolution)
{
  string line; // this is a temp var to hold the current line from the file
  int count = 0;
  bool rflag = false;
  PDB model;
  int modelnum=1;

  // For each line in the file
  while( getline(PDBfile, line) && !failure)
    {
      count++;

      // Check to see if there are more than one model
      // If there is, we will just skip the file altogether because
      // it is less of a hassle to deal with.
      // So here, we check the model number, if it is greater
      // than 1, we will just break from this loop and print
      // a note
      int found = line.find("MODEL");
      if( found == 0 )
        {
          int model_number;
          if(line.substr(10,4) == "    ")
            {
              continue;
            }
          else if( !from_string<int>(model_number, line.substr(10,4), dec) )
            {
              failure = true;
              failflag = MODEL_TO_NUMBER_FAILED;
              continue;
            }
          else if( model_number > 1 )
            {
              model.model_number = modelnum;
              models.push_back(model);
              model.clear();
              modelnum++;
            }
        }

      // Check the resolution of the PDB
      found = line.find("REMARK   2 RESOLUTION.");
      if( found == 0 )
        {
          this->resolution = 0;
          if(line.substr(23, 7) == "NOT APP")
            {
              //failure = true;
              //failflag = RESOLUTION_NOT_APPLICABLE;
              this->resolution = -1;
            }
          else if(line.substr(23, 7) == "       ")
            {
              failure = true;
              failflag = RESOLUTION_BLANK;
            }
          else
            {
              if(!from_string<float>(this->resolution, line.substr(23,7), dec))
                {
                  failure = true;
                  failflag = RESOLUTION_TO_NUMBER_FAILED;
                }
              if(this->resolution > resolution)
                {
                  failure = true;
                  failflag  = RESOLUTION_TOO_HIGH;
                }
            }
          continue;
        }

      // Parse the line if we are on an ATOM line
      found = line.find("ATOM");
      if( found == 0 )
        {
          Atom a(line, count);
          failure = a.fail();
          atoms.push_back(a);
          model.atoms.push_back(a);
          continue;
        }

      // Parse the line if we are on a HETATM line
      // Used to find ligands
      found = line.find("HETATM");
      if( found == 0 )
        {
          Atom h(line, count);
          failure = h.fail();
          hetatms.push_back(h);
          model.hetatms.push_back(h);
          continue;
        }

      // Parse the line if we are on a CONECT line
      // Used to find ligands
      found = line.find("CONECT");
      if( found == 0 )
        {
          conect.push_back(line);
          model.conect.push_back(line);
          continue;
        }
    }

  // This checks to see if we even had a resolution line
  // in the PDB.  If we don't, we are skipping it
  if( this->resolution == -2 )
    {
      failure = true;
      failflag = NO_RESOLUTION;
    }

  // This ensures that we have the last model
  // we parsed in out model list
  model.model_number = modelnum;
  models.push_back(model);
}

#ifndef NO_BABEL
// This function will call the Babel library to add 
// hydrogens to the residues
void PDB::addHydrogensToPair(AminoAcid& a, AminoAcid& b, int cd1, int cd2)
{
  OBMol mol;
  string addedH;
  istringstream tempss;
  bool ligand;
  if(b.atom[0]->line.find("HETATM") != string::npos)
    {
      ligand = true;
    }
  else
    {
      ligand = false;
    }

  // This section is just to suppress all of the 
  // warning message that aren't important to us
  {
    OBConversion apiConv;
    OBFormat* pAPI = OBConversion::FindFormat("obapi");
    if(pAPI)
      {
        apiConv.SetOutFormat(pAPI);
        apiConv.AddOption("errorlevel", OBConversion::GENOPTIONS, "0");
        apiConv.Write(NULL, &std::cout);
      }
  }

  // Now, let's pack up the information into a string
  string packedFile="";
  for(unsigned int i=0; i < a.altlocs[cd1].size(); i++)
    {
      if( !a.altlocs[cd1][i]->skip )
        {
          packedFile += a.altlocs[cd1][i]->line + "\n";
        }
    }
  
  int cd2_al = cd2;
  if(b.residue == "ASP" || b.residue == "GLU")
    {
      cd2_al = cd2%(b.altlocs.size());
    }

  for(unsigned int i=0; i < b.altlocs[cd2_al].size(); i++)
    {
      if( !b.altlocs[cd2_al][i]->skip )
        {
            packedFile += b.altlocs[cd2_al][i]->line + "\n";
        }
    }
  packedFile += a.makeConect(cd1);
  packedFile += b.makeConect(cd2_al);


  // Now, let's set up some Babel information
  // First, we get the PDB format to tell
  // Babel how to read the information and 
  // how to output it
  OBFormat* pdbformat = this->conv.FindFormat("pdb");
  this->conv.SetInFormat(pdbformat);
  this->conv.SetOutFormat(pdbformat);

  // Here is where Babel reads everything
  // and adds hydrogens to the pair
  // TO ADD: option to set pH
  this->conv.ReadString(&mol,packedFile);
  mol.AddHydrogens(false,true,PH_LEVEL);

  // Let's write the newly written hydrogens to 
  // a string and parse it
  addedH = this->conv.WriteString(&mol);
  tempss.str(addedH);

  // This ensures that the ligand hydrogens are labeled as
  // HETATM instead of ATOM just for the sake of STAAR. 
  // This may be wrong, but it should be fine since we are
  // stripping out that information later when we write 
  // the GAMESS inp files
  if( ligand )
    {
      string line;
      string f = "";
      while( getline(tempss,line) )
        {
          if( line.find(b.residue) != string::npos )
            {
              line.replace(0,6,"HETATM");
            }
          f += line + "\n";
        }
      tempss.seekg(ios_base::beg);
      tempss.clear();
      tempss.str(f);
    }

  this->failure = false;
  this->parsePDB(tempss,99999.99);

  // This is just to ensure that all of the atoms
  // are grouped together because Babel just 
  // appends the H to the end of the file
  if( !ligand )
      this->sortAtoms();

  // Split the atoms up into amino acids and chains
  this->populateChains(true);
}
#endif

// Search the chain by id
//    - Not used, right now
vector<Chain>::iterator PDB::findChainNumber(char id)
{
  Chain c(id);
  return find(chains.begin(), chains.end(), c);
}

void PDB::setResiduesToFind(vector<string>* r1,
                            vector<string>* r2)
{
  residue1 = r1;
  residue2 = r2;
}

void PDB::setLigandsToFind(vector<string>* l)
{
  ligandsToFind = l;
}

// Organizes the the data by chains
void PDB::populateChains(bool center)
{
  // Flag to hold the current chain id
  char chainID =  '-';
  // Index for the chain
  int chainIndex = -1;
  vector<char>chainIDs;
  
  // Go through each atom
  for(unsigned int i = 0; i < atoms.size(); i++)
    {
      if(atoms[i].chainID != chainID )
        {
          chainID = atoms[i].chainID;
          chainIndex = distance(chainIDs.begin(),find(chainIDs.begin(),chainIDs.end(),chainID));
          if( chainIndex == chainIDs.size() )
            {
              chainIDs.push_back(chainID);
              chains.resize(chainIndex+1); 
              chains[chainIndex].id = chainID;
            }
        }

      // Separate the atoms into amino acids
      AminoAcid aa;
      unsigned int residue_number = atoms[i].resSeq;
      char iCode = atoms[i].iCode;
      vector<char> altloc_ids;
      while(residue_number == atoms[i].resSeq &&
            atoms[i].chainID == chainID &&
            atoms[i].iCode == iCode )
        {
          // Determine how many alternate locations there are
          if(atoms[i].altLoc != ' ')
            {
              vector<char>::iterator found = find(altloc_ids.begin(), altloc_ids.end(), atoms[i].altLoc);
              if(found == altloc_ids.end())
                {
                  altloc_ids.push_back(atoms[i].altLoc);
                }
            }
          aa.atom.push_back(&atoms[i]);
          i++;
          if( i == atoms.size() )
            {
              break;
            }
        }
      i--;
      aa.determineAltLoc(altloc_ids);
      aa.residue = atoms[i].residueName;
      vector<string>::iterator found1 = find(residue1->begin(), residue1->end(), aa.residue);
      vector<string>::iterator found2 = find(residue2->begin(), residue2->end(), aa.residue);
      if( found1 != residue1->end() || found2 != residue2->end())
        {
          aa.calculateCenter(center);
        }
      else
        {
          aa.skip = true;
        }

      // This is a little hack to skip over residue locations that have 
      // insertion codes.  Mainly because this adds more complexity than
      // we need.
      if( aa.atom[0]->iCode != ' ' )
        {
          aa.skip = true;
        }

      // Store the reference of the atom in the corresponding chain info
      chains[chainIndex].addAminoAcid(aa);
    }
  
  //reset the flags
  chainID = '-';
  chainIndex = -1;

  // Go through each hetatm
  if( ligandsToFind )
    {
      for(int i = 0; i < hetatms.size(); i++)
        {
          // Check if this is a new chain
          if(hetatms[i].chainID != chainID )
            {
              chainID = hetatms[i].chainID;
              chainIndex = distance(chainIDs.begin(),find(chainIDs.begin(),chainIDs.end(),chainID));
              if( chainIndex == chainIDs.size() )
                {
                  chainIDs.push_back(chainID);
                  chains.resize(chainIndex+1); 
                  chains[chainIndex].id = chainID;
                }
            }
      
          // Separate the hetatms into Residues
          Residue r;
          unsigned int residue_number = hetatms[i].resSeq;
          while(residue_number == hetatms[i].resSeq &&
                hetatms[i].chainID == chainID)
            {
              r.atom.push_back(&hetatms[i]);
              i++;
              if( i == hetatms.size() )
                {
                  break;
                }
            }
          i--;
          r.residue = hetatms[i].residueName;
          vector<string>::iterator found1 = find(ligandsToFind->begin(), ligandsToFind->end(), r.residue);
          if(found1 != ligandsToFind->end())
            {
              r.calculateCenter(center);
            }
      
          // Store the reference of the hetatm in the corresponding chain info
          chains[chainIndex].addHetatm(r);
        }
    }
}

// this function just puts the benzene in r1 and formate or 
// ligand in r2
void PDB::getPair(int& resSeq1, 
                  int& resSeq2, 
                  Residue* r1, 
                  Residue* r2,
                  bool ligand)
{
  // The following is just to put the benzene in aa1h
  // and the formate in aa2h just so that I can keep
  // this straight in my head.  This may not be necessary
  // but it helps me when I am looking through the code.
  if( !ligand )
    {
      // If these residues were in different chains
      if( this->chains.size() != 1 )
        {
          // If the benzene was naturally first,
          // it is in the first chain while the formate
          // is in the second chain
          if(this->chains[0].aa[0].residue == "PHE")
            {
              *r1 = this->chains[0].aa[0];
              *r2 = this->chains[1].aa[0];
            }
          // otherwise they are in the opposite order
          else
            {
              *r1 = this->chains[1].aa[0];
              *r2 = this->chains[0].aa[0];
            }
        }
      // If these residues were in the same chain
      else
        {
          // If the benzene was naturally first,
          // it is in the first chain while the formate
          // is in the second chain
          if(this->chains[0].aa[0].residue == "PHE")
            {
              *r1 = this->chains[0].aa[0];
              *r2 = this->chains[0].aa[1];
            }
          // otherwise they are in the opposite order
          else
            {
              //cout << *this << endl;
              *r1 = this->chains[0].aa[1];
              *r2 = this->chains[0].aa[0];
            }
        }
      if( (*r2).removeExcessHydrogens(this->conect) )
        {
          cout << brown << "Corrected" << reset << ": " << filename
               << " | " << (*r1).residue << (*r1).atom[0]->resSeq << " Chain " << (*r1).atom[0]->chainID
               << " - " << (*r2).residue << (*r2).atom[0]->resSeq << " Chain " << (*r2).atom[0]->chainID << endl;
        }
    }
  else
    {
      this->findLigands(*ligandsToFind);
      *r1 = this->chains[0].aa[0];
      *r2 = *(this->ligands[0]);
    }
}

void PDB::findLigands(vector<string> ligandsToFind)
{
  int ligsize = ligandsToFind.size();

  for(int i=0; i<chains.size(); i++)
    {
      int hetsize = this->chains[i].hetatms.size();
      for(int j=0; j<hetsize; j++)
        {
          for(int k=0; k<ligsize; k++)
            {
              if(this->chains[i].hetatms[j].residue == ligandsToFind[k])
                {
                  this->ligands.push_back(&(this->chains[i].hetatms[j]));
                }
            }
        }
    }
}

void PDB::sortAtoms()
{
  sort(atoms.begin(),atoms.end());
}


ostream& operator<<(ostream& output, const PDB& p) 
{
  for(unsigned int i = 0; i<p.atoms.size(); i++)
    {
      output << p.atoms[i] << endl;
    }
  return output;
}
