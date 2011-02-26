/****************************************************************************************************/
//  COPYRIGHT 2011, University of Tennessee
//  Author: David Jenkins (david.d.jenkins@gmail.com)
//  File: PDB.hpp
//  Date: 12 Jan 2011
//  Version: 1.0
//  Description: Contains the class definition and implementation for PDB
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


#ifndef __PDB_HPP__
#define __PDB_HPP__

#include <vector>
#include <algorithm>
#include <sys/stat.h>
#include <openbabel/obconversion.h>
#include <openbabel/mol.h>
#include "Options.hpp"
#include "AminoAcid.hpp"
#include "Atom.hpp"
#include "Seqres.hpp"
#include "Utils.hpp"
#include "Chain.hpp"

using namespace OpenBabel;

class PDB
{

private:
  // Finds the chain id in the chains vector
  // Returns an iterator to the chain
  vector<Chain>::iterator findChainNumber(char id);  

  // Indicates parsing success or failure
  bool failure;

  void parsePDBstream(istream& PDBfile);
  bool atomsCompare();
public:
  // Default constructor that ensures everything is empty
  PDB();

  // Constructor that parses the file pointed to by 
  // supplied filename
  PDB(char* fn);

  // Constructor that parses the supplied file
  PDB(istream& file);

  // Destructed that empties everything
  ~PDB();

  // Returns true if the parsing failed
  bool fail();
  
  // Parses the file and stores the data
  void parsePDB(char* fn);
  void parsePDB(istream& file);

  // Calls Babel to add the hydrogens and inputs them into the PDB
  void addHydrogensToPair(AminoAcid& a, AminoAcid& b);

  // Organizes the data read from parsePDB into chains
  void populateChains(bool center);

  void sortAtoms();

  vector<Chain>         chains;         // Variable to hold the chain information
  vector<Atom>          atoms;          // Vector hold all the atom lines
  vector<Atom>          hetatms;        // Vector holding all the hetatm lines
  vector<Seqres>        seqres;         // Vector holding all the seqres lines
  char*                 filename;       // Holds the filename, if needed
  OBConversion          conv;           // Holds OpenBabel reading of important
                                        //  residue pairs adding H. Will also 
                                        //  be used to output to GAMESS format
};
  


#endif
