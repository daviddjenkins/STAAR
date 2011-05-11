/****************************************************************************************************/
//  COPYRIGHT 2011, University of Tennessee
//  Author: David Jenkins (david.d.jenkins@gmail.com)
//  File: Options.hpp
//  Date: 10 Jan 2011
//  Version: 1.0
//  Description: Header file for the command line parser and Options class
//
//  Updates: Added the ability to pick same chain or not (11 Feb 2011)
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

#ifndef __OPTIONS_HPP__
#define __OPTIONS_HPP__

#include <cstdlib>
#include <cstdio>
#include <iostream>
#include <cstring>
#include <string>
#include <unistd.h>
#include <getopt.h>
#include "Utils.hpp"


using namespace std;

class Options{
private:
  // True if parsing failed, false otherwise
  bool failure;
public:
  bool   center;                // Flag indicating to perform center calculations
  char  *pdbfile;               // Directory that the PDB files are stored
  char  *outputfile;            // File that the results are written to
  float  threshold;             // Holds the cutoff for what is considered a close
                                // distance between residues
  bool sameChain;               // Flag indicating whether to look only on same 
                                // chain for interactions or not
  vector<string>residue1;       // first group of residues to be matches with...
  vector<string>residue2;       // ...these!
  vector<string>ligands;        // Ligands that residue1 will be matches with
  int numLigands;               // Number of ligands
  char* gamessfolder;           // Directory in which all of the GAMESS INP files
                                // will be stored.
  bool outputGamessINP;         // True if GAMESS INP are to be outputted
  char* pdblist;                // List of PDB files to parse
  string extension;             // extension to use to append to pdb list files
  float resolution;             // Max resolution cut-off
  char* chain_list;             // like pdblist, but contains chains to search in

  // Constructor that sets everything to empty stuff
  Options();  
  // Constructor to parse the command line args
  Options(int argc, char **argv);
  // Empty Destructor
  ~Options();

  // Parse the command line options
  void parseCmdline(int argc, char **argv);
  // Return true of cmd line parsing failed, false otherwise
  bool fail();

};

#endif
