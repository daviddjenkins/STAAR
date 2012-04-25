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
#ifndef NO_BABEL
#include <openbabel/obconversion.h>
#include <openbabel/mol.h>
#endif
#include "AminoAcid.hpp"
#include "Atom.hpp"
#include "Seqres.hpp"
#include "Utils.hpp"
#include "Chain.hpp"


static char INPheader[] = \
  " $CONTRL SCFTYP=RHF RUNTYP=EDA ICHARG=-1 MULT=1 COORD=CART MAXIT=200 $END\n" \
  " $SYSTEM TIMLIM=1000 $END\n"                                         \
  " $BASIS GBASIS=TZV $END\n"                                           \
  " $GUESS GUESS=HUCKEL $END\n"                                         \
  " $SCF SOSCF=.F. DAMP=.T. SHIFT=.T. DEM=.F. $END\n";


  // " $CONTRL SCFTYP=RHF  RUNTYP=EDA ICHARG=-1 MULT=1 COORD=CART MAXIT=200 $END\n" \
  // " $SYSTEM TIMLIM=1000 $END\n"                                         \
  // " $BASIS GBASIS=TZV $END\n"                                           \
  // " $GUESS GUESS=HUCKEL $END\n"                                         \
  // " $SCF SOSCF=.F. DAMP=.T. SHIFT=.T. DEM=.F. $END";

#define PH_LEVEL 7.4

// Error conditions
#define FAILED_TO_OPEN_FILE         -1
#define RESOLUTION_NOT_APPLICABLE   -2
#define RESOLUTION_BLANKS           -3
#define RESOLUTION_TOO_HIGH         -4
#define RESOLUTION_TO_NUMBER_FAILED -5
#define RESOLUTION_BLANK            -6
#define NO_RESOLUTION               -7
#define MODEL_TO_NUMBER_FAILED      -8
#define MULTIPLE_MODELS_SKIP        -9
class PDB
{

private:
  // Finds the chain id in the chains vector
  // Returns an iterator to the chain
  vector<Chain>::iterator findChainNumber(char id);  

  // Indicates parsing success or failure
  bool failure;

  void parsePDBstream(istream& PDBfile, float resolution);
  bool atomsCompare();
public:
  // Default constructor that ensures everything is empty
  PDB();

  // Constructor that parses the file pointed to by 
  // supplied filename
  PDB(const char* fn, float resolution);

  // Constructor that parses the supplied file
  PDB(istream& file, float resolution);

  // Destructed that empties everything
  ~PDB();

  // Deletes everything 
  void clear();

  // Returns true if the parsing failed
  bool fail();

  // Print a failure message
  void printFailure();
  
  // Parses the file and stores the data
  void parsePDB(const char* fn, float resolution);
  void parsePDB(istream& file, float resolution);

#ifndef NO_BABEL
  // Calls Babel to add the hydrogens and inputs them into the PDB
  void addHydrogensToPair(AminoAcid& a, AminoAcid& b, int cd1, int cd2);
#endif

  // Organizes the data read from parsePDB into chains
  void populateChains(bool center);

  // Organizes ligands into an array
  void findLigands(vector<string> ligandsToFind);
  
  void getPair(int& resSeq1, 
               int& resSeq2, 
               Residue* r1, 
               Residue* r2, 
               bool ligand);

  void setResiduesToFind(vector<string>* r1,
                         vector<string>* r2);

  void setLigandsToFind(vector<string>* l);
  // Puts the atoms in order by their sequence number
  void sortAtoms();

  vector<string>* ligandsToFind;
  vector<string>* residue1;
  vector<string>* residue2;

  vector<Chain>           chains;         // Variable to hold the chain information
  vector<Atom>            atoms;          // Vector hold all the atom lines
  vector<Atom>            hetatms;        // Vector holding all the hetatm lines
  vector<Residue*>        ligands;        // Vector holding all the ligand lines
  vector<Seqres>          seqres;         // Vector holding all the seqres lines
  vector<string>          conect;         // Vector holding all the CONECT lines

  const char*             filename;       // Holds the filename, if needed
#ifndef NO_BABEL
  OpenBabel::OBConversion conv;           // Holds OpenBabel reading of important
                                          //  residue pairs adding H. Will also 
                                          //  be used to output to GAMESS format
#endif
  int failflag;
  float resolution;

  int model_number;
  vector<PDB> models;                     // Vector of models in the PDB
  // NOTE: the way I handle models is not intelligent!  I did a quick and dirty 
  // implementation just to get it done.  What makes it stupid is that I just
  // creates duplicates of stuff in memory because I didn't want to spend
  // must more time on this to make it pretty.  Someone else can do that if they
  // want, but they will have to change a bunch of code.  Had I known about how
  // the faculty wanted to handle multiple models in the beginning, it would have
  // been trivial to code it up.  But now that there is just so much code centered
  // on a single model, it's a lot of work to add in now.  What should be done is 
  // create a model class that will essentially take the place of what the PDB class
  // is right now.  Then the PDB class will just contain a vector of all the models
  // But having it the way it is now is not much of a problem anyway because the 
  // PDBs aren't that large anyway so it isn't going to completely dominate the 
  // system's memory...hopefully.

  friend ostream& operator<<(ostream& output, const PDB& p);
};
  


#endif
