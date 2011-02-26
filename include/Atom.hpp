/****************************************************************************************************/
//  COPYRIGHT 2011, University of Tennessee
//  Author: David Jenkins (david.d.jenkins@gmail.com)
//  File: Atom.hpp
//  Date: 12 Jan 2011
//  Version: 1.0
//  Description: Contains the class definition and implementation for Atom
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


#ifndef __ATOM_HPP__
#define __ATOM_HPP__

#include "Coordinates.hpp"
#include "Utils.hpp"

#define MASS_C 12.0107
#define MASS_N 14.00674
#define MASS_O 15.9994
#define CHARGE_C 0.898571
#define CHARGE_O -0.83462
#define CHARGE_H -0.229336

class Atom
{
private:

public:
  // Default constructor to reset everything
  Atom();
  // Constructor that will parse an ATOM or HETATM line
  Atom(string line);
  // Constructor that will parse an ATOM or HETATM line
  Atom(char* line);
  // Destructor to reset everything
  ~Atom();
  // Parses an ATOM or HETATM line
  void parseAtom(string line);
  // Prints out all the values space delimited
  void print();
  // Prints out atom in pdb format to the given FILE*
  void print(FILE* output);
  // Returns true if parsing failed, false otherwise
  bool fail();

  unsigned int  serialNumber;  //Atom serial number: 7-11
  string        name;          //Atom name:          13-16
  char          altLoc;        //Alt. location:      17
  string        residueName;   //Residue name:       18-20
  char          chainID;       //Chain identifier:   22
  unsigned int  resSeq;        //Residue sequence #: 23-26
  char          iCode;         //Insertion code:     27
  Coordinates   coord;         //Coordinates:        31-54
  double        occupancy;     //Occupancy:          55-60
  double        tempFactor;    //Temperature Factor: 61-66
  string        element;       //Element Symbol:     77-78
  string        charge;        //Atom charge:        79-80
  bool failure;

  string line;

  friend ostream& operator<<(ostream& output, const Atom& p);
  bool operator<(const Atom& rhs) const;

};


#endif
