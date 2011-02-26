/****************************************************************************************************/
//  COPYRIGHT 2011, University of Tennessee
//  Author: David Jenkins (david.d.jenkins@gmail.com)
//  File: AminoAcid.hpp
//  Date: 16 Jan 2011
//  Version: 1.0
//  Description: Holds the class declaration for AminoAcid
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


#ifndef __AMINOACID_HPP__
#define __AMINOACID_HPP__

#include <cstdlib>
#include <cstdio>
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include "Atom.hpp"

// Defines to set indexes for the vector that holds plane information
#define CG_PLANE_COORD_PTT  0 
#define CD1_PLANE_COORD_PTT 1
#define CD2_PLANE_COORD_PTT 2

#define C__PLANE_COORD_AG  2
#define O_1_PLANE_COORD_AG 0
#define O_2_PLANE_COORD_AG 1


class AminoAcid{
private:
  // calculates the centers of the AA and sets the plane
  // information order to calculate the angle later
  void centerPHEorTYR();
  void centerTRP();
  void centerASP();
  void centerGLU();
  // however these does set the plane information correctly
  void centerASP_nocharge();
  void centerGLU_nocharge();
  void printPHEorTYR(FILE* output);
  void printASP(FILE* output);
  void printGLU(FILE* output);

public:
  // constructor
  AminoAcid();

  // destructor
  ~AminoAcid();

  // calculate the center of the AA. Calls the individual
  // functions above depending on AA
  void calculateCenter(bool center);

  void printNeededAtoms(FILE* output);

  // this holds pointers to the ATOM strings
  vector<Atom*> atom;
  // holds the coordinates of the atoms to calculate the plane
  vector< vector <Atom*> > plane_info;
  // holds the calculated centers that will be examined
  vector<Coordinates>  center;
  // holds the residue name
  string residue;
  // true of this AA isn't important, false otherwise
  bool skip;
  // true if there is an alt loc in ATOM line
  bool altLoc;

  string line;

  friend ostream& operator<<(ostream& output, const AminoAcid& p);

};

#endif
