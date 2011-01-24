/****************************************************************************************************/
//  COPYRIGHT 2011, University of Tennessee
//  Author: David Jenkins (david.d.jenkins@gmail.com)
//  File: Seqres.hpp
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


#ifndef __SEQRES_HPP__
#define __SEQRES_HPP__

#include <vector>
#include "Utils.hpp"

class Seqres
{
private:
  // True when parsing failed, false otherwise
  bool failure;

public:
  // Constructor to reset everything
  Seqres();
  // Constructor to parse the SEQRES line
  Seqres(string line);
  // Destructor to reset everything
  ~Seqres();
  // Parses the inputted SEQRES line
  int parseLine(string line);
  // Returns true when parsing failed, false otherwise
  bool fail();

  unsigned int serialNumber;     // Serial number     : 8-10
  char chainID;                  // Chain ID          : 12
  unsigned int numberOfResidues; // Number of residues: 14-17
};


#endif
