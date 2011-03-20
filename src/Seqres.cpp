/****************************************************************************************************/
//  COPYRIGHT 2011, University of Tennessee
//  Author: David Jenkins (david.d.jenkins@gmail.com)
//  File: Seqres.cpp
//  Date: 12 Jan 2011
//  Version: 1.0
//  Description: Essentially just parses the SEQRES lines
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
#include <string>
#include "Seqres.hpp"
#include "CoutColors.hpp"

using namespace std;

// Default constructor to reset everything
Seqres::Seqres()
{
  serialNumber = 0;
  chainID = '-';
  numberOfResidues = 0;
  failure = false;
}

// Constructor used to call the parsing function
Seqres::Seqres(string line)
{
  failure = false;
  parseLine(line);
}

// Destructor to reset everything
Seqres::~Seqres()
{
  serialNumber = 0;
  chainID = '-';
  numberOfResidues = 0;
  failure = false;
  //residues.clear();
}

// Returns true when parsing failed, false otherwise
bool Seqres::fail()
{  
  return failure;
}

// Parses the SEQRES line
int Seqres::parseLine(string line)
{
  if(line[line.length()-1] == '\r')
    {
      line.erase(line.length()-1,1);
    }

  // Error check to ensure the file is formatted correctly
  if(line.length() != 80)
    {
      cerr << red << "Error" << reset << ":Possible malformed PDB file on SEQRES line." << endl;
      failure = true;
    }

  // Grab the serial number, chain id, and number of residues
  if(!from_string<unsigned int>(serialNumber,line.substr(7,3),dec))
    {
      cerr << red << "Error" << reset << ":failed to convert serial number into an unsigned int in Seqres::parseLine" << endl;
      failure = true;
      return -1;
    }

  chainID = line[11];

  if(!from_string<unsigned int>(numberOfResidues,line.substr(13,4),dec))
    {
      cerr << red << "Error" << reset << ":failed to convert serial number into an unsigned int in Seqres::parseLine" << endl;
      failure = true;
      return -1;
    }
  return 0;
}
