/****************************************************************************************************/
//  COPYRIGHT 2011, University of Tennessee
//  Author: David Jenkins (david.d.jenkins@gmail.com)
//  File: Chain.cpp
//  Date: 12 Jan 2011
//  Version: 1.0
//  Description: Implements some functions for the Chain class.  Functions include ones that
//               add Atom* and Seqres* values to different vectors.
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
#include "Chain.hpp"

// Constructor to reset everything
Chain::Chain()
{
  id = '-';
  aa.clear();
  hetatms.clear();
  seqres.clear();
}

// Constructor that sets the chain id
Chain::Chain(char i)
{
  id = i;
}

// Destructor to reset everything
Chain::~Chain()
{
  aa.clear();
  hetatms.clear();
  seqres.clear();
  id = '-';
}

// Adds reference to an ATOM to a vector of Atom*
void Chain::addAminoAcid(AminoAcid a)
{
  aa.push_back(a);
}

// Adds reference to a HETATM to a vector of Atom*
void Chain::addHetatm(Atom* h)
{
  hetatms.push_back(h);
}

// Adds reference to a SEQRES to a vector of Seqres*
void Chain::addSeqres(Seqres* s)
{
  seqres.push_back(s);
}


