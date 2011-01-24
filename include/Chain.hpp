/****************************************************************************************************/
//  COPYRIGHT 2011, University of Tennessee
//  Author: David Jenkins (david.d.jenkins@gmail.com)
//  File: Chain.hpp
//  Date: 12 Jan 2011
//  Version: 1.0
//  Description: Contains the class definition and implementation for Chain
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


#ifndef __CHAIN_HPP__
#define __CHAIN_HPP__

#include "AminoAcid.hpp"
#include "Atom.hpp"
#include "Seqres.hpp"

class Chain
{
public:
  // Constructor to reset everything
  Chain();
  // Constructor that sets the chain id
  Chain(char i);
  // Destructor to reset everything
  ~Chain();
  // Adds an AA to a vector
  void addAminoAcid(AminoAcid a);
  // Adds reference to a HETATM to a vector of Atom*
  void addHetatm(Atom* h);
  // Adds reference to a SEQRES to a vector of Seqres*
  void addSeqres(Seqres* s);
  
  char                  id;             // Chain id
  vector<AminoAcid>     aa;             // Vector of atoms in this chain
  vector<Atom*>         hetatms;        // Vector of hetatms in this chain
  vector<Seqres*>       seqres;         // Vector of seqres in this chain

  // Overloads the == operator
  inline bool operator==(const Chain &rhs) const
  {
    return this->id == rhs.id;
  }

};


#endif
