/****************************************************************************************************/
//  COPYRIGHT 2011, University of Tennessee
//  Author: David Jenkins (david.d.jenkins@gmail.com)
//  File: AminoAcid.cpp
//  Date: 16 Jan 2011
//  Date Modified: 7 Feb 2011
//  Version: 1.0
//  Description: Class implementations for AminoAcid
//
//  Updates: Fixed the non-center of charge calculations for ASP and GLU (7 Feb 2011)
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

#include "AminoAcid.hpp"
#include "Geometry.hpp"
#include "CoutColors.hpp"

AminoAcid::AminoAcid()
{
  atom.clear();
  center.clear();
  skip = false;
  altLoc = false;
  corrected = false;
}

AminoAcid::~AminoAcid()
{
  atom.clear();
  center.clear();
  skip = false;
  altLoc = false;
  corrected = false;
}

// All combinations of centers are calculated for every possible
// combination of locations.  For instance, if we have something like
// the following:
// ATOM    816  N   GLU    99      30.629  13.024  26.769 ...
// ATOM    817  CA  GLU    99      29.580  12.059  27.155 ...
// ATOM    818  C   GLU    99      30.037  10.639  26.931 ...
// ATOM    819  O   GLU    99      29.810   9.732  27.750 ...
// ATOM    820  CB  GLU    99      28.301  12.393  26.389 ...
// ATOM    821  CG  GLU    99      27.763  13.769  26.785 ...
// ATOM    822  CD AGLU    99      26.791  14.298  25.745 ...
// ATOM    823  CD BGLU    99      27.482  13.867  28.279 ...
// ATOM    824  OE1AGLU    99      26.375  13.533  24.848 ...
// ATOM    825  OE1BGLU    99      26.910  12.902  28.837 ...
// ATOM    826  OE2BGLU    99      27.834  14.887  28.915 ...
// ATOM    827  OE2AGLU    99      26.437  15.492  25.828 ...
// There would be a total of 8 combinations of possible centers:
// CG, CD A, OE1 A, OE2 A
// CG, CD A, OE1 A, OE2 B
// CG, CD A, OE1 B, OE2 A
// CG, CD A, OE1 B, OE2 B
// CG, CD B, OE1 A, OE2 A
// CG, CD B, OE1 A, OE2 B
// CG, CD B, OE1 B, OE2 A
// CG, CD B, OE1 B, OE2 B
// Thus, we store all of these possible centers in a vector to 
// examined later when we calculate the distances.
//////////////////////////////////////////////////////////////////////

// Calculate the centers for PHE or TYR
// these AA need the following atoms:
//   CG, CZ, CD1, CD2, CE1, CE2
// all of the other ones don't affect the center
// This also sets the plane coordinates used to calculate 
// the angle later in the program
void AminoAcid::centerPHEorTYR()
{
  // This makes a vector of a vector of ATOM* 
  // since we are looking for 6 atoms, the size of the outer vector
  // is 6 while the size of the inner vectors are dependent on the 
  // number of alternate locations for each one. Indexes are:
  // CG = 0 ; CZ = 1 ; CD1 = 2 ; CD2 = 3 ; CE1 = 4 ; CE2 = 5
  vector< vector <Atom*> > temp;
  temp.resize(6);

  // Push all of the important atoms in their respective vectors
  for(unsigned int i =0; i< atom.size(); i++)
    {
      if(atom[i]->altLoc != ' ' )
        {
          altLoc = true;
        }
      if(atom[i]->name == " CG ")
        {
          temp[0].push_back(atom[i]);
        }
      else if(atom[i]->name == " CZ ")
        {
          temp[1].push_back(atom[i]);
        }
      else if(atom[i]->name == " CD1")
        {
          temp[2].push_back(atom[i]);
        }
      else if(atom[i]->name == " CD2")
        {
          temp[3].push_back(atom[i]);
        }
      else if(atom[i]->name == " CE1")
        {
          temp[4].push_back(atom[i]);
        }
      else if(atom[i]->name == " CE2")
        {
          temp[5].push_back(atom[i]);
        }
      else
        {
          // This means that this atom is not useful so we flag it 
          atom[i]->skip = true;
        }
    }
  
  // Error check.  If we don't have at least one of each
  // of the atoms, we throw a warning, set an ignore flag
  // for future reference, and leave the function
  if(temp[0].size() == 0 || temp[1].size() == 0 ||
     temp[2].size() == 0 || temp[3].size() == 0 ||
     temp[4].size() == 0 || temp[5].size() == 0 )
    {
      skip = true;
#ifndef DISABLE_WARNING
      cout << cyan << "WARNING" << reset << ": Could not find all atoms in the TRP ring at "
           << atom[0]->resSeq << endl;
#endif
      return;
    }

  // Allocate size for all of the possible centers
  center.resize(temp[0].size() * temp[1].size() * temp[2].size() *
                temp[3].size() * temp[4].size() * temp[5].size());

  // And go through all combinations of the alternate locations
  unsigned int index = 0;
  for(unsigned int i = 0; i < temp[0].size(); i++) // go through all the CG
    {
      for(unsigned int j = 0; j < temp[1].size(); j++) // go through all the CZ
        {
          for(unsigned int k = 0; k < temp[2].size(); k++) // go through all the CD1
            {
              for(unsigned int l = 0; l < temp[3].size(); l++) // go through all the CD2
                {
                  for(unsigned int m = 0; m < temp[4].size(); m++) // go through all the CE1
                    {
                      for(unsigned int n = 0; n < temp[5].size(); n++, index++) // go through all the CE2
                        {
                          center[index].plane_info.resize(6);
                          center[index].plane_info[ CG_PLANE_COORD_PTT] = &(temp[0][i]->coord);
                          center[index].plane_info[CD1_PLANE_COORD_PTT] = &(temp[2][k]->coord);
                          center[index].plane_info[CD2_PLANE_COORD_PTT] = &(temp[3][l]->coord);
                          center[index].plane_info[3] = &(temp[1][i]->coord);
                          center[index].plane_info[4] = &(temp[4][m]->coord);
                          center[index].plane_info[5] = &(temp[5][n]->coord);

                          // This is just an average of the all the coordinates
                          center[index] = ( temp[0][i]->coord +
                                            temp[1][j]->coord +
                                            temp[2][k]->coord +
                                            temp[3][l]->coord +
                                            temp[4][m]->coord +
                                            temp[5][n]->coord )/(6);

                          // And this is just so we know what combination of
                          // alternate locations (may be able to take out later)
                          string t;
                          t.insert(t.end(),1,temp[0][i]->altLoc);
                          t.insert(t.end(),1,temp[1][j]->altLoc);
                          t.insert(t.end(),1,temp[2][k]->altLoc);
                          t.insert(t.end(),1,temp[3][l]->altLoc);
                          t.insert(t.end(),1,temp[4][m]->altLoc);
                          t.insert(t.end(),1,temp[5][n]->altLoc);
                          center[index].altLoc = t;
#ifdef DEBUG
                          cout << "CHECK: " << residue << atom[0]->resSeq << " : " << center[index] << endl;
#endif
                        }
                    }
                }
            }
        }
    }
}

// Calculate the centers for TRP
// these AA need the following atoms:
//   CG, CH2, CD1, CD2, NE1, CE2, CE3, CZ2, CZ3
// all of the other ones don't affect the center
// This also sets the plane coordinates used to calculate 
// the angle later in the program
void AminoAcid::centerTRP()
{
  // This makes a vector of a vector of ATOM* 
  // since we are looking for 9 atoms, the size of the outer vector
  // is 9 while the size of the inner vectors are dependent on the 
  // number of alternate locations for each one. Indexes are:
  // CG = 0  ; CH1 = 1 ; CD1 = 2 ; CD2 = 3 ; NE1 = 4 ; CE2 = 5 ;
  // CE3 = 6 ; CZ2 = 7 ; CZ3 = 8
  vector< vector <Atom*> > temp;
  temp.resize(9);

  // Push all of the important atoms in their respective vectors
  for(unsigned int i =0; i< atom.size(); i++)
    {
      if(atom[i]->altLoc != ' ' )
        {
          altLoc = true;
        }
      if(atom[i]->name == " CG ")
        {
          temp[0].push_back(atom[i]);
        }
      else if(atom[i]->name == " CH2")
        {
          temp[1].push_back(atom[i]);
        }
      else if(atom[i]->name == " CD1")
        {
          temp[2].push_back(atom[i]);
        }
      else if(atom[i]->name == " CD2")
        {
          temp[3].push_back(atom[i]);
        }
      else if(atom[i]->name == " NE1")
        {
          temp[4].push_back(atom[i]);
        }
      else if(atom[i]->name == " CE2")
        {
          temp[5].push_back(atom[i]);
        }
      else if(atom[i]->name == " CE3")
        {
          temp[6].push_back(atom[i]);
        }
      else if(atom[i]->name == " CZ2")
        {
          temp[7].push_back(atom[i]);
        }
      else if(atom[i]->name == " CZ3")
        {
          temp[8].push_back(atom[i]);
        }
      else
        {
          // This means that this atom is not useful so we flag it 
          atom[i]->skip = true;
        }
    }

  // Error check.  If we don't have at least one of each
  // of the atoms, we throw a warning, set an ignore flag
  // for future reference, and leave the function
  if(temp[0].size() == 0 || temp[1].size() == 0 ||
     temp[2].size() == 0 || temp[3].size() == 0 ||
     temp[4].size() == 0 || temp[5].size() == 0 ||
     temp[6].size() == 0 || temp[7].size() == 0 ||
     temp[8].size() == 0)
    {
      skip = true;
#ifndef DISABLE_WARNING
      cout << cyan << "WARNING" << reset << ": Could not find all atoms in the TRP ring at "
           << atom[0]->resSeq << endl;
#endif
      return;
    }
  // Allocate size for all of the possible centers
  center.resize( temp[0].size() * temp[1].size() * temp[2].size() * temp[3].size()*
                 temp[4].size() * temp[5].size() * temp[6].size() * temp[7].size() * temp[8].size());

  // And go through all combinations of the alternate locations
  // (Sorry that it is really freaking ugly)
  unsigned int index = 0;
  for(unsigned int i = 0; i < temp[0].size(); i++) // go through all the CG
    {
      for(unsigned int j = 0; j < temp[1].size(); j++) // go through all the CH2
        {
          for(unsigned int k = 0; k < temp[2].size(); k++) // go through all the CD1
            {
              for(unsigned int l = 0; l < temp[3].size(); l++) // go through all the CD2
                {
                  for(unsigned int m = 0; m < temp[4].size(); m++) // go through all the NE1
                    {
                      for(unsigned int n = 0; n < temp[5].size(); n++) // go through all the CE2
                        {
                          for(unsigned int o = 0; o < temp[6].size(); o++) // go through all the CE3
                            {
                              for(unsigned int p = 0; p < temp[7].size(); p++) // go through all the CZ2
                                {
                                  for(unsigned int q = 0; q < temp[8].size(); q++, index++) // go through all the CZ2
                                    {
                                      center[index].plane_info.resize(3);
                                      center[index].plane_info[CG_PLANE_COORD_PTT]  = &temp[0][i]->coord;
                                      center[index].plane_info[CD1_PLANE_COORD_PTT] = &temp[2][k]->coord;
                                      center[index].plane_info[CD2_PLANE_COORD_PTT] = &temp[3][l]->coord;

                                      // This is a weighted average of the atoms weighted
                                      // by their mass
                                      center[index] = ( temp[0][i]->coord * MASS_C +
                                                        temp[1][j]->coord * MASS_C +
                                                        temp[2][k]->coord * MASS_C +
                                                        temp[3][l]->coord * MASS_C +
                                                        temp[4][m]->coord * MASS_N +
                                                        temp[5][n]->coord * MASS_C +
                                                        temp[6][o]->coord * MASS_C +
                                                        temp[7][p]->coord * MASS_C +
                                                        temp[8][q]->coord * MASS_C )/(MASS_C*8+MASS_N);

                                      // And this is just so we know what combination of
                                      // alternate locations (may be able to take out later)
                                      string t;
                                      t.insert(t.end(),1,temp[0][i]->altLoc);
                                      t.insert(t.end(),1,temp[1][j]->altLoc);
                                      t.insert(t.end(),1,temp[2][k]->altLoc);
                                      t.insert(t.end(),1,temp[3][l]->altLoc);
                                      t.insert(t.end(),1,temp[4][m]->altLoc);
                                      t.insert(t.end(),1,temp[5][n]->altLoc);
                                      t.insert(t.end(),1,temp[6][o]->altLoc);
                                      t.insert(t.end(),1,temp[7][p]->altLoc);
                                      t.insert(t.end(),1,temp[8][q]->altLoc);
                                      center[index].altLoc = t;
#ifdef DEBUG
                                      cout << "CHECK: TRP" << atom[0]->resSeq << " : " << center[index] << endl;
#endif
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }
}

// Calculate the centers for LYS
// these AA need the following atoms:
//   CE, NZ, and (for now) CD
// all of the other ones don't affect the center
// This also sets the plane coordinates used to calculate 
// the angle later in the program
void AminoAcid::centerLYS()
{
  // This makes a vector of a vector of ATOM*
  // since we are looking for 3 atoms, the size of the outer vector
  // is 3 while the size of the inner vectors are dependent on the
  // number of alternate locations for each one. Indexes are:
  // CD = 0 ; CE = 1 ; NZ = 2
  vector< vector <Atom*> > temp;
  temp.resize(3);

  // Push all of the important atoms in their respective vectors
  for (unsigned int i = 0; i < atom.size(); i++)
    {
      if (atom[i]->altLoc != ' ')
	{
	  altLoc = true;
	}
// The " CD " atom is only here to give three points needed
// for a plane to be formed. Later, it will be used to find
// the locations of the hydrogens and have it's position be
// a location for a hydrogen.
      if (atom[i]->name == " CD ")
	{
	  temp[2].push_back(atom[i]);
	  atom[i]->skip = true;
	}
      else if (atom[i]->name == " CE ")
	{
	  temp[0].push_back(atom[i]);
	}
      else if (atom[i]->name == " NZ ")
	{
	  temp[1].push_back(atom[i]);
	}
      else
	{
	  // This means that this atom is not useful so we flag it
	  atom[i]->skip = true;
	}
    }

  // Error check. If we don't have at least one of each
  // of the atoms, we throw a warning, set an ignore flag
  // for future reference, and leave the function
  if (temp[0].size() == 0 || temp[1].size() == 0 || temp[2].size() == 0)
    {
      skip = true;
#ifndef DISABLE_WARNING
      cout << cyan << "WARNING" << reset << ": Could not find all atoms in the TRP ring at "
	   << atom[0]->resSeq << endl;
#endif
      return;
    }
  // Allocate size for all of the possible centers
  center.resize( temp[0].size() * temp[1].size() * temp[2].size());

  // And go through all combinations of the alternate locations
  unsigned int index = 0;
  for (unsigned int i = 0; i < temp[0].size(); i++) // Go through all the CE
    {
      for (unsigned int j = 0; j < temp[1].size(); j++) // Go through all the NZ
	{
	  for (unsigned int k = 0; k < temp[2].size(); k++, index++) // Go through all the CD
	    {
	      center[index].plane_info.resize(3);
	      center[index].plane_info[0] = &temp[0][i]->coord;
	      center[index].plane_info[1] = &temp[1][j]->coord;
	      center[index].plane_info[2] = &temp[2][k]->coord;

	      // This is a weighted average of the NZ and CE atoms 
	      // weighted by their mass
	        center[index] = ( temp[0][i]->coord * MASS_C +
	  		          temp[1][j]->coord * MASS_N ) / (MASS_C + MASS_N);

	      // And this is just so we know what combination of
	      // alternate locations (may be able to take out later)
	      string t;
	      t.insert(t.end(), 1, temp[0][i]->altLoc);
	      t.insert(t.end(), 1, temp[1][j]->altLoc);
	      center[index].altLoc = t;
#ifdef DEBUG
	      cout << "CHECK: LYS" << atom[0]->resSeq << " : " << center[index] << endl;
#endif
	    }
	}
    }
}

// Calculate the centers for ARG
// these AA need the following atoms:
//   CD, NE, CZ, NH1, and NH2
// all of the other ones don't affect the center
// This also sets the plane coordinates used to calculate 
// the angle later in the program
void AminoAcid::centerARG()
{
  // This makes a vector of a vector of ATOM*
  // since we are looking for 5 atoms, the size of the outer vector
  // is 5 while the size of the inner vectors are dependent on the
  // number of alternate locations for each one. Indexes are:
  // CD = 0 ; NE = 1 ; CZ = 2 ; NH1 = 3 ; NH2 = 4
  vector< vector <Atom*> > temp;
  temp.resize(5);

  // Push all of the important atoms in their respective vectors
  for (unsigned int i = 0; i < atom.size(); i++)
    {
      if (atom[i]->altLoc != ' ' )
	{
	  altLoc = true;
	}
      if (atom[i]->name == " CD ")
	{
	  temp[0].push_back(atom[i]);
	}
      else if (atom[i]->name == " NE ")
	{
	  temp[1].push_back(atom[i]);
	}
      else if (atom[i]->name == " CZ ")
	{
	  temp[2].push_back(atom[i]);
	}
      else if (atom[i]->name == " NH1")
	{
	  temp[3].push_back(atom[i]);
	}
      else if (atom[i]->name == " NH2")
	{
	  temp[4].push_back(atom[i]);
	}
      else
	{
	  // This means that this atom is not useful so we flag it
	  atom[i]->skip = true;
	}
    }

  // Error check. If we don't have at least one of each
  // of the atoms, we throw a warning, set an ignore flag
  // for future reference, and leave the function
  if (temp[0].size() == 0 || temp[1].size() == 0 ||
      temp[2].size() == 0 || temp[3].size() == 0 ||
      temp[4].size() == 0)
    {
      skip = true;
#ifndef DISABLE_WARNING
      cout << cyan << "WARNING" << reset << ": Could not find all atoms in the ARG ring at "
	   << atom[0]->resSeq << endl;
#endif
      return;
    }

  // Allocate size for all of the possible centers
  center.resize( temp[0].size() * temp[1].size() * temp[2].size() * temp[3].size() *
		 temp[4].size());

  // And go through all combinations of the alternate locations
  unsigned int index = 0;
  for (unsigned int i = 0; i < temp[0].size(); i++) // Go through all the CD
    {
      for (unsigned int j = 0; j < temp[1].size(); j++) // Go through all the NE
	{
	  for (unsigned int k = 0; k < temp[2].size(); k++) // Go through all the CZ
	    {
	      for (unsigned int l = 0; l < temp[3].size(); l++) // Go through all the NH1
		{
		  for (unsigned int m = 0; m < temp[4].size(); m++, index++) // Go through all the NH2
		    {
		      center[index].plane_info.resize(5);
		      center[index].plane_info[0] = &temp[0][i]->coord;
		      center[index].plane_info[1] = &temp[1][j]->coord;
		      center[index].plane_info[2] = &temp[2][k]->coord;
		      center[index].plane_info[3] = &temp[3][l]->coord;
		      center[index].plane_info[4] = &temp[4][m]->coord;

		      // This is a weighted average of the atoms weighted
		      // by their mass
		      center[index] = ( temp[0][i]->coord * MASS_C +
					temp[1][j]->coord * MASS_N +
					temp[2][k]->coord * MASS_C +
					temp[3][l]->coord * MASS_N +
					temp[4][m]->coord * MASS_N )/((MASS_C*2) + (MASS_N*3));

		      // And this is just so we know what combination of
		      // alternate locations (may be able to take out later)
		      string t;
		      t.insert(t.end(), 1, temp[0][i]->altLoc);
		      t.insert(t.end(), 1, temp[1][j]->altLoc);
		      t.insert(t.end(), 1, temp[2][k]->altLoc);
		      t.insert(t.end(), 1, temp[3][l]->altLoc);
		      t.insert(t.end(), 1, temp[4][m]->altLoc);
		      center[index].altLoc = t;
#ifdef DEBUG
		      cout << "CHECK: ARG" << atom[0]->resSeq << " : " << center[index] << endl;
#endif
		    }
		}
	    }
	}
    }
}

// DO NOT USE THIS!!!! If you use this, a fairy will lose its wings!
// This is only here for legacy reasons because this is what the old
// STAAR code used.  It is wrong.
// Calculate the centers of charge for ASP
// these AA need the following atoms:
//   CG, CB, OD1, OD2
// all of the other ones don't affect the center
// This also sets the plane coordinates used to calculate 
// the angle later in the program
void AminoAcid::centerASP()
{
  // This makes a vector of a vector of ATOM* 
  // since we are looking for 4 atoms, the size of the outer vector
  // is 4 while the size of the inner vectors are dependent on the 
  // number of alternate locations for each one. Indexes are:
  // CG = 0 ; CB = 1 ; OD1 = 2 ; OD2 = 3 ;
  vector< vector <Atom*> > temp;
  temp.resize(4);

  // Push all of the important atoms in their respective vectors
  for(unsigned int i =0; i< atom.size(); i++)
    {
      if(atom[i]->altLoc != ' ')
        {
          altLoc = true;
        }
      if(atom[i]->name == " CG ")
        {
          temp[0].push_back(atom[i]);
        }
      else if(atom[i]->name == " H  ")
        {
          temp[1].push_back(atom[i]);
        }
      else if(atom[i]->name == " OD1")
        {
          temp[2].push_back(atom[i]);
        }
      else if(atom[i]->name == " OD2")
        {
          temp[3].push_back(atom[i]);
        }
      else
        {
          // This means that this atom is not useful so we flag it 
          atom[i]->skip = true;
        }
    }

  // Error check.  If we don't have at least one of each
  // of the atoms, we throw a warning, set an ignore flag
  // for future reference, and leave the function
  if(temp[0].size() == 0 || temp[1].size() == 0 ||
     temp[2].size() == 0 || temp[3].size() == 0)
    {
      skip = true;
#ifndef DISABLE_WARNING
      cout << cyan << "WARNING" << reset << ": Could not find all atoms in the ASP side chain at "
           << atom[0]->resSeq << endl;
#endif
      return;
    }

  // Allocate size for all of the possible centers
  center.resize(temp[0].size() * temp[1].size() * temp[2].size() * temp[3].size());

  // And go through all combinations of the alternate locations
  unsigned int index = 0;
  for(unsigned int i = 0; i < temp[0].size(); i++) // go through all the CG
    {
      for(unsigned int j = 0; j < temp[1].size(); j++) // go through all the CB
        {
          for(unsigned int k = 0; k < temp[2].size(); k++) // go through all the OD1
            {
              for(unsigned int l = 0; l < temp[3].size(); l++, index++) // go through all the OD2
                {
                  center[index].plane_info.resize(3);
                  center[index].plane_info[C__PLANE_COORD_AG]  = &temp[0][i]->coord;
                  center[index].plane_info[O_1_PLANE_COORD_AG] = &temp[2][k]->coord;
                  center[index].plane_info[O_2_PLANE_COORD_AG] = &temp[3][l]->coord;

                  // This is just an average of the all the coordinates
                  center[index] = ( temp[0][i]->coord * CHARGE_C +
                                    temp[1][j]->coord * CHARGE_H + // treating CB as H in formic acid
                                    temp[2][k]->coord * CHARGE_O +
                                    temp[3][l]->coord * CHARGE_O ) * -1;
                  
                  // And this is just so we know what combination of
                  // alternate locations (may be able to take out later)
                  string t;
                  t.insert(t.end(),1,temp[0][i]->altLoc);
                  t.insert(t.end(),1,temp[1][j]->altLoc);
                  t.insert(t.end(),1,temp[2][k]->altLoc);
                  t.insert(t.end(),1,temp[3][l]->altLoc);
                  center[index].altLoc = t;
#ifdef DEBUG
                  cout << "CHECK: ASP" << atom[0]->resSeq << " : " << center[index] << endl;
#endif
                }
            }
        }
    }
}

// DO NOT USE THIS!!!! If you use this, puppies will be harmed!
// This is only here for legacy reasons because this is what the old
// STAAR code used.  It is wrong.  
// Calculate the centers of charge for GLU
// these AA need the following atoms:
//   CG, CD, OE1, OE2
// all of the other ones don't affect the center
// This also sets the plane coordinates used to calculate 
// the angle later in the program
void AminoAcid::centerGLU()
{
  // This makes a vector of a vector of ATOM* 
  // since we are looking for 4 atoms, the size of the outer vector
  // is 4 while the size of the inner vectors are dependent on the 
  // number of alternate locations for each one. Indexes are:
  // CG = 0 ; CD = 1 ; OE1 = 2 ; OE2 = 3
  vector< vector <Atom*> > temp;
  temp.resize(4);

  // Push all of the important atoms in their respective vectors  
  for(unsigned int i =0; i< atom.size(); i++)
    {
      if(atom[i]->altLoc != ' ')
        {
          altLoc = true;
        }
      if(atom[i]->name == " H  ")
        {
          temp[0].push_back(atom[i]);
        }
      else if(atom[i]->name == " CD ")
        {
          temp[1].push_back(atom[i]);
        }
      else if(atom[i]->name == " OE1")
        {
          temp[2].push_back(atom[i]);
        }
      else if(atom[i]->name == " OE2")
        {
          temp[3].push_back(atom[i]);
        }
      else
        {
          // This means that this atom is not useful so we flag it 
          atom[i]->skip = true;
        }
    }

  // Error check.  If we don't have at least one of each
  // of the atoms, we throw a warning, set an ignore flag
  // for future reference, and leave the function
  if(temp[0].size() == 0 || temp[1].size() == 0 ||
     temp[2].size() == 0 || temp[3].size() == 0)
    {
      skip = true;
#ifndef DISABLE_WARNING
      cout << cyan << "WARNING" << reset << ": Could not find all atoms in the GLU side chain at "
           << atom[0]->resSeq << endl;
#endif
      return;
    }

  // Allocate size for all of the possible centers
  center.resize(temp[0].size() * temp[1].size() * temp[2].size() * temp[3].size());

    // And go through all combinations of the alternate locations
  unsigned int index = 0;
  for(unsigned int i = 0; i < temp[0].size(); i++) // go through all the CG
    {
      for(unsigned int j = 0; j < temp[1].size(); j++) // go through all the CD
        {
          for(unsigned int k = 0; k < temp[2].size(); k++) // go through all the OE1
            {
              for(unsigned int l = 0; l < temp[3].size(); l++, index++) // go through all the OE2
                {

                  center[index].plane_info.resize(3);
                  center[index].plane_info[C__PLANE_COORD_AG]  = &(temp[1][j]->coord);
                  center[index].plane_info[O_1_PLANE_COORD_AG] = &(temp[2][k]->coord);
                  center[index].plane_info[O_2_PLANE_COORD_AG] = &(temp[3][l]->coord);
                  
                  // This is a weighted average weighted by the charges.
                  // Since the charges sum to -1, we just need to negate the signs
                  center[index] = ( temp[0][i]->coord * CHARGE_H + // CG treated as H in formic acid
                                    temp[1][j]->coord * CHARGE_C +
                                    temp[2][k]->coord * CHARGE_O +
                                    temp[3][l]->coord * CHARGE_O ) * -1;
                  
                  // And this is just so we know what combination of
                  // alternate locations (may be able to take out later)
                  string t;
                  t.insert(t.end(),1,temp[0][i]->altLoc);
                  t.insert(t.end(),1,temp[1][j]->altLoc);
                  t.insert(t.end(),1,temp[2][k]->altLoc);
                  t.insert(t.end(),1,temp[3][l]->altLoc);
                  center[index].altLoc = t;
#ifdef DEBUG
                  cout << "CHECK: GLU" << atom[0]->resSeq << " : " << center[index] << endl;
#endif
                }
            }
        }
    }
}

// Calculate the centers for ASP
// these AA need the following atoms:
//   OD1 and OD2
// all of the other ones don't affect the center
// This is is really simple: all it does is add the
// OD1 and OD2 coords to the center vector.
// This also sets the plane coordinates used to calculate 
// the angle later in the program
void AminoAcid::centerASP_oxygen()
{
  // This makes a vector of a vector of ATOM* 
  // since we are looking for 3 atoms, the size of the outer vector
  // is 3 while the size of the inner vectors are dependent on the 
  // number of alternate locations for each one. Indexes are:
  // CG = 0 ; OD1 = 1 ; OD2 = 2
  vector< vector <Atom*> > temp;
  temp.resize(3);

  // Push all of the important atoms in their respective vectors  
  for(unsigned int i =0; i< atom.size(); i++)
    {
      if(atom[i]->altLoc != ' ')
        {
          altLoc = true;
        }
      if(atom[i]->name == " CG ")
        {
          temp[0].push_back(atom[i]);
        }
      else if(atom[i]->name == " OD1")
        {
          temp[1].push_back(atom[i]);
        }
      else if(atom[i]->name == " OD2")
        {
          temp[2].push_back(atom[i]);
        }
      else
        {
          // This means that this atom is not useful so we flag it 
          atom[i]->skip = true;
        }
    }

  // Error check.  If we don't have at least one of each
  // of the atoms, we throw a warning, set an ignore flag
  // for future reference, and leave the function
  if( temp[1].size() == 0 ||
      temp[2].size() == 0)
    {
      skip = true;
#ifndef DISABLE_WARNING
      cout << cyan << "WARNING" << reset << ": Could not find all atoms in the ASP side chain at "
           << atom[0]->resSeq << endl;
#endif
      return;
    }
  unsigned int offset = temp[0].size() * temp[1].size() * temp[2].size();
  // Allocate size for all of the possible centers
  center.resize(offset * 2);

  // And go through all combinations of the alternate locations
  unsigned int index = 0;
  for(unsigned int k = 0; k < temp[1].size(); k++) // go through all the OD1
    {      
      for(unsigned int j = 0; j < temp[0].size(); j++) // go through all the CD
        {
          for(unsigned int l = 0; l < temp[2].size(); l++, index++) // go through all the OD2
            {
              
              center[index].plane_info.resize(3);
              center[index].plane_info[C__PLANE_COORD_AG]  = &(temp[0][j]->coord);
              center[index].plane_info[O_1_PLANE_COORD_AG] = &(temp[1][k]->coord);
              center[index].plane_info[O_2_PLANE_COORD_AG] = &(temp[2][l]->coord);
              center[index + offset].plane_info.resize(3);
              center[index + offset].plane_info[C__PLANE_COORD_AG]  = &(temp[0][j]->coord);
              center[index + offset].plane_info[O_1_PLANE_COORD_AG] = &(temp[1][k]->coord);
              center[index + offset].plane_info[O_2_PLANE_COORD_AG] = &(temp[2][l]->coord);
                  
              center[index] = temp[1][k]->coord;
              center[index + offset] = temp[2][l]->coord;
                  
              // And this is just so we know what combination of
              // alternate locations (may be able to take out later)
              string t;
              t.insert(t.end(),1,temp[0][j]->altLoc);
              t.insert(t.end(),1,temp[1][k]->altLoc);
              t.insert(t.end(),1,temp[2][l]->altLoc);
              center[index].altLoc = t;
              center[index + offset].altLoc = t;
#ifdef DEBUG
              cout << "CHECK: ASP" << atom[0]->resSeq << " : " << center[index] << endl;
#endif
            }
        }
    }
}

// Calculate the centers for GLU
// these AA need the following atoms:
//   OE1 and OE2
// all of the other ones don't affect the center
// This is is really simple: all it does is add the
// OE1 and OE2 coords to the center vector.
// This also sets the plane coordinates used to calculate 
// the angle later in the program
void AminoAcid::centerGLU_oxygen()
{
  // This makes a vector of a vector of ATOM* 
  // since we are looking for 3 atoms, the size of the outer vector
  // is 3 while the size of the inner vectors are dependent on the 
  // number of alternate locations for each one. Indexes are:
  // CD = 0 ; OE1 = 1 ; OE2 = 2
  vector< vector <Atom*> > temp;
  temp.resize(3);

  // Push all of the important atoms in their respective vectors  
  for(unsigned int i =0; i< atom.size(); i++)
    {
      if(atom[i]->altLoc != ' ')
        {
          altLoc = true;
        }
      if(atom[i]->name == " CD ")
        {
          temp[0].push_back(atom[i]);
        }
      else if(atom[i]->name == " OE1")
        {
          temp[1].push_back(atom[i]);
        }
      else if(atom[i]->name == " OE2")
        {
          temp[2].push_back(atom[i]);
        }
      else
        {
          // This means that this atom is not useful so we flag it 
          atom[i]->skip = true;
        }
    }

  // Error check.  If we don't have at least one of each
  // of the atoms, we throw a warning, set an ignore flag
  // for future reference, and leave the function
  if(temp[0].size() == 0 || temp[1].size() == 0 ||
     temp[2].size() == 0)
    {
      skip = true;
#ifndef DISABLE_WARNING
      cout << cyan << "WARNING" << reset << ": Could not find all atoms in the GLU side chain at "
           << atom[0]->resSeq << endl;
#endif
      return;
    }
  unsigned int offset = temp[0].size() * temp[1].size() * temp[2].size();
  // Allocate size for all of the possible centers
  center.resize(offset * 2);

  // And go through all combinations of the alternate locations
  unsigned int index = 0;
  for(unsigned int k = 0; k < temp[1].size(); k++) // go through all the OE1
    {
      for(unsigned int j = 0; j < temp[0].size(); j++) // go through all the CD
        {
          for(unsigned int l = 0; l < temp[2].size(); l++, index++) // go through all the OE2
            {
              
              center[index].plane_info.resize(3);
              center[index].plane_info[C__PLANE_COORD_AG]  = &(temp[0][j]->coord);
              center[index].plane_info[O_1_PLANE_COORD_AG] = &(temp[1][k]->coord);
              center[index].plane_info[O_2_PLANE_COORD_AG] = &(temp[2][l]->coord);
              center[index + offset].plane_info.resize(3);
              center[index + offset].plane_info[C__PLANE_COORD_AG]  = &(temp[0][j]->coord);
              center[index + offset].plane_info[O_1_PLANE_COORD_AG] = &(temp[1][k]->coord);
              center[index + offset].plane_info[O_2_PLANE_COORD_AG] = &(temp[2][l]->coord);
                  
              center[index] = temp[1][k]->coord;
              center[index + offset] = temp[2][l]->coord;
                  
              // And this is just so we know what combination of
              // alternate locations (may be able to take out later)
              string t;
              t.insert(t.end(),1,temp[0][j]->altLoc);
              t.insert(t.end(),1,temp[1][k]->altLoc);
              t.insert(t.end(),1,temp[2][l]->altLoc);
              center[index].altLoc = t;
              center[index + offset].altLoc = t;
#ifdef DEBUG
              cout << "CHECK: GLU" << atom[0]->resSeq << " : " << center[index] << endl;
#endif
            }
        }
    }
}

void AminoAcid::centerPHEorTYR_simplified()
{
  center.resize(1);

  // Temp variable that will have the center coordinates
  Coordinates tempCenter(0.0, 0.0, 0.0);

  // Here, we are just getting 2 C atoms. Doesn't matter
  // which ones we choose as long as they are
  // diametrically opposed. Based off of 
  // http://comp.chem.nottingham.ac.uk/dichrocalc/files/atomlabels-sidechains.png
  // the following are diametically opposed:
  //   CG  - CZ
  //   CD2 - CE1
  //   CD1 - CE2
  center[0].set(0.0, 0.0, 0.0);
  center[0].plane_info.resize(3);
  for(unsigned int i =0; i< atom.size(); i++)
    {
      if(atom[i]->name == " CE2")
        {
          center[0] += atom[i]->coord;
          center[0].plane_info[CE2_PLANE_COORD_PTT] = &atom[i]->coord;
        }
      else if(atom[i]->name == " CD1")
        {
          center[0] += atom[i]->coord;
          center[0].plane_info[CD1_PLANE_COORD_PTT] = &atom[i]->coord;
        }
      else if(atom[i]->name == " CG ")
        {
          center[0].plane_info[CG_PLANE_COORD_PTT] = &atom[i]->coord;
        }
    }
  center[0] /= 2;
}

void AminoAcid::centerASP_charge()
{
  center.resize(1);
  Coordinates tempCenter(0.0, 0.0, 0.0);
  for(unsigned int i =0; i< atom.size(); i++)
    {
      if(atom[i]->name == " CG ")
        {
          center[0] = atom[i]->coord;
          tempCenter += atom[i]->coord;
        }
      else if(atom[i]->name == " H  ")
        {
          tempCenter -= atom[i]->coord;
        }
    }
  center[0] += (tempCenter * HYDROGEN_BOND_DISTANCE) / tempCenter.norm();
}

void AminoAcid::centerLYS_charge()
{
  center.resize(1);
  Coordinates tempA;
  Coordinates tempB;
  for (unsigned int i = 0; i < atom.size(); i++)
    {
      if (atom[i]->name == " NZ " && !atom[i]->skip)
	{
	  tempA = atom[i]->coord;
	}
      if (atom[i]->name == " CE " && !atom[i]->skip)
	{
	  tempB = atom[i]->coord;
	}
    }
  center[0] = ((tempA * 82.7) + (tempB * 17.3)) / 100;
}

void AminoAcid::centerARG_charge()
{
  center.resize(1);
  Coordinates tempCenter(0.0, 0.0, 0.0);
  for (unsigned int i = 0; i < atom.size(); i++)
    {
      if (atom[i]->name == " CZ ")
	{
	  center[0] = atom[i]->coord;
	}
    }
}

void AminoAcid::centerGLU_charge()
{
  center.resize(1);
  Coordinates tempCenter(0.0, 0.0, 0.0);
  for(unsigned int i =0; i< atom.size(); i++)
    {
      if(atom[i]->name == " CD ")
        {
          center[0] = atom[i]->coord;
          tempCenter += atom[i]->coord;
        }
      else if(atom[i]->name == " H  ")
        {
          tempCenter -= atom[i]->coord;
        }
    }
  center[0] += (tempCenter * HYDROGEN_BOND_DISTANCE) / tempCenter.norm();
}


// Calculate the centers for PO4, 2HP, PI
// these AA need the following atoms:
//   P, O1, O2, O3, and O4
void AminoAcid::centerPO4or2HPorPI()
{
  // This makes a vector of a vector of ATOM* 
  // since we are looking for 5 atoms, the size of the outer vector
  // is 5 while the size of the inner vectors are dependent on the 
  // number of alternate locations for each one. Indexes are:
  // P = 0 ; O1 = 1 ; O2 = 2 ; O3 = 3 ; O4 = 4
  vector< vector <Atom*> > temp;
  temp.resize(5);

  // Push all of the important atoms in their respective vectors  
  for(unsigned int i =0; i< atom.size(); i++)
    {
      if(atom[i]->altLoc != ' ')
        {
          altLoc = true;
        }
      if(atom[i]->name == " P  ")
        {
          temp[0].push_back(atom[i]);
        }
      else if(atom[i]->name == " O1 ")
        {
          temp[1].push_back(atom[i]);
        }
      else if(atom[i]->name == " O2 ")
        {
          temp[2].push_back(atom[i]);
        }
      else if(atom[i]->name == " O3 ")
        {
          temp[3].push_back(atom[i]);
        }
      else if(atom[i]->name == " O4 ")
        {
          temp[4].push_back(atom[i]);
        }
      else
        {
          // This means that this atom is not useful so we flag it 
          atom[i]->skip = true;
        }
    }

  // Error check.  If we don't have at least one of each
  // of the atoms, we throw a warning, set an ignore flag
  // for future reference, and leave the function
  if(temp[0].size() == 0 || temp[1].size() == 0 ||
     temp[2].size() == 0 || temp[3].size() == 0 ||
     temp[4].size() == 0 )
    {
      skip = true;
#ifndef DISABLE_WARNING
      cout << cyan << "WARNING" << reset << ": Could not find all atoms in the " << residue <<  " at "
           << atom[0]->resSeq << " chain " << atom[0]->chainID << endl;
#endif
      return;
    }

  // Allocate size for all of the possible centers
  center.resize( temp[0].size() * temp[1].size() * temp[2].size() * temp[3].size() * temp[4].size() );

  // And go through all combinations of the alternate locations
  unsigned int index = 0;
  float totalmass = MASS_P + MASS_O * 4;
  for(unsigned int j = 0; j < temp[0].size(); j++) // go through all the P
    {
      for(unsigned int k = 0; k < temp[1].size(); k++) // go through all the O1
        {
          for(unsigned int l = 0; l < temp[2].size(); l++) // go through all the O2
            {
              for(unsigned int m = 0; m < temp[3].size(); m++) // go through all the O3
                {
                  for(unsigned int n = 0; n < temp[4].size(); n++, index++) // go through all the O4
                    {
                      center[index] = temp[0][j]->coord * MASS_P +
                        ( temp[1][k]->coord + temp[2][l]->coord +
                          temp[3][m]->coord + temp[4][n]->coord ) * MASS_O;

                      center[index] /= totalmass;

                      center[index].plane_info.resize(5);
                      center[index].plane_info[0] = &(temp[0][j]->coord);
                      center[index].plane_info[1] = &(temp[1][k]->coord);
                      center[index].plane_info[2] = &(temp[2][l]->coord);
                      center[index].plane_info[3] = &(temp[3][m]->coord);
                      center[index].plane_info[4] = &(temp[4][n]->coord);

#ifdef DEBUG
                      cout << "CHECK: " << residue  << " " << atom[0]->resSeq << " : " << center[index] << endl;
#endif
                    }
                }
            }
        }
    }
}

void AminoAcid::centerPO4or2HPorPI_charge()
{
  // This makes a vector of a vector of ATOM* 
  // since we are looking for 6 atoms, the size of the outer vector
  // is 3 while the size of the inner vectors are dependent on the 
  // number of alternate locations for each one. Indexes are:
  // P = 0 ; O1 = 1 ; O2 = 2 ; O3 = 3 ; O4 = 4
  Atom* P = NULL;
  vector< Atom* > O;
  vector< Atom* > H;

  // Push all of the important atoms in their respective vectors  
  for(unsigned int i =0; i< atom.size(); i++)
    {
      if(atom[i]->altLoc != ' ')
        {
          altLoc = true;
        }
      if(atom[i]->name == " P  ")
        {
          P = atom[i];
        }
      else if(atom[i]->name == " O1 ")
        {
          O.push_back(atom[i]);
        }
      else if(atom[i]->name == " O2 ")
        {
          O.push_back(atom[i]);
        }
      else if(atom[i]->name == " O3 ")
        {
          O.push_back(atom[i]);
        }
      else if(atom[i]->name == " O4 ")
        {
          O.push_back(atom[i]);
        }
      else if(atom[i]->name == " H  ")
        {
          H.push_back(atom[i]);
        }
      else
        {
          // This means that this atom is not useful so we flag it 
          atom[i]->skip = true;
        }
    }

  // Error check to see if we have 1 P, 4 O, and at least 1 H
  if(P == NULL || O.size() != 4 || H.size() == 0 )
    {
      skip = true;
#ifndef DISABLE_WARNING
      cout << cyan << "WARNING" << reset << ": Could not find all atoms in the " << residue <<  " at "
           << atom[0]->resSeq << " chain " << atom[0]->chainID << endl;
#endif
      cout << cyan << "Skipping " << reset << " because there weren't enough hydrogens in " 
           << residue << " " << atom[0]->resSeq << " chain " << atom[0]->chainID << endl;
      return;
    }

  center.resize( 1 );
  // Here is an estimate of the center by doing a charge weighted average
  // This could be wrong, but it is just something to put in here for the
  // time being until I talk with Dr. Hinde.
  center[0] = P->coord * CHARGE_P;
  Coordinates temp(0.0, 0.0, 0.0);
  temp  = O[0]->coord;
  temp += O[1]->coord;
  temp += O[2]->coord;
  temp += O[3]->coord;
  center[0] += temp*CHARGE_O;
  temp.set(0.0,0.0,0.0);
  for(unsigned int i=0; i<H.size(); i++)
    {
      temp += H[i]->coord;
    }
  center[0] += temp * CHARGE_H;

  // And since the formal charge of PO4 is -3, divide by it
  if(residue == "PO4")
    center[0] /= -3;
  else if(residue == "2HP")
    center[0] *= -1;
  else if(residue == " PI")
    {
      center[0] /= -2;
    }
}

// Calculate the centers for 2PO and PO3
// these AA need the following atoms:
//   P, O1(P), O2(P), O3(P)
void AminoAcid::center2POorPO3()
{
  // This makes a vector of a vector of ATOM* 
  // since we are looking for 5 atoms, the size of the outer vector
  // is 5 while the size of the inner vectors are dependent on the 
  // number of alternate locations for each one. Indexes are:
  // P = 0 ; O1(P) = 1 ; O2(P) = 2 ; O3(P) = 3 
  vector< vector <Atom*> > temp;
  temp.resize(4);

  // Push all of the important atoms in their respective vectors  
  for(unsigned int i =0; i< atom.size(); i++)
    {
      if(atom[i]->altLoc != ' ')
        {
          altLoc = true;
        }
      if(atom[i]->name == " P  ")
        {
          temp[0].push_back(atom[i]);
        }
      else if(atom[i]->name == " O1 " || atom[i]->name == " O1P")
        {
          temp[1].push_back(atom[i]);
        }
      else if(atom[i]->name == " O2 " || atom[i]->name == " O2P")
        {
          temp[2].push_back(atom[i]);
        }
      else if(atom[i]->name == " O3 " || atom[i]->name == " O3P")
        {
          temp[3].push_back(atom[i]);
        }
      else
        {
          // This means that this atom is not useful so we flag it 
          atom[i]->skip = true;
        }
    }

  // Error check.  If we don't have at least one of each
  // of the atoms, we throw a warning, set an ignore flag
  // for future reference, and leave the function
  if(temp[0].size() == 0 || temp[1].size() == 0 ||
     temp[2].size() == 0 || temp[3].size() == 0 )
    {
      skip = true;
#ifndef DISABLE_WARNING
      cout << cyan << "WARNING" << reset << ": Could not find all atoms in the " << residue <<  " at "
           << atom[0]->resSeq << endl;
#endif
      return;
    }

  // Allocate size for all of the possible centers
  center.resize( temp[0].size() * temp[1].size() * temp[2].size() * temp[3].size() );

  // And go through all combinations of the alternate locations
  unsigned int index = 0;
  float totalmass = MASS_P + MASS_O * 4;
  for(unsigned int j = 0; j < temp[0].size(); j++) // go through all the P
    {
      for(unsigned int k = 0; k < temp[1].size(); k++) // go through all the O1
        {
          for(unsigned int l = 0; l < temp[2].size(); l++) // go through all the O2
            {
              for(unsigned int m = 0; m < temp[3].size(); m++) // go through all the O3
                {
                  center[index] = temp[0][j]->coord * MASS_P +
                    ( temp[1][k]->coord + temp[2][l]->coord +
                      temp[3][m]->coord  ) * MASS_O;

                  center[index] /= totalmass;

                  center[index].plane_info.resize(4);
                  center[index].plane_info[0] = &(temp[0][j]->coord);
                  center[index].plane_info[1] = &(temp[1][k]->coord);
                  center[index].plane_info[2] = &(temp[2][l]->coord);
                  center[index].plane_info[3] = &(temp[3][m]->coord);

#ifdef DEBUG
                  cout << "CHECK: " << residue  << " " << atom[0]->resSeq << " : " << center[index] << endl;
#endif
                }
            }
        }
    }
}

void AminoAcid::center2POorPO3_charge()
{
  // This makes a vector of a vector of ATOM* 
  // since we are looking for 6 atoms, the size of the outer vector
  // is 3 while the size of the inner vectors are dependent on the 
  // number of alternate locations for each one. Indexes are:
  // P = 0 ; O1 = 1 ; O2 = 2 ; O3 = 3 ; O4 = 4
  Atom* P = NULL;
  vector< Atom* > O;
  vector< Atom* > H;

  // Push all of the important atoms in their respective vectors  
  for(unsigned int i =0; i< atom.size(); i++)
    {
      if(atom[i]->altLoc != ' ')
        {
          altLoc = true;
        }
      if(atom[i]->name == " P  ")
        {
          P = atom[i];
        }
      else if(atom[i]->name == " O1 " || atom[i]->name == " O1P")
        {
          O.push_back(atom[i]);
        }
      else if(atom[i]->name == " O2 " || atom[i]->name == " O2P")
        {
          O.push_back(atom[i]);
        }
      else if(atom[i]->name == " O3 " || atom[i]->name == " O3P")
        {
          O.push_back(atom[i]);
        }
      else if(atom[i]->name == " H  ")
        {
          H.push_back(atom[i]);
        }
      else
        {
          // This means that this atom is not useful so we flag it 
          atom[i]->skip = true;
        }
    }

  // Error check to see if we have 1 P, 4 O, and at least 1 H
  if(P == NULL || O.size() != 3 || H.size() == 0 )
    {
      skip = true;
#ifndef DISABLE_WARNING
      cout << cyan << "WARNING" << reset << ": Could not find all atoms in the " << residue <<  " at "
           << atom[0]->resSeq << endl;
#endif
      return;
    }

  center.resize( 1 );
  // Here is an estimate of the center by doing a charge weighted average
  // This could be wrong, but it is just something to put in here for the
  // time being until I talk with Dr. Hinde.
  center[0] = P->coord * CHARGE_P;
  Coordinates temp(0.0, 0.0, 0.0);
  temp  = O[0]->coord;
  temp += O[1]->coord;
  temp += O[2]->coord;
  center[0] += temp*CHARGE_O;
  temp.set(0.0,0.0,0.0);
  for(unsigned int i=0; i<H.size(); i++)
    {
      temp += H[i]->coord;
    }
  center[0] += temp * CHARGE_H;

  // And since the formal charge of PO4 is -3, divide by it
  if(residue == "2PO")
    center[0] /= -2;
  else if(residue == "PO3")
    center[0] /= -3;
}

// Calculates the center of the amino acid
void AminoAcid::calculateCenter(bool centerOfCharge)
{
  if( !centerOfCharge )
    {
      // TRP
      if ( residue == "TRP" )
        {
          centerTRP();
        }
      // PHE or TYR
      else if ( residue == "PHE" || residue == "TYR" )
        {
          centerPHEorTYR();
        }
      // ASP
      else if (residue == "ASP")
        {
          centerASP_oxygen();
        }
      // GLU
      else if (residue == "GLU")
        {
          centerGLU_oxygen();
        }
      // LYS
      else if (residue == "LYS")
	{
	  centerLYS();
	}
      // ARG
      else if (residue == "ARG")
	{
	  centerARG();
	}
      else if (residue == "PO4" || residue == "2HP" || residue == " PI")
        {
          centerPO4or2HPorPI();
        }
      else if (residue == "2PO" || residue == "PO3")
        {
          center2POorPO3();
        }
      // This is for all of other amino acids out there that we 
      // don't support
      else
        {
          //cerr << red << "ERROR" << reset << ": " << residue << " is not yet supported." << endl;
          skip = true;
        }
    }
  // Now we are going to calculate the center of charges for the formate
  else
    {
      // PHE or TYR
      if ( residue == "PHE" || residue == "TYR" )
        {
          // So, this does a simplified center of mass calculation that
          // Dr. Hinde used.  See function notes above for details.
          centerPHEorTYR_simplified();
        }
      // ASP
      else if (residue == "ASP")
        {
          centerASP_charge();
        }
       // LYS
      else if (residue == "LYS")
	{
	  centerLYS_charge();
	}
       // ARG
      else if (residue == "ARG")
	{
	  centerARG_charge();
	}
       // GLU
      else if (residue == "GLU")
        {
          centerGLU_charge();
        }
      else if (residue == "PO4" || residue == "2HP" || residue == " PI")
        {
          centerPO4or2HPorPI_charge();
        }
      else if (residue == "2PO" || residue == "PO3")
        {
          center2POorPO3();
        }
    }
}


void AminoAcid::calculateAnglesPreHydrogens(AminoAcid aa2,
                                            int index1,
                                            int index2,
                                            float* angle,
                                            float* angle1,
                                            float* angleP)
{
  AminoAcid aa1 = *this;
  Coordinates planeP;
  Coordinates planeProject;

  // Calculate the equation of the plane
  float det = getPlaneEquation( *(aa1.center[index1].plane_info[ CG_PLANE_COORD_PTT]),
                                *(aa1.center[index1].plane_info[CD1_PLANE_COORD_PTT]),
                                *(aa1.center[index1].plane_info[C_2_PLANE_COORD_PTT]),
                                &planeP);

  // Calculate the angle between a plane and a line
  *angle = angleBetweenPlaneAndLine( planeP,
                                     aa1.center[index1],
                                     aa2.center[index2]);
  
  // Calculate the "plane project" coordinates
  planeProjectCoordinate( planeP, 
                          aa2.center[index2],
                          -1 * det, 
                          &planeProject);
  
  // find the second angle betwen the CG ATOM and the "plane project"
  *angle1 = findAngle(*aa1.center[index1].plane_info[ CG_PLANE_COORD_PTT],
                      aa1.center[index1],
                      planeProject);
  
  // finally, calculate the angle between the planes of AA1 and AA2
  *angleP = calculateAngleBetweenPlanes( planeP,
                                         aa2,
                                         index2 );

}

bool AminoAcid::calculateDistancesAndAnglesPostHydrogens(AminoAcid aa2,
                                                         Coordinates closestOxygen,
                                                         float* dist,
                                                         float* distOxy,
                                                         float* distOxy2,
                                                         float* angle,
                                                         float* angleOxy,
                                                         float* angleOxy2)
{
  AminoAcid aa1 = *this;

  // These are the 3 points in the benzene ring determined in centerPHEorTYR_simplified()
  Coordinates dBenzene1 = *aa1.center[0].plane_info[1] - *aa1.center[0].plane_info[0];
  Coordinates dBenzene2 = *aa1.center[0].plane_info[2] - *aa1.center[0].plane_info[0];

  // These values are just to match the Perl script
  float a = dBenzene1.x;
  float b = dBenzene1.y;
  float c = dBenzene1.z;
  float d = dBenzene2.x;
  float e = dBenzene2.y;
  float f = dBenzene2.z;

  // Get the perpendicular vector
  float xp = b * f - c * e;
  float yp = c * d - a * f;
  float zp = a * e - b * d;
  Coordinates perp(xp, yp, zp);  

  // Calculate the distance between the centers
  // This is the vector pointing from the benzene center to formate center of charge
  Coordinates distance = aa2.center[0] - aa1.center[0];
  
  float num = dotProduct(perp, distance);
  float perpnorm = perp.norm();
  float distFromMassToChg = distance.norm();
  float denom = perpnorm * distFromMassToChg;

  if(denom == 0)
    {
      cerr << red << "Error" << reset << ": denom is zero.  Skipping residue" << endl;
      return false;
    }  

  // We already have one of the oxygens as an input param
  // Now let's find the other one
  Coordinates otherOxygen;
  for(int i = 0; i < aa2.atom.size(); i++)
    {
      // Only look at the oxygen atoms
      if(aa2.atom[i]->element == "O")
        {
          // Make sure this oxygen is different than the closest one
          if(aa2.atom[i]->coord.x != closestOxygen.x ||
             aa2.atom[i]->coord.y != closestOxygen.y ||
             aa2.atom[i]->coord.z != closestOxygen.z)
            {
              otherOxygen = aa2.atom[i]->coord;
            }
        }
    }

  // Vector between benzene center of mass and closest oxygen
  Coordinates oD  = closestOxygen - aa1.center[0];
  // Vector between benzene center of mass and the other oxygen
  Coordinates oD2 = otherOxygen - aa1.center[0];
  
  // Oxygen dot product
  float oxy_numerator  = dotProduct(perp, oD); 
  float oxy_numerator2 = dotProduct(perp, oD2); 

  // distance from beneze center and the oxygens
  float distFromCenterToOxy  = oD.norm();
  float distFromCenterToOxy2 = oD2.norm();

  // Denominators
  float oxy_denom  = perpnorm * distFromCenterToOxy;
  float oxy_denom2 = perpnorm * distFromCenterToOxy2;

  if(oxy_denom == 0)
    {
      cerr << red << "Error" << reset << ": oxy_denom are zero.  Skipping residue" << endl;
      return false;
    }
  
  float u = num / denom;
  float uOxy  = oxy_numerator  / oxy_denom;
  float uOxy2 = oxy_numerator2 / oxy_denom2;

  // Force u to be [-1,1]
  if(u > 1)       u =  1;
  else if(u < -1) u = -1;

  if(uOxy > 1)       uOxy =  1;
  else if(uOxy < -1) uOxy = -1;

  if(uOxy2 > 1)       uOxy2 =  1;
  else if(uOxy2 < -1) uOxy2 = -1;

  // Get the angle and change it to degrees
  // take absolute vale to put it in [0,90]
  *dist      = distFromMassToChg;
  *distOxy   = distFromCenterToOxy;
  *distOxy2  = distFromCenterToOxy2;
  *angle     = fabs( 90 - acos(u)     * 180/3.14159 );
  *angleOxy  = fabs( 90 - acos(uOxy)  * 180/3.14159 );
  *angleOxy2 = fabs( 90 - acos(uOxy2) * 180/3.14159 );

}

void AminoAcid::markAltLocAtomsPHEorTYR(int index)
{
  for(int i=0; i < this->atom.size(); i++)
    {
      if( this->atom[i]->name == " CG " &&
          &this->atom[i]->coord != (this->center[index].plane_info[CG_PLANE_COORD_PTT]))
        {
          this->atom[i]->skip = true;
        }
      else if( this->atom[i]->name == " CD1" &&
               &this->atom[i]->coord != (this->center[index].plane_info[CD1_PLANE_COORD_PTT]))
        {
          this->atom[i]->skip = true;
        }
      else if( this->atom[i]->name == " CD2" &&
               &this->atom[i]->coord != (this->center[index].plane_info[CD2_PLANE_COORD_PTT]))
        {
          this->atom[i]->skip = true;
        }
      else if( this->atom[i]->name == " CZ " &&
               &this->atom[i]->coord != (this->center[index].plane_info[3]))
        {
          this->atom[i]->skip = true;
        }
      else if( this->atom[i]->name == " CE1" &&
               &this->atom[i]->coord != (this->center[index].plane_info[4]))
        {
          this->atom[i]->skip = true;
        }
      else if( this->atom[i]->name == " CE2" &&
               &this->atom[i]->coord != (this->center[index].plane_info[5]))
        {
          this->atom[i]->skip = true;
        }
    }  
}

void AminoAcid::markAltLocAtomsLYS(int index)
{
  for (int i = 0; i < this->atom.size(); i++)
    {
      if (this->atom[i]->name == " CE " &&
	  &this->atom[i]->coord != (this->center[index].plane_info[0]))
	{
	  this->atom[i]->skip = true;
	}
// This was, again, going to be used to add hydrogens
// to the side chain without the use of openbabel and
// I had to make sure the atom was not marked as skip.
      else if (this->atom[i]->name == " CD " &&
	       &this->atom[i]->coord != (this->center[index].plane_info[1]))
	{
	  this->atom[i]->skip = true;
	}
      else if (this->atom[i]->name == " NZ " &&
	       &this->atom[i]->coord != (this->center[index].plane_info[2]))
	{
	  this->atom[i]->skip = true;
	}
    }
}

void AminoAcid::markAltLocAtomsARG(int index)
{
  for (int i = 0; i < this->atom.size(); i++)
    {
      if (this->atom[i]->name == " CD " &&
	  &this->atom[i]->coord != (this->center[index].plane_info[0]))
	{
	  this->atom[i]->skip = true;
	}
      else if (this->atom[i]->name == " NE " &&
	  &this->atom[i]->coord != (this->center[index].plane_info[1]))
	{
	  this->atom[i]->skip = true;
	}
      else if (this->atom[i]->name == " CZ " &&
	  &this->atom[i]->coord != (this->center[index].plane_info[2]))
	{
	  this->atom[i]->skip = true;
	}
      else if (this->atom[i]->name == " NH1" &&
	  &this->atom[i]->coord != (this->center[index].plane_info[3]))
	{
	  this->atom[i]->skip = true;
	}
      else if (this->atom[i]->name == "NH2" &&
	  &this->atom[i]->coord != (this->center[index].plane_info[4]))
	{
	  this->atom[i]->skip = true;
	}
    }
}

void AminoAcid::markAltLocAtomsASP(int index)
{
  for(int i=0; i < this->atom.size(); i++)
    {
      if( this->atom[i]->name == " CG " &&
          &this->atom[i]->coord != (this->center[index].plane_info[C__PLANE_COORD_AG]))
        {
          this->atom[i]->skip = true;
        }
      else if( this->atom[i]->name == " OD1" &&
               &this->atom[i]->coord != (this->center[index].plane_info[O_1_PLANE_COORD_AG]))
        {
          this->atom[i]->skip = true;
        }
      else if( this->atom[i]->name == " OD2" &&
               &this->atom[i]->coord != (this->center[index].plane_info[O_2_PLANE_COORD_AG]))
        {
          this->atom[i]->skip = true;
        }
    }  
}

void AminoAcid::markAltLocAtomsGLU(int index)
{
  for(int i=0; i < this->atom.size(); i++)
    {
      if( this->atom[i]->name == " CD " &&
          &this->atom[i]->coord != (this->center[index].plane_info[C__PLANE_COORD_AG]))
        {
          this->atom[i]->skip = true;
        }
      else if( this->atom[i]->name == " OE1" &&
               &this->atom[i]->coord != (this->center[index].plane_info[O_1_PLANE_COORD_AG]))
        {
          this->atom[i]->skip = true;
        }
      else if( this->atom[i]->name == " OE2" &&
               &this->atom[i]->coord != (this->center[index].plane_info[O_2_PLANE_COORD_AG]))
        {
          this->atom[i]->skip = true;
        }
    }  
}

void AminoAcid::markAltLocAtomsPO4or2HPorPI(int index)
{
  for(int i=0; i < this->atom.size(); i++)
    {
      if( this->atom[i]->name == " P  " &&
          &this->atom[i]->coord != (this->center[index].plane_info[0]))
        {
          this->atom[i]->skip = true;
        }
      else if( this->atom[i]->name == " O1 " &&
               &this->atom[i]->coord != (this->center[index].plane_info[1]))
        {
          this->atom[i]->skip = true;
        }
      else if( this->atom[i]->name == " O2 " &&
               &this->atom[i]->coord != (this->center[index].plane_info[2]))
        {
          this->atom[i]->skip = true;
        }
      else if( this->atom[i]->name == " O3 " &&
               &this->atom[i]->coord != (this->center[index].plane_info[3]))
        {
          this->atom[i]->skip = true;
        }
      else if( this->atom[i]->name == " O4 " &&
               &this->atom[i]->coord != (this->center[index].plane_info[4]))
        {
          this->atom[i]->skip = true;
        }
    }  
}

void AminoAcid::markAltLocAtoms2POorPO3(int index)
{
  for(int i=0; i < this->atom.size(); i++)
    {
      if( this->atom[i]->name == " P  " &&
          &this->atom[i]->coord != (this->center[index].plane_info[0]))
        {
          this->atom[i]->skip = true;
        }
      else if( (this->atom[i]->name == " O1 " || this->atom[i]->name == " O1P") &&
               &this->atom[i]->coord != (this->center[index].plane_info[1]))
        {
          this->atom[i]->skip = true;
        }
      else if( (this->atom[i]->name == " O2 "  || this->atom[i]->name == " O2P") &&
               &this->atom[i]->coord != (this->center[index].plane_info[2]))
        {
          this->atom[i]->skip = true;
        }
      else if( (this->atom[i]->name == " O3 "  || this->atom[i]->name == " O2P") &&
               &this->atom[i]->coord != (this->center[index].plane_info[3]))
        {
          this->atom[i]->skip = true;
        }
    }  
}

void AminoAcid::markAltLocAtoms(int index)
{
  if(residue == "PHE" || residue == "TYR")
    {
      markAltLocAtomsPHEorTYR(index);
    }
  else if (residue == "LYS")
    {
      markAltLocAtomsLYS(index);
    }
  else if (residue == "ARG")
    {
      markAltLocAtomsARG(index);
    }
  else if(residue == "ASP")
    {
      markAltLocAtomsASP(index);
    }
  else if(residue == "GLU")
    {
      markAltLocAtomsGLU(index);
    }
  else if(residue == "PO4" || residue == "2HP" || residue == " PI")
    {
      markAltLocAtomsPO4or2HPorPI(index);
    }
  else if(residue == "2PO" || residue == "PO3")
    {
      markAltLocAtoms2POorPO3(index);
    }
}

void AminoAcid::unmarkAltLocAtomsPHEorTYR()
{
  for(int i=0; i < this->atom.size(); i++)
    {
      if( this->atom[i]->name == " CG " )
        {
          this->atom[i]->skip = false;
        }
      else if( this->atom[i]->name == " CD1" )
        {
          this->atom[i]->skip = false;
        }
      else if( this->atom[i]->name == " CD2" )
        {
          this->atom[i]->skip = false;
        }
      else if( this->atom[i]->name == " CZ " )
        {
          this->atom[i]->skip = false;
        }
      else if( this->atom[i]->name == " CE1" )
        {
          this->atom[i]->skip = false;
        }
      else if( this->atom[i]->name == " CE2" )
        {
          this->atom[i]->skip = false;
        }
    }  
}

void AminoAcid::unmarkAltLocAtomsLYS()
{
  for (int i = 0; i < this->atom.size(); i++)
    {
      if (this->atom[i]->name == " CD " )
	{
	  this->atom[i]->skip = false;
	}
      else if (this->atom[i]->name == " CE " )
	{
	  this->atom[i]->skip = false;
	}
      else if (this->atom[i]->name == " NZ " )
	{
	  this->atom[i]->skip = false;
	}
    }
}

void AminoAcid::unmarkAltLocAtomsARG()
{
  for (int i = 0; i < this->atom.size(); i++)
    {
      if (this->atom[i]->name == " CG ")
	{
	  this->atom[i]->skip = false;
	}
      else if (this->atom[i]->name == " CD ")
	{
	  this->atom[i]->skip = false;
	}
      else if (this->atom[i]->name == " NE ")
	{
	  this->atom[i]->skip = false;
	}
      else if (this->atom[i]->name == " CZ ")
	{
	  this->atom[i]->skip = false;
	}
      else if (this->atom[i]->name == " NH1")
	{
	  this->atom[i]->skip = false;
	}
      else if (this->atom[i]->name == " NH2")
	{
	  this->atom[i]->skip = false;
	}
    }
}

void AminoAcid::unmarkAltLocAtomsASP()
{
  for(int i=0; i < this->atom.size(); i++)
    {
      if( this->atom[i]->name == " CG " )
        {
          this->atom[i]->skip = false;
        }
      else if( this->atom[i]->name == " OD1" )
        {
          this->atom[i]->skip = false;
        }
      else if( this->atom[i]->name == " OD2" )
        {
          this->atom[i]->skip = false;
        }
    }  
}

void AminoAcid::unmarkAltLocAtomsGLU()
{
  for(int i=0; i < this->atom.size(); i++)
    {
      if( this->atom[i]->name == " CD " )
        {
          this->atom[i]->skip = false;
        }
      else if( this->atom[i]->name == " OE1" )
        {
          this->atom[i]->skip = false;
        }
      else if( this->atom[i]->name == " OE2" )
        {
          this->atom[i]->skip = false;
        }
    }  
}

void AminoAcid::unmarkAltLocAtomsPO4or2HPorPI()
{
  for(int i=0; i < this->atom.size(); i++)
    {
      if( this->atom[i]->name == " P  " )
        {
          this->atom[i]->skip = false;
        }
      else if( this->atom[i]->name == " O1 " )
        {
          this->atom[i]->skip = false;
        }
      else if( this->atom[i]->name == " O2 " )
        {
          this->atom[i]->skip = false;
        }
      else if( this->atom[i]->name == " O3 " )
        {
          this->atom[i]->skip = false;
        }
      else if( this->atom[i]->name == " O4 " )
        {
          this->atom[i]->skip = false;
        }
    }  
}

void AminoAcid::unmarkAltLocAtoms2POorPO3()
{
  for(int i=0; i < this->atom.size(); i++)
    {
      if( this->atom[i]->name == " P  " )
        {
          this->atom[i]->skip = false;
        }
      else if( this->atom[i]->name == " O1 "  || this->atom[i]->name == " O1P" )
        {
          this->atom[i]->skip = false;
        }
      else if( this->atom[i]->name == " O2 "  || this->atom[i]->name == " O2P" )
        {
          this->atom[i]->skip = false;
        }
      else if( this->atom[i]->name == " O3 "  || this->atom[i]->name == " O3P" )
        {
          this->atom[i]->skip = false;
        }
    }  
}

void AminoAcid::unmarkAltLocAtoms()
{
  if(residue == "PHE" || residue == "TYR")
    {
      unmarkAltLocAtomsPHEorTYR();
    }
  else if (residue == "LYS")
    {
      unmarkAltLocAtomsLYS();
    }
  else if (residue == "ARG")
    {
      unmarkAltLocAtomsARG();
    }
  else if(residue == "ASP")
    {
      unmarkAltLocAtomsASP();
    }
  else if(residue == "GLU")
    {
      unmarkAltLocAtomsGLU();
    }
  else if(residue == "PO4" || residue == "2HP" || residue == " PI")
    {
      unmarkAltLocAtomsPO4or2HPorPI();
    }
  else if(residue == "2PO" || residue == "PO3")
    {
      unmarkAltLocAtoms2POorPO3();
    }
}

string AminoAcid::makeConectPHEorTYR()
{
  string serials[6];
  for(int i=0; i<this->atom.size(); i++)
    {
      if( this->atom[i]->name == " CG " && !this->atom[i]->skip )
        {
          serials[0] = this->atom[i]->line.substr(6,5);
        }
      else if( this->atom[i]->name == " CD1"  && !this->atom[i]->skip )
        {
          serials[1] = this->atom[i]->line.substr(6,5);
        }
      else if( this->atom[i]->name == " CD2"  && !this->atom[i]->skip )
        {
          serials[2] = this->atom[i]->line.substr(6,5);
        }
      else if( this->atom[i]->name == " CZ " && !this->atom[i]->skip )
        {
          serials[3] = this->atom[i]->line.substr(6,5);
        }
      else if( this->atom[i]->name == " CE1" && !this->atom[i]->skip )
        {
          serials[4] = this->atom[i]->line.substr(6,5);
        }
      else if( this->atom[i]->name == " CE2" && !this->atom[i]->skip )
        {
          serials[5] = this->atom[i]->line.substr(6,5);
        }
    }
  string conect = "CONECT" + serials[0] + serials[2] + serials[1] + "                                                 \n";
  conect += "CONECT" + serials[1] + serials[4] + serials[0] + "                                                 \n";
  conect += "CONECT" + serials[2] + serials[5] + serials[0] + "                                                 \n";
  conect += "CONECT" + serials[3] + serials[5] + serials[4] + "                                                 \n";
  conect += "CONECT" + serials[4] + serials[1] + serials[3] + "                                                 \n";
  conect += "CONECT" + serials[5] + serials[2] + serials[3] + "                                                 \n";
  return conect;
}

string AminoAcid::makeConectLYS()
{
  string serials[2];
  for (int i = 0; i < this->atom.size(); i++)
    {
      if (this->atom[i]->name == " NZ " && !this->atom[i]->skip)
	{
	  serials[0] = this->atom[i]->line.substr(6,5);
	}
      else if (this->atom[i]->name == " CE " && !this->atom[i]->skip)
	{
	  serials[1] = this->atom[i]->line.substr(6,5);
	}
    }
  string conect = "CONECT" + serials[0] + serials[1] + "                                                                \n";
  return conect;
}

string AminoAcid::makeConectARG()
{
  string serials[5];
  for (int i = 0; i < this->atom.size(); i++)
    {
      if (this->atom[i]->name == " CD " && !this->atom[i]->skip)
	{
	  serials[0] = this->atom[i]->line.substr(6,5);
	}
      else if (this->atom[i]->name == " NE " && !this->atom[i]->skip)
	{
	  serials[1] = this->atom[i]->line.substr(6,5);
	}
      else if (this->atom[i]->name == " CZ " && !this->atom[i]->skip)
	{
	  serials[2] = this->atom[i]->line.substr(6,5);
	}
      else if (this->atom[i]->name == " NH1" && !this->atom[i]->skip)
	{
	  serials[3] = this->atom[i]->line.substr(6,5);
	}
      else if (this->atom[i]->name == " NH2" && !this->atom[i]->skip)
	{
	  serials[4] = this->atom[i]->line.substr(6,5);
	}
    }
  string conect = "CONECT" + serials[0] + serials[1] + "                                                                \n";
  conect += "CONECT" + serials[1] + serials[0] + serials[2] + "                                                           \n";
  conect += "CONECT" + serials[2] + serials[1] + serials[3] + serials[4] + "                                                      \n";
  conect += "CONECT" + serials[3] + serials[2] + "                                                                \n";
  conect += "CONECT" + serials[4] + serials[2] + "                                                                \n";
  return conect;
}

string AminoAcid::makeConectGLU()
{
  string serials[3];
  for(int i=0; i<this->atom.size(); i++)
    {
      if( this->atom[i]->name == " CD " && !this->atom[i]->skip )
        {
          serials[0] = this->atom[i]->line.substr(6,5);
        }
      else if( this->atom[i]->name == " OE1"  && !this->atom[i]->skip )
        {
          serials[1] = this->atom[i]->line.substr(6,5);
        }
      else if( this->atom[i]->name == " OE2"  && !this->atom[i]->skip )
        {
          serials[2] = this->atom[i]->line.substr(6,5);
        }
    }
  string conect = "CONECT" + serials[0] + serials[1] + serials[2] + "                                                 \n";
  conect += "CONECT" + serials[1] + serials[0] + "                                                       \n";
  conect += "CONECT" + serials[2] + serials[0] + "                                                       \n";
  return conect;
}

string AminoAcid::makeConectASP()
{
  string serials[3];
  for(int i=0; i<this->atom.size(); i++)
    {
      if( this->atom[i]->name == " CG " && !this->atom[i]->skip )
        {
          serials[0] = this->atom[i]->line.substr(6,5);
        }
      else if( this->atom[i]->name == " OD1"  && !this->atom[i]->skip )
        {
          serials[1] = this->atom[i]->line.substr(6,5);
        }
      else if( this->atom[i]->name == " OD2"  && !this->atom[i]->skip )
        {
          serials[2] = this->atom[i]->line.substr(6,5);
        }
    }
  string conect = "CONECT" + serials[0] + serials[1] + serials[2] + "                                                 \n";
  conect += "CONECT" + serials[1] + serials[0] + "                                                       \n";
  conect += "CONECT" + serials[2] + serials[0] + "                                                       \n";
  return conect;
}

string AminoAcid::makeConectPO4or2HPorPI()
{
  string serials[5];
  for(int i=0; i<this->atom.size(); i++)
    {
      if( this->atom[i]->name == " P  " && !this->atom[i]->skip )
        {
          serials[0] = this->atom[i]->line.substr(6,5);
        }
      else if( this->atom[i]->name == " O1 "  && !this->atom[i]->skip )
        {
          serials[1] = this->atom[i]->line.substr(6,5);
        }
      else if( this->atom[i]->name == " O2 "  && !this->atom[i]->skip )
        {
          serials[2] = this->atom[i]->line.substr(6,5);
        }
      else if( this->atom[i]->name == " O3 "  && !this->atom[i]->skip )
        {
          serials[3] = this->atom[i]->line.substr(6,5);
        }
      else if( this->atom[i]->name == " O4 "  && !this->atom[i]->skip )
        {
          serials[4] = this->atom[i]->line.substr(6,5);
        }
    }
  string conect = "CONECT" + serials[0] + serials[1] + serials[2] + serials[3] + serials[4] + "\n";
  conect += "CONECT" + serials[1] + serials[0] + "                                                       \n";
  conect += "CONECT" + serials[2] + serials[0] + "                                                       \n";
  conect += "CONECT" + serials[3] + serials[0] + "                                                       \n";
  conect += "CONECT" + serials[4] + serials[0] + "                                                       \n";
  return conect;
}

string AminoAcid::makeConect2POorPO3()
{
  string serials[5];
  for(int i=0; i<this->atom.size(); i++)
    {
      if( (this->atom[i]->name == " P  ") && !this->atom[i]->skip )
        {
          serials[0] = this->atom[i]->line.substr(6,5);
        }
      else if( (this->atom[i]->name == " O1 " || this->atom[i]->name == " O1P" )  && !this->atom[i]->skip )
        {
          serials[1] = this->atom[i]->line.substr(6,5);
        }
      else if( (this->atom[i]->name == " O2 " || this->atom[i]->name == " O1P" )  && !this->atom[i]->skip )
        {
          serials[2] = this->atom[i]->line.substr(6,5);
        }
      else if( (this->atom[i]->name == " O3 " || this->atom[i]->name == " O1P" )  && !this->atom[i]->skip )
        {
          serials[3] = this->atom[i]->line.substr(6,5);
        }
    }
  string conect = "CONECT" + serials[0] + serials[1] + serials[2] + serials[3] + "      \n";
  conect += "CONECT" + serials[1] + serials[0] + "                                                       \n";
  conect += "CONECT" + serials[2] + serials[0] + "                                                       \n";
  conect += "CONECT" + serials[3] + serials[0] + "                                                       \n";
  return conect;
}
string AminoAcid::makeConect()
{
  if(residue == "PHE" || residue == "TYR")
    {
      return makeConectPHEorTYR();
    }
  else if (residue == "LYS")
    {
      return makeConectLYS();
    }
  else if (residue == "ARG")
    {
      return makeConectARG();
    }
  else if(residue == "ASP")
    {
      return makeConectASP();
    }
  else if(residue == "GLU")
    {
      return makeConectGLU();
    }
  else if(residue == "PO4" || residue == "2HP" || residue == " PI")
    {
      return makeConectPO4or2HPorPI();
    }
  else if(residue == "2PO" || residue == "PO3")
    {
      return makeConect2POorPO3();
    }
}

// This checks the validity of the hydrogens for 
// GLU and ASP residues.  Sometimes babel will add 4 hydrogens:
// 2 to C and 1 to each O.  We can figure out why, but we came
// up with this hackish fix.  We take the 2 that are connected
// to the carbon, average them, and use that as the coordinates 
// for the hydrogen that we are looking for.  We then throw all
// of the other hydrogens away.
bool AminoAcid::removeExcessHydrogens(vector<string> conect)
{
  if( !(residue == "GLU" || residue == "ASP") )
    return false;
  
  vector<Atom*>::iterator it;
  int carbonSerialNumber;
  int hydrogenCount = 0;
  for(it = atom.begin(); it != atom.end(); ++it)
    {
      if( (*it)->name == " CG " || (*it)->name == " CD " )
        {
          carbonSerialNumber = (*it)->serialNumber;
        }
      else if( (*it)->name == " H  " )
        {
          hydrogenCount++;
        }
    }
  corrected = false;
  if( hydrogenCount > 1 )
    {
      corrected = true;
      Coordinates avg(0,0,0);
      vector<string>::iterator conectit;
      Atom* lastHydrogen;
      for(it = atom.begin(); it < atom.end(); ++it)
        {
          for(conectit = conect.begin(); conectit < conect.end(); ++conectit)
            {
              int number;
              int cnumber;
              if(!from_string<int>(number,(*conectit).substr(7,5),dec))
                {
                  cerr << red << "Error" << reset << ": failed to convert serial number into an int"  << endl;
                }
              string tempstr = (*conectit).substr(11,5);
              if(tempstr != "     ")
                {
                  if(!from_string<int>(cnumber,(*conectit).substr(11,5),dec))
                    {
                      cerr << red << "Error" << reset << ": failed to convert serial number into an int "  << endl;
                    }

                  if( (*it)->name == " H  " && (*it)->serialNumber == number 
                      && cnumber == carbonSerialNumber )
                    {
                      avg += (*it)->coord;
                    }
                }
            }
          if((*it)->name == " H  " )
            {
              lastHydrogen = (*it);
              vector<Atom*>::iterator tempit = it-1;
              atom.erase(it);
              it = tempit;
            }
        }
      avg /= 2;
      lastHydrogen->coord = avg;
      char cstr[25];
      sprintf(cstr,"%8.3lf%8.3lf%8.3lf",avg.x, avg.y, avg.z);
      string temp(cstr);
      lastHydrogen->line = lastHydrogen->line.substr(0,30) + 
        temp.substr(0,24) + lastHydrogen->line.substr(54,26);
      atom.push_back(lastHydrogen);
    }
  return corrected;
}

void AminoAcid::printPHEorTYR(FILE* output)
{
  for(int i=0; i < this->atom.size(); i++)
    {
      if((atom[i]->name == " CD1" || 
          atom[i]->name == " CD2" || 
          atom[i]->name == " CE1" || 
          atom[i]->name == " CE2" || 
          atom[i]->name == " CZ " ||
          atom[i]->name == " CG ") && !atom[i]->skip)
        {
          this->atom[i]->print(output);
        }
    }
}

void AminoAcid::printLYS(FILE* output)
{
  for(int i=0; i < this->atom.size(); i++)
    {
      if((atom[i]->name == " CE " || 
          atom[i]->name == " NZ ") && !atom[i]->skip)
        {
          this->atom[i]->print(output);
        }
    }
}

void AminoAcid::printARG(FILE* output)
{
  for(int i=0; i < this->atom.size(); i++)
    {
      if((atom[i]->name == " CG " || 
          atom[i]->name == " CD " || 
          atom[i]->name == " NE " || 
          atom[i]->name == " CZ " || 
          atom[i]->name == " NH1" ||
          atom[i]->name == " NH2") && !atom[i]->skip)
        {
          this->atom[i]->print(output);
        }
    }
}

void AminoAcid::printASP(FILE* output)
{
  for(int i=0; i < this->atom.size(); i++)
    {
      if((atom[i]->name == " OD1" || 
         atom[i]->name == " OD2" || 
         atom[i]->name == " CG " ) && !atom[i]->skip)
        {
          this->atom[i]->print(output);
        }
    }
}

void AminoAcid::printGLU(FILE* output)
{
  for(int i=0; i < this->atom.size(); i++)
    {
      if((atom[i]->name == " OE1" || 
         atom[i]->name == " OE2" || 
         atom[i]->name == " CD " )&& !atom[i]->skip)
        {
          this->atom[i]->print(output);
        }
    }
}

void AminoAcid::printNeededAtoms(FILE* output)
{

  if(residue == "PHE" || residue == "TYR")
    {
      printPHEorTYR(output);
    }
  else if (residue == "LYS")
    {
      printLYS(output);
    }
  else if (residue == "ARG")
    {
      printARG(output);
    }
  else if(residue == "ASP")
    {
      printASP(output);
    }
  else if(residue == "GLU")
    {
      printGLU(output);
    }
  else
    {
      cerr << red << "Error" << reset << ": Unsupported residue" << endl;
    }
}

ostream& operator<<(ostream& output, const AminoAcid& p) 
{
  for(int i=0; i < p.atom.size(); i++)
    {
      output << *(p.atom[i]) << endl;
    }
  return output;  // for multiple << operators.
}
