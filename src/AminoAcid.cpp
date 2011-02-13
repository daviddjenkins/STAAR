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

AminoAcid::AminoAcid()
{
  atom.clear();
  center.clear();
  plane_info.clear();
  skip = false;
  altLoc = false;
}

AminoAcid::~AminoAcid()
{
  atom.clear();
  center.clear();
  plane_info.clear();
  skip = false;
  altLoc = false;
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
    }
  
  // Error check.  If we don't have at least one of each
  // of the atoms, we throw a warning, set an ignore flag
  // for future reference, and leave the function
  if(temp[0].size() == 0 || temp[1].size() == 0 ||
     temp[2].size() == 0 || temp[3].size() == 0 ||
     temp[4].size() == 0 || temp[5].size() == 0 )
    {
      skip = true;
      cout << "WARNING: Could not find all atoms in the TRP ring at "
	   << atom[0]->resSeq << endl;
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
			  center[index].plane_info.resize(3);
			  center[index].plane_info[ CG_PLANE_COORD_PTT] = &(temp[0][i]->coord);
			  center[index].plane_info[CD1_PLANE_COORD_PTT] = &(temp[2][k]->coord);
			  center[index].plane_info[CD2_PLANE_COORD_PTT] = &(temp[3][l]->coord);
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
			  cout << "CHECK: TRP" << atom[0]->resSeq << " : " << center[index] << endl;
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
      cout << "WARNING: Could not find all atoms in the TRP ring at "
	   << atom[0]->resSeq << endl;
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
      else if(atom[i]->name == " CB ")
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
    }

  // Error check.  If we don't have at least one of each
  // of the atoms, we throw a warning, set an ignore flag
  // for future reference, and leave the function
  if(temp[0].size() == 0 || temp[1].size() == 0 ||
     temp[2].size() == 0 || temp[3].size() == 0)
    {
      skip = true;
      cout << "WARNING: Could not find all atoms in the ASP side chain at "
	   << atom[0]->resSeq << endl;
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
      if(atom[i]->name == " CG ")
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
    }

  // Error check.  If we don't have at least one of each
  // of the atoms, we throw a warning, set an ignore flag
  // for future reference, and leave the function
  if(temp[0].size() == 0 || temp[1].size() == 0 ||
     temp[2].size() == 0 || temp[3].size() == 0)
    {
      skip = true;
      cout << "WARNING: Could not find all atoms in the GLU side chain at "
	   << atom[0]->resSeq << endl;
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
// OD1 and OD2 coords to the center charge.
// This also sets the plane coordinates used to calculate 
// the angle later in the program
void AminoAcid::centerASP_nocharge()
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
    }

  // Error check.  If we don't have at least one of each
  // of the atoms, we throw a warning, set an ignore flag
  // for future reference, and leave the function
  if( temp[1].size() == 0 ||
      temp[2].size() == 0)
    {
      skip = true;
      cout << "WARNING: Could not find all atoms in the ASP side chain at "
	   << atom[0]->resSeq << endl;
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

// Calculate the centers for ASP
// these AA need the following atoms:
//   OE1 and OE2
// all of the other ones don't affect the center
// This is is really simple: all it does is add the
// OE1 and OE2 coords to the center charge.
// This also sets the plane coordinates used to calculate 
// the angle later in the program
void AminoAcid::centerGLU_nocharge()
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
    }

  // Error check.  If we don't have at least one of each
  // of the atoms, we throw a warning, set an ignore flag
  // for future reference, and leave the function
  if(temp[0].size() == 0 || temp[1].size() == 0 ||
     temp[2].size() == 0)
    {
      skip = true;
      cout << "WARNING: Could not find all atoms in the GLU side chain at "
	   << atom[0]->resSeq << endl;
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

// Calculates the center of the amino acid
void AminoAcid::calculateCenter(bool center)
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
      // Center of charge
      if( center )
	{
	  centerASP();
	}
      else
	{
	  centerASP_nocharge();
	}
    }
  // GLU
  else if (residue == "GLU")
    {
      // Center of charge
      if( center )
	{
	  centerGLU();
	}
      else
	{
	  centerGLU_nocharge();
	}
    }
  // This is for all of other amino acids out there that we 
  // don't support
  else
    {
      //cerr << "ERROR: " << residue << " is not yet supported." << endl;
    }
}

ostream& operator<<(ostream& output, const AminoAcid& p) {
  for(int i=0; i < p.atom.size(); i++)
    {
      output << *(p.atom[i]) << endl;
    }
  return output;  // for multiple << operators.
}
