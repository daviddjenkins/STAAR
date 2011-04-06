/****************************************************************************************************/
//  COPYRIGHT 2011, University of Tennessee
//  Author: David Jenkins (david.d.jenkins@gmail.com)
//  File: Coordinates.hpp
//  Date: 11 Jan 2011
//  Version: 1.0
//  Description: Contains the class definition and implementation for Coordinates
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


#ifndef __COORDINATES_HPP__
#define __COORDINATES_HPP__

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cmath>

using namespace std;

class Coordinates{
public:
  // Obviously, the x, y, and z coordinates
  float x;
  float y;
  float z;

  // The following is a string that contains the
  // alternate location sequence
  string altLoc;

  vector<Coordinates*> plane_info;

  // Simple constructor that sets everything to 0
  Coordinates();
  
  // Constructor that will set the values
  Coordinates(float x1,
              float y1,
              float z1);

  // Sets the coordinate values
  void set(float x1,
           float y1,
           float z1);

  // Computes the distance between this and another point
  float distance( Coordinates& point );

  // Computes the length/norm of a vector
  float norm();

  // Empty destructor
  ~Coordinates();

  // Operator overloads
  Coordinates& operator=(const Coordinates &rhs);
  bool operator==(const Coordinates &rhs) const;
  bool operator!=(const Coordinates &rhs) const;
  Coordinates & operator+=(const Coordinates &rhs);
  Coordinates & operator-=(const Coordinates &rhs);
  Coordinates & operator*=(const Coordinates &rhs);
  const Coordinates operator+(const Coordinates &other) const;
  const Coordinates operator-(const Coordinates &other) const;
  const Coordinates operator*(const Coordinates &other) const;
  Coordinates& operator*=(const float &num);
  const Coordinates operator*(const float &num) const;
  Coordinates& operator/=(const float &num);
  const Coordinates operator/(const float &num) const;

  friend ostream& operator<<(ostream& output, const Coordinates& p);

};


#endif
