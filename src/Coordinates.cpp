/****************************************************************************************************/
//  COPYRIGHT 2011, University of Tennessee
//  Author: David Jenkins (david.d.jenkins@gmail.com)
//  File: Coordinates.cpp
//  Date: 11 Jan 2011
//  Version: 1.0
//  Description: Contains the class definition for Coordinates
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

#include "Coordinates.hpp"

// Initialize coordinates to 0
Coordinates::Coordinates()
{
  x=0.0;
  y=0.0;
  z=0.0;    
}

// Set the coordinates to something else during construction
Coordinates::Coordinates(float x1,
                         float y1,
                         float z1)
{
  x=x1;
  y=y1;
  z=z1;
}

// Set the coordinates to something else
void Coordinates::set(float x1,
                      float y1,
                      float z1)
{
  this->x=x1;
  this->y=y1;
  this->z=z1;
}

float Coordinates::distance(Coordinates& point)
{
  Coordinates t = (*this-point);
  t *= t;
  return sqrt(t.x + t.y + t.z);
}

float Coordinates::norm()
{
  return sqrt(x * x + y * y + z * z);
}

// Empty destructor
Coordinates::~Coordinates()
{
  plane_info.clear();
}

// Overload the = operator
Coordinates& Coordinates::operator=(const Coordinates &rhs) 
{
  if (this != &rhs) 
    {
      this->x = rhs.x;
      this->y = rhs.y;
      this->z = rhs.z;
      return *this; 
    }
  return *this;
}

// Overload the += operator
Coordinates & Coordinates::operator+=(const Coordinates &rhs) 
{
  this->x += rhs.x;
  this->y += rhs.y;
  this->z += rhs.z;
  return *this;
}

// Overload the -= operator
Coordinates & Coordinates::operator-=(const Coordinates &rhs) 
{
  this->x -= rhs.x;
  this->y -= rhs.y;
  this->z -= rhs.z;
  return *this;
}

// Overload the *= operator
Coordinates & Coordinates::operator*=(const Coordinates &rhs) 
{
  this->x *= rhs.x;
  this->y *= rhs.y;
  this->z *= rhs.z;
  return *this;
}

// Overload the + operator
const Coordinates Coordinates::operator+(const Coordinates &other) const 
{
  return Coordinates(*this) += other;
}

// Overload the - operator
const Coordinates Coordinates::operator-(const Coordinates &other) const 
{
  return Coordinates(*this) -= other;
}

// Overload the * operator
const Coordinates Coordinates::operator*(const Coordinates &other) const 
{
  return Coordinates(*this) *= other;
}

// Overload the Coordinates *= float operatory
Coordinates& Coordinates::operator*=(const float &num)
{
  this->x *= num;
  this->y *= num;
  this->z *= num;
  return *this;
}

// Overload the Coordinates * float operatory
const Coordinates Coordinates::operator*(const float &num) const 
{
  return Coordinates(*this) *= num;
}

// Overload the Coordinates /= float operatory
Coordinates& Coordinates::operator/=(const float &num)
{
  this->x /= num;
  this->y /= num;
  this->z /= num;
  return *this;
}

// Overload the Coordinates / float operatory
const Coordinates Coordinates::operator/(const float &num) const 
{
  return Coordinates(*this) /= num;
}


ostream& operator<<(ostream& output, const Coordinates& p) {
  output <<  p.x << "\t" << p.y << "\t" << p.z;
  return output;  // for multiple << operators.
}
