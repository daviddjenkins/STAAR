/****************************************************************************************************/
//  COPYRIGHT 2011, University of Tennessee
//  Author: David Jenkins (david.d.jenkins@gmail.com)
//  File: Geometry.hpp
//  Date: 20 Jan 2011
//  Version: 1.0
//  Description: Function declarations for a number of analytical geometry values
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


#ifndef __GEOMETRY_HPP__
#define __GEOMETRY_HPP__

#include <cstdlib>
#include <cstdio>
#include <iostream>
#include <fstream>
#include <cmath>

#include "Coordinates.hpp"
#include "AminoAcid.hpp"

// Find the dot product between 2 points
float dotProduct(Coordinates& point1,
                 Coordinates& point2);

// Finds the determinant of 3 points with their
// indices are column entries of a row. So:
//      x1 y1 z1
//      x2 y2 z2
//      x3 y3 z3
float determinant(Coordinates& point1,
                  Coordinates& point2,
                  Coordinates& point3);

// Gets the coordinates for the plain
float getPlaneEquation(Coordinates& point1,
		       Coordinates& point2,
		       Coordinates& point3,
		       Coordinates* result);

// Calculate the angle between a plane and a line
float angleBetweenPlaneAndLine(Coordinates& plane,
                               Coordinates& point1,
                               Coordinates& point2);

// calculates the plane-line intercept constant
float constantForPlaneLineIntercept(Coordinates& plane,
                                    Coordinates& point,
                                    float intercept);

// I guess projects a plane onto a coordinate
// taken straight from the orignal STAAR code
// the answer is stored in result
void planeProjectCoordinate(Coordinates& plane,
                            Coordinates& point,
                            float intercept,
                            Coordinates* result);

// Given 2 points and the coordinates where they intersect,
// this function finds the angle between them
float findAngle(Coordinates& point1,
                Coordinates& intercept,
                Coordinates& point2);

// this calculates the angle between two planes.  In this case,
// we take as arguments a precalculated set of coordinates for 
// one plane (to reduce the number of computations), the second
// residue in question, and the index for the center that we
// are curious about
float calculateAngleBetweenPlanes( Coordinates& planeP,
				   AminoAcid& aa2,
				   int index2 );

// Converts the vector to a unit vector
// USE ONLY ON VECTORS, NOT COORDINATES!
// Return Coordinates or void?
Coordinates unitVector(Coordinates& point);

// Computes a cross product of point1 x point2
// USE ONLY ON VECTORS, NOT COORDINATES!
Coordinates crossProduct(Coordinates& point1,
						 Coordinates& point2);


// Define a vector from point1 to point2
Coordinates defineVector(Coordinates& point1,
						 Coordinates& point2);

#endif
