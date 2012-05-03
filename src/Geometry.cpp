/****************************************************************************************************/
//  COPYRIGHT 2011, University of Tennessee
//  Author: David Jenkins (david.d.jenkins@gmail.com)
//  File: Geometry.cpp
//  Date: 20 Jan 2011
//  Version: 1.0
//  Description: Contains the function implementations for the declarations in 
//               Geometry.hpp.  These functions to things such as calculate
//               angles, evaluate planes, calculate dot product, etc etc
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

#include "float.h"
#include "Geometry.hpp"
#include "CoutColors.hpp"

// Find the dot product between 2 points
float dotProduct( Coordinates& point1,
                  Coordinates& point2 )
{
  Coordinates t = point1 * point2;
  return t.x + t.y + t.z;
}

// Finds the determinant of 3 points with their
// indices are column entries of a row. So:
//      x1 y1 z1
//      x2 y2 z2
//      x3 y3 z3
float determinant( Coordinates& point1,
                   Coordinates& point2,
                   Coordinates& point3 )
{
  float t11, t12, t13;
  float t21, t22, t23;
  float t31, t32, t33;

  // This is only written like to make sure it is clear
  // what each column and row consists of
  t11 = point1.x;    t12 = point1.y;    t13 = point1.z;
  t21 = point2.x;    t22 = point2.y;    t23 = point2.z;
  t31 = point3.x;    t32 = point3.y;    t33 = point3.z;
  
  // Standard 3x3 determinant equation
  // (unless I made a typo)
  return ( t11 * ( t22*t33 - t23*t32 ) - 
           t12 * ( t21*t33 - t23*t31 ) +
           t13 * ( t21*t32 - t22*t31 ));
}

// Gets the coordinates for the plain
float getPlaneEquation( Coordinates& point1,
                        Coordinates& point2,
                        Coordinates& point3,
                        Coordinates* result )
{
  Coordinates t1, t2, t3;
  t1.set(1, point1.y, point1.z);
  t2.set(1, point2.y, point2.z);
  t3.set(1, point3.y, point3.z);

  // Get the determinant of a matrix with the 
  // x values set to 1
  result->x = determinant( t1, 
                           t2, 
                           t3 );

  t1.set(point1.x, 1, point1.z);
  t2.set(point2.x, 1, point2.z);
  t3.set(point3.x, 1, point3.z);

  // Get the determinant of a matrix with the 
  // y values set to 1
  result->y = determinant( t1, 
                           t2, 
                           t3 );

  t1.set(point1.x, point1.y, 1);
  t2.set(point2.x, point2.y, 1);
  t3.set(point3.x, point3.y, 1);

  // Get the determinant of a matrix with the 
  // z values set to 1
  result->z = determinant( t1, 
                           t2, 
                           t3 );

  // and because we will need it later, return the
  // determinant of the 3 points
  return determinant( point1, 
                      point2, 
                      point3 );
}

// Calculate the angle between a plane and a line
float angleBetweenPlaneAndLine ( Coordinates& plane,
                                 Coordinates& point1,
                                 Coordinates& point2 )
{
  float normalv;
  float normalp;
  float dotProd;
  float angle;
  float cosValue;

  Coordinates tempCoord;
  
  // Calculate the norm of the plane coordinates
  normalv = plane.norm();

  // subtract the 2 points together and get their norm
  tempCoord = point2 - point1;
  normalp = tempCoord.norm();

  // get the dot producted between the plane and difference of the points
  dotProd = dotProduct(plane, tempCoord);

  // Now we can get the cosine value with this information
  cosValue = dotProd/(normalv * normalp);

  // double check that we aren't doing anything illegal with arccos
  if( cosValue < -1.0 || cosValue > 1.0 )
    {
      cerr << red << "Error" << reset << ":Error: Cosine of angle, " << cosValue << ", lies outside of -1 and 1" << endl;
      angle = 1000;
    }
  else
    {
      // and actually calculate the angle
      angle = abs ( 90 - ( acos(cosValue) * 180 / 3.14159  ) );
    }
  return angle;
}

// calculates the plane-line intercept constant
float constantForPlaneLineIntercept(Coordinates& plane,
                                    Coordinates& point,
                                    float intercept)
{

  //denom = plane.x*plane.x + plane.y*plane.y + plane.z*plane.z;
  float denom = dotProduct( plane, 
                            plane );
  
  // Make sure we aren't dividing by 0
  // because that baaaaaaaaaaaaad
  if(denom != 0)
    {
       return -1*( intercept + 
                   dotProduct(plane, point) ) / denom;
    }
  else
    {
      cout << "Warning: Plane vector does not exist" << endl;
    }
  return FLT_MAX;
}

// I guess projects a plane onto a coordinate
// taken straight from the orignal STAAR code
// the answer is stored in result
void planeProjectCoordinate(Coordinates& plane,
                            Coordinates& point,
                            float intercept,
                            Coordinates* result)
{
  float t = constantForPlaneLineIntercept(plane, point, intercept);
  *result = point + plane * t;
}

// Given 2 points and the coordinates where they intersect,
// this function finds the angle between them
float findAngle(Coordinates& point1,
                Coordinates& intersect,
                Coordinates& point2)
{
  Coordinates vec1;
  Coordinates vec2;

  float normVec1;
  float normVec2;
  float dotProd;
  float cosValue;

  // make vectors from these points
  vec1 = point1 - intersect;
  vec2 = point2 - intersect;

  // and find their norms
  normVec1 = vec1.norm();
  normVec2 = vec2.norm();

  // now let's find the dot product between these 2 vectors
  dotProd = dotProduct(vec1, vec2);
  
  // and calculate the cosine of this angle
  cosValue = dotProd / ( normVec1 * normVec2 );

  // error check to make sure we don't do something retarded
  if( cosValue < -1 || cosValue >1 )
    {
      cerr << red << "Error" << reset << ":Error: Cosine of angle, " << cosValue << ", lies outside of -1 and 1" << endl;
      return 1000;      
    }
  else
    {
      // and return the angle
      return acos(cosValue) * 180 / 3.14159;
    }

}

// this calculates the angle between two planes.  In this case,
// we take as arguments a precalculated set of coordinates for 
// one plane (to reduce the number of computations), the second
// residue in question, and the index for the center that we
// are curious about
float calculateAngleBetweenPlanes( Coordinates& planeP,
                                   AminoAcid& aa2,
                                   int index2 )
{
  Coordinates planeQ;
  
  float normV;
  float normP;
  float dotProd;
  float cosValue;

  // get the plane equation for the center we are curious about
  getPlaneEquation( *(aa2.center[index2].plane_info[ C__PLANE_COORD_AG]),
                    *(aa2.center[index2].plane_info[O_1_PLANE_COORD_AG]),
                    *(aa2.center[index2].plane_info[O_2_PLANE_COORD_AG]),
                    &planeQ);

  // calculate the norms of these planes
  normV = planeP.norm();
  normP = planeQ.norm();

  // and their dot product
  dotProd = dotProduct( planeP,
                        planeQ );

  // and use them to calculate the value of the cosine of the angle
  cosValue = dotProd / ( normV * normP );

  // error check
  if( cosValue < -1 || cosValue > 1 )
    {
      cerr << red << "Error" << reset << ":Error: Cosine of angle, " << cosValue << ", lies outside of -1 and 1" << endl;
      return 1000;            
    }
  else
    {
      // and actually calculate the angle value
      float angle = acos(cosValue) * 180 / 3.14159;
      if(angle > 90)
        {
          angle = 180.0 - angle;
        }
      return angle;
    }
}

// Converts the vector to a unit vector
// USE ONLY ON VECTORS, NOT COORDINATES!
// Return Coordinates or void?
Coordinates unitVector( Coordinates& point )
{
	if (abs(point.norm() - 1) < .00001)
	{
		return point;
	}
	float norm = point.norm();
	Coordinates newPoint((point.x / norm),
						 (point.y / norm),
						 (point.y / norm));
	return newPoint;
}

// Computes a cross product of point1 x point2
// USE ONLY ON VECTORS, NOT COORDINATES!
Coordinates crossProduct(Coordinates& point1,
						 Coordinates& point2)
{
	Coordinates newPoint((point1.y * point2.z - point1.z * point2.y),
						 (point1.z * point2.x - point1.x * point2.z),
						 (point1.x * point2.y - point1.y * point2.x));
	return newPoint;
}

// Define a vector from point1 to point2
Coordinates defineVector(Coordinates& point1,
						 Coordinates& point2)
{
	return point2 - point1;
}
