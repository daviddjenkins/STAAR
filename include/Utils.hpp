/****************************************************************************************************/
//  COPYRIGHT 2011, University of Tennessee
//  Author: David Jenkins (david.d.jenkins@gmail.com)
//  File: Utils.hpp
//  Date: 11 Jan 2011
//  Version: 1.0
//  Description: Contains a number of function declarations that will be used to do various,
//               non-application specific tasks such as converting a string to any type
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

#ifndef __UTILS_HPP__
#define __UTILS_HPP__

#include <cstdlib>
#include <cstdio>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cmath>
#include <getopt.h>
#include <unistd.h>
#include <dirent.h>
#include <vector>
#include <sys/time.h>
#include <sys/stat.h>

using namespace std;

// Simply makes a system call to count the number of files a directory
int countFilesInDirectory( char* path );

// checks to see if the path given is a directory
bool isDirectory( char* path );

// returns the time in seconds
double getTime();

// Splits a line up by a given delimiter.  Used many for parsing the
// residue command line option
vector<string> split(const string &s, char delim);

//from_string template also included with Utils.hpp
// Converts a string to any class specified by T
template <class T>
bool from_string(T&                     t,
                 const std::string&     s,
                 std::ios_base& (*f)(std::ios_base&))
{
  std::istringstream    iss(s);
  return !(iss >> f >> t).fail();
}


#endif
