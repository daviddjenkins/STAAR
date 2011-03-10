/****************************************************************************************************/
//  COPYRIGHT 2011, University of Tennessee
//  Author: David Jenkins (david.d.jenkins@gmail.com)
//  File: Utils.cpp
//  Date: 11 Jan 2011
//  Version: 1.0
//  Description: Contains a number of functions that will be used to do various, non application
//               specific tasks such as converting a string to any type
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

#include "Utils.hpp"
#include "CoutColors.hpp"


// Simply makes a system call to count the number of files a directory
int countFilesInDirectory( char* path )
{
  char command[strlen(path)+12];
  int count;

  // Write the command out
  sprintf(command, "ls %s | wc -l", path);

  // Send the system command while piping the output into this program
  FILE* fp = popen( command, "r" );

  // Try to read the command output
  if ( !fgets(command, strlen(path), fp) ){
    cerr << red << "Error" << reset << ": Failed to get the number of files in given directory" << endl;
    return -1;
  }

  // Try to convert the output to an int
  if(!from_string<int>(count,command,dec)){
    cerr << red << "Error" << reset << ": Failed to convert number into integer to get number of files in given directory" << endl;
    return -1;
  }
  fclose(fp);
  return count;
}

// Returns true if the path given is a directory,
// false if it is a file, and exits if either
// it doesn't exists or is some other type of path
bool isDirectory( char* path )
{
  int status;
  struct stat st_buf;

  status = stat (path, &st_buf);
  if (status != 0) {
    return false;
  }

  // Just a regular ol' file...
  if (S_ISREG (st_buf.st_mode)) 
    {
      return false;
    }
  // this is a directory...
  else  if (S_ISDIR (st_buf.st_mode)) 
    {
      return true;
    }
  // and this is something else!
  else 
    {
      cerr << red << "Error" << reset << ": Inputted file must be either a single file or a directory." << endl;
      exit(1);
    }
}

// returns the current time in seconds.
// used for timing the program
double getTime()
{
  struct timeval t;
  gettimeofday(&t,NULL);
  return t.tv_sec+((double)t.tv_usec)/1000000.0;
}

vector<string> split(const string &s, char delim) 
{
  vector<string> elems;
  stringstream ss(s);
  string item;
  while(getline(ss, item, delim))
    {
      elems.push_back(item);
    }
  return elems;
}

void printHeader(ostream & os)
{
  // This is just the STAAR logo
  os << brown << " #####    #######      #         #      ######  " << endl;
  os << brown << "#     #      #        # #       # #     #     # " << endl;
  os << brown << "#            #       #   #     #   #    #     # " << endl;
  os << brown << " #####       #      #     #   #     #   ######  " << endl;
  os << brown << "      #      #      #######   #######   #   #   " << endl;
  os << brown << "#     #      #      #     #   #     #   #    #  " << endl;
  os << brown << " #####       #      #     #   #     #   #     # " << endl;
  os << brown << " STAAR: \e[33mST\e[matistical \e[33mA\e[mnalysis of "
     <<"\e[33mA\e[mromatic \e[33mR\e[mings " << endl;
  os << endl;
  // This is just some other basic info
  os << "Version 2.0" << endl;
  os << "University of Tennessee (C) 2011" << endl;
  os << endl;
}

