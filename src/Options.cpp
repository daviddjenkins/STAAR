/****************************************************************************************************/
//  COPYRIGHT 2011, University of Tennessee
//  Author: David Jenkins (david.d.jenkins@gmail.com)
//  File: Options.cpp
//  Date: 10 Jan 2011
//  Version: 1.0
//  Description: Implementation file for the command line parser and Options class
//
//  TODO: after fixing center of charge stuff, change the default back to false
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

#include "Options.hpp"

// Initialize the Options to empty stuff
Options::Options()
{
  center = true;
  pdbfile = NULL;
  outputfile = NULL;
  failure = false;
  threshold = 15.0;
}

// Intialize options then parse the cmd line arguments
Options::Options( int argc, char **argv )
{
  center = true;
  pdbfile = NULL;
  outputfile = NULL;
  failure = false;
  threshold = 15.0;
  parseCmdline( argc, argv );
}

// Empty destructor
Options::~Options(){}

// Print a help page
void printHelp()
{
  cerr << "-h or --help          " << "Displays this message" << endl;
  cerr << "-p or --pdbdir        " << "Specifies the folder for PDB files" << endl;
  cerr << "-o or --out           " << "Specifies the output file" << endl;
  cerr << "-c or --center        " << "Specifies whether to calculate center of charge" << endl;
  cerr << "-t or --threshold     " << "Set the distance threshold between amino acids" << endl;
}

// Return true of cmd line parsing failed, false otherwise
bool Options::fail()
{
  return failure;
}

// Parses the cmd line options
void Options::parseCmdline( int argc, char **argv )
{
  
  int c;
  // Set the options that are possible with this program
  static struct option long_options[] = 
    {
      {"help",		required_argument, 0, 'h'},
      {"pdbdir",	required_argument, 0, 'p'},
      {"out",		required_argument, 0, 'o'}, 
      {"center",	required_argument, 0, 'c'},
      {"threshold",	required_argument, 0, 't'},
      {0, 0, 0, 0}
    };
  int option_index;

  // Go through the options and set them to variables
  while( !( ( c = getopt_long(argc, argv, "hp:o:ct:", long_options, &option_index) ) < 0 ) )
    {
    switch(c)
      {
      case 'h':
	// the user is asking for help here
        printHelp();
        break;

      case 'p':
	// the user is specifying the PDB file/directory
        this->pdbfile = optarg;
        break;

      case 'o':
	// the user is wanting to set the output file, but we need to make
	// sure that is ia a file and not a directory
	if( !isDirectory( optarg ) )
	  {
	    this->outputfile = optarg;
	  }
	else
	  {
	    cerr << "Output file must be a file, not a directory" << endl;
	    printHelp();
	    exit(1);
	  }
        break;

      case 'c':
	// user wants to do the center of charge calculations, but really
	// there is no option at this point because we need to fix the 
	// setting of the plane variables
        this->center = true;
        break;

      case 't':
	// user wants to set the threshold value but we need to make
	// sure that an actual number was inputted and not just a string
	if ( sscanf(optarg, "%f", &(this->threshold)) != 1 )
	  {
	    cerr << "Must insert a valid number for the threshold" << endl;
	    printHelp();
	    exit(1);
	  }
        break;

      default:
	printHelp();
	exit(1);
	break;
      }
  }

  // Error check to make sure there wasn't extra crap on the command line
  // and all necessary arguments were set
  if( optind < argc ){
    cerr << "There are some extra commands inputted. Please cleanup the command line!" << endl;
    // print some help info
    printHelp();
    failure=true;
  }else if( !pdbfile ){
    cerr << "Must specify the PDB list file with -p or --pdblist" <<  endl;
    printHelp();
    failure=true;
  }else if( !outputfile ){
    cerr << "Must specify the op file with -o or --op" <<  endl;
    printHelp();
    failure=true;
  }

}


