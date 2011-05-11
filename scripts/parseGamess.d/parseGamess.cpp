/****************************************************************************************************/
//  COPYRIGHT 2011, University of Tennessee
//  Author: David Jenkins (david.d.jenkins@gmail.com)
//  File: parseGamess.cpp
//  Date: 29 Mar 2011
//  Version: 1.0
//  Description: Parses through the GAMESS output and appends the MK results onto the STAAR
//               output table
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

#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <iostream>
#include <sstream>
#include <fstream>
#include <string>
#include <vector>

using namespace std;

// Separates a string by whitespace
vector<string> get_fields(const string &s)
{
  vector<string> fields;
  string field;
  istringstream iss(s);
  while (iss >> field)
    {
      fields.push_back(field);
    }
  return fields;
}

// Splits a string by a given character delimiter
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

// Converts a string to any class specified by T
template <class T>
bool from_string(T&                     t,
                 const std::string&     s,
                 std::ios_base& (*f)(std::ios_base&))
{
  std::istringstream    iss(s);
  return !(iss >> f >> t).fail();
}


int main(int argc, char* argv[]){
  if(argc != 3){
    fprintf(stderr,"%s STAAR_table GAMESS_outdir > new_STAAR_table_output\n",argv[0]);
    exit(1);
  }

  ifstream staar_table(argv[1]);
  string line;
  string gout_dir(argv[2]);

  // Read in the STAAR output line by line
  while( getline(staar_table, line) )
    {
      // If this is a comment line, we can skip it
      if(line[0] == '#')
        {
          cout << line << endl;
          continue;
        }

      // Since the table is a CSV, separate it up by commas
      vector<string> fields = split(line,',');

      // This is the filepath of the gamess input file for this pair
      string gamess_inp = fields[10];

      // this points to the name of the input file
      size_t pos = fields[10].rfind("/");
      string gname;

      // This will get the filename of the gamess output file
      if( pos == string::npos )
        {
          pos = 0;
          while((pos = fields[10].find("inp",pos+1)) != string::npos)
            {
              fields[10].replace(pos, 3, "out");
            }
          gname = gout_dir + "/" + fields[10];
        }
      else
        {
          int pos2 = 0;
          while((pos2 = fields[10].find("inp",pos2+1)) != string::npos)
            {
              fields[10].replace(pos2, 3, "out");
            }

          gname = gout_dir + "/" + fields[10].substr(pos + 1);
        }

      // open up the GAMESS output file
      ifstream gamess_out( gname.c_str() );
      if(!gamess_out)
        {
          perror(gname.c_str());
          exit(1);
        }

      // Append the filename of the GAMESS output file
      line += "," + gname;
      string line2;
      int count=0;

      // Go through each line of the file to find the results we want
      while( getline(gamess_out, line2) )
        {
          if( line2.find("ELECTROSTATIC ENERGY")       != string::npos ||
              line2.find("POLARIZATION ENERGY")        != string::npos ||
              line2.find("EXCHANGE REPULSION ENERGY")  != string::npos ||
              line2.find("CHARGE TRANSFER ENERGY ")    != string::npos ||
              line2.find("TOTAL INTERACTION ENERGY")   != string::npos ||
              line2.find("HIGH ORDER COUPLING ENERGY") != string::npos   )
            {
              string har  = get_fields(line2.substr(37,12))[0];
              string kcal = get_fields(line2.substr(52,8))[0];
              if( kcal == "********")
                {
                  count = -2;
                  break;
                }
              // And append it to the end of the line
              line += "," + har 
                + "," + kcal;
              count ++;
            }
          // Cut out results that have a high mix term
          if( line2.find("HIGH ORDER COUPLING ENERGY") != string::npos )
            {
              float mix;
              if( !from_string<float>(mix, get_fields(line2.substr(52,8))[0], dec) )
                {
                  count = -2;
                  break;
                }
              if( fabs(mix) > 0.25 )
                {
                  count = -1;
                  break;
                }
            }
        }

      // output the new line or print error if the results weren't present
      if(count == -1)
        cerr << "Skipped because mix energy was too high : "<< gamess_inp << " : " 
             << fields[0] << fields[6] << " - " << fields[1] << fields[7] << " pair in " << fields[9] << endl;
      else if(count == -2)
        cerr << "Skipping because the run did not converge within 200 iterations: " << gamess_inp << endl;
      else if(count != 0)
        cout << line << endl;
      else 
        cerr << "No results found for " << gamess_inp << " : " 
             << fields[0] << fields[6] << " - " << fields[1] << fields[7] << " pair in " << fields[9] << endl;

      // Close the gamess output file
      gamess_out.close();
    }

  // close the STAAR results file
  staar_table.close();

  return 0;
}
