/****************************************************************************************************/
//  COPYRIGHT 2011, University of Tennessee
//  Author: David Jenkins (david.d.jenkins@gmail.com)
//  File: CoutColors.cpp
//  Date: 09 Mar 2011
//  Version: 1.0
//  Description: Makes outputting color easier
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

#include "CoutColors.hpp"

using namespace std;

static int _COLOR_FORMAT_FLAG;

ostream& reset(ostream &os)
{
#ifndef DISABLE_COUTCOLORS
  if(_COLOR_FORMAT_FLAG==1)
    {
      os << "\e[m";
      _COLOR_FORMAT_FLAG = 0;
    }
#endif
  return os;
}

ostream& black(ostream &os)
{
#ifndef DISABLE_COUTCOLORS
  if(_COLOR_FORMAT_FLAG==1)
    os << "\e[m";

  _COLOR_FORMAT_FLAG=1;
  os << "\e[30m";
#endif
  return os;
}

ostream& red(ostream &os)
{
#ifndef DISABLE_COUTCOLORS
  if(_COLOR_FORMAT_FLAG==1)
    os << "\e[m";

  _COLOR_FORMAT_FLAG=1;
  os << "\e[31m";
#endif
  return os;
}

ostream& green(ostream &os)
{
#ifndef DISABLE_COUTCOLORS
  if(_COLOR_FORMAT_FLAG==1)
    os << "\e[m";

  _COLOR_FORMAT_FLAG=1;
  os << "\e[32m";
#endif
  return os;
}

ostream& brown(ostream &os)
{
#ifndef DISABLE_COUTCOLORS
  if(_COLOR_FORMAT_FLAG==1)
    os << "\e[m";

  _COLOR_FORMAT_FLAG=1;
  os << "\e[33m";
#endif
  return os;
}

ostream& blue(ostream &os)
{
#ifndef DISABLE_COUTCOLORS
  if(_COLOR_FORMAT_FLAG==1)
    os << "\e[m";

  _COLOR_FORMAT_FLAG=1;
  os << "\e[34m";
#endif
  return os;
}

ostream& purple(ostream &os)
{
#ifndef DISABLE_COUTCOLORS
  if(_COLOR_FORMAT_FLAG==1)
    os << "\e[m";

  _COLOR_FORMAT_FLAG=1;
  os << "\e[35m";
#endif
  return os;
}

ostream& cyan(ostream &os)
{
#ifndef DISABLE_COUTCOLORS
  if(_COLOR_FORMAT_FLAG==1)
    os << "\e[m";

  _COLOR_FORMAT_FLAG=1;
  os << "\e[36m";
#endif
  return os;
}

ostream& gray(ostream &os)
{
#ifndef DISABLE_COUTCOLORS
  if(_COLOR_FORMAT_FLAG==1)
    os << "\e[m";

  _COLOR_FORMAT_FLAG=1;
  os << "\e[37m";
#endif
  return os;
}

