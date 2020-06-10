/***************************************************************************************************
*                                                                                                  *
* Definitions for functions of basic purpose                                                       *
*                                                                                                  *
* Author: Ivan Antonov (antonov.ivan@gatech.edu)                                                   *
*                                                                                                  *
***************************************************************************************************/

// $Id$

#include <string>

using namespace std;

#ifndef _BASE_UTIL_H_INCLUDED
#define _BASE_UTIL_H_INCLUDED

string uc_string(string &str);
int trim(string &str, const string symbols, const int mode = 0);

#endif /* _BASE_UTIL_H_INCLUDED */
