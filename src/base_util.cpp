/***************************************************************************************************
*                                                                                                  *
* Author: Ivan Antonov (antonov.ivan@gatech.edu)                                                   *
*                                                                                                  *
***************************************************************************************************/

// $Id$

#include "base_util.h"

/***************************************************************************************************
*                                                                                                  *
* uc_string() -- modifies the given string translating it in uppercase                             *
*                                                                                                  *
* Arguments:                                                                                       *
*     str -- string to translate in uppercase                                                      *
*                                                                                                  *
* Returns:                                                                                         *
*     string str -- initial string translated in uppercase                                         *
*                                                                                                  *
***************************************************************************************************/
string uc_string(string &str)
{
	for(int i=0; i < (int)str.length(); ++i)  // for each character in string
		str[i] = toupper(str[i]);   // convert char to uppercase
	return str;
}

int trim_begin(string &str, const string symbols)
{
	int  str_i, symb_i, len;
	bool to_trim;

	len = 0; // Length of the chunk to be erased
	for(str_i = 0; str_i < (int)str.length(); str_i++)
	{
		to_trim = false;
		for(symb_i = 0; symb_i < (int)symbols.length(); symb_i++)
		{
			if( str[str_i] == symbols[symb_i] )
			{
				to_trim = true;
				break;
			}
		}
		
		if(to_trim)
			len++;
		else
			break;
	}
	
	if(len > 0)
		str.erase(str.begin(), str.begin() + len);

	return len;
}

int trim_end(string &str, const string symbols)
{
	int  str_i, symb_i, len;
	bool to_trim;

	len = 0; // Length of the chunk to be erased
	for(str_i = (int)str.length()-1; str_i >= 0; str_i--)
	{
		to_trim = false;
		for(symb_i = 0; symb_i < (int)symbols.length(); symb_i++)
		{
			if( str[str_i] == symbols[symb_i] )
			{
				to_trim = true;
				break;
			}
		}
		
		if(to_trim)
			len++;
		else
			break;
	}
	
	if(len > 0)
		str.erase(str.end() - len, str.end());

	return len;
}

/***************************************************************************************************
*                                                                                                  *
* trip() -- removes all specified symbols from the beginning and / or end of the string            *
*                                                                                                  *
* Arguments:                                                                                       *
*     0 -- trim both beginning and end of the string                                               *
*     1 -- trim beginning of the string only                                                       *
*     2 -- trim end of the string only                                                             *
*                                                                                                  *
* Returns:                                                                                         *
*     int num -- total number of character has been removed                                        *
*                                                                                                  *
***************************************************************************************************/
int trim(string &str, const string symbols, const int mode)
{
	int num = 0;

	if(mode == 0 || mode == 1)   // Trim the beginning of the string
		num += trim_begin(str, symbols);

	if(mode == 0 || mode == 2)   // Trim the end of the string
		num += trim_end(str, symbols);
	
	return num;
}
