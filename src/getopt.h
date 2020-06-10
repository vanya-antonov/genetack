/***************************************************************************************************
*                                                                                                  *
* Author: Aliksandr Kravchenko (alexk@genebee.msu.su)                                              *
*                                                                                                  *
***************************************************************************************************/

// $Id$

#include <iostream>

using namespace std;

class GetOpt
{
	int	optind;
	char*	next;

	int argc;
	char* const *argv;
	const char* optstring;
public:
	GetOpt(int argc, char* const *argv, const char* optstring);
	
	int next_opt();

	const char* optarg;	// key value 
};
