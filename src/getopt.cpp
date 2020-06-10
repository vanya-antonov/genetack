/***************************************************************************************************
*                                                                                                  *
* Author: Aliksandr Kravchenko (alexk@genebee.msu.su)                                              *
*                                                                                                  *
***************************************************************************************************/

// $Id$

#include <stdio.h>
#include <string.h>
#include <stdlib.h>     // atoi()   , exit()
#include "getopt.h"

GetOpt::GetOpt(int _argc, char* const *_argv, const char* _optstring)
{
	optind = 0;
	optarg = NULL;
	next   = NULL;

	argc = _argc;
	argv = _argv;
	optstring = _optstring;

	if (argc < 0 || !argv || !optstring) argc = 0;
}

int GetOpt::next_opt () 
{
	if (optind == 0) next = NULL;
	optarg = NULL;
	if (next == NULL || *next == '\0')
	{
		if (optind == 0) optind++;
		if (optind >= argc) return EOF; 

		if (argv[optind][0] != '-') // Isn't key
		{
			// I don't know what Kravchenko meant here, so I just die. You can take a look at his code through SVN.
			cerr << "Someting wrong inside GetOpt: optind " << optind << ", argv[optind] = '" << argv[optind] << "'\n";
			exit(1);
		}
		
		if (argv[optind][1] == '\0')
		{
			if (optind < argc) optarg = argv[optind];
			else optarg = NULL;
			return EOF;
		}
		if (strcmp(argv[optind], "--") == 0)
		{
			optind++;
			if (optind < argc) optarg = argv[optind];
			else optarg = NULL;
			return EOF;
		}
		next = argv[optind];
		next++;		// skip past -
		optind++;
	}
	
	char  c = *next++;
	char* cp = (char*)strchr(optstring, c);
	if (cp == NULL || c == ':') return '?';

	if (cp[1] == ':') // Need argument for option c
	{
		if (*next != '\0') {
			optarg = next;
			next = NULL;
		} else if (optind < argc) {
			optarg = argv[optind];
			optind++;
		} else
			return '?';
	}
	return c;
}
