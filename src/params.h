/***************************************************************************************************
*                                                                                                  *
* Definitions for the Params class.                                                                *
* Params class is designed to handle command line parameters for FSMark program.                   *
* This header file contains also version # and version date of the program.                        *
*                                                                                                  *
* Author: Ivan Antonov (ivan.antonov@gatech.edu)                                                   *
*                                                                                                  *
***************************************************************************************************/

// $Id$

#ifndef _PARAMS_H_INCLUDED
#define _PARAMS_H_INCLUDED

#include <string>
#include <iostream>
#include <stdlib.h>
#include <ctime>

using namespace std;

class Params
{
	private:
		string prg_name;      // Program name (executable fiel)
		string descr;         // One line version description
		
		void init();
		void parse_svn_str(const string svn_str);
	public:
		string command_line;
		string svn_date;
		int    svn_revision;
		string version;           // Program version #
		string display_name;      // Program name to display
		string seq_fn;
		string hmm_fn;
		string table_fn;
		bool   seq_starts_with_ATG;
		
		void print();
		void usage();

		Params(const int argc, char* argv[]);
};

#endif /* _PARAMS_H_INCLUDED */

