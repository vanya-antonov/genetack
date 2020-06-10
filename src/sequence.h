/***************************************************************************************************
*                                                                                                  *
* Definitions for the Sequence class. Sequence class is designed to work with biological           *
* sequencies of nucleotides or amino acids                                                         *
*                                                                                                  *
* Author: Ivan Antonov (ivan.antonov@gatech.edu)                                                   *
*                                                                                                  *
***************************************************************************************************/

// $Id$

#ifndef _SEQUENCE_H_INCLUDED
#define _SEQUENCE_H_INCLUDED

#include <string>
#include <map>
#include <iostream>
#include <fstream>
#include <stdlib.h>

#include "base_util.h"

using namespace std;

class Sequence
{
	private:
		string next_seqname;
		map    <string, int> all_seqnames;
	
	public:
		string seqname;
		string seq;
		
		Sequence(const string fn);
		bool read_next_seq( istream* in );
		int length();
};

#endif /* _SEQUENCE_H_INCLUDED */

