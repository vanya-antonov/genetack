/***************************************************************************************************
*                                                                                                  *
* Author: Ivan Antonov (ivan.antonov@gatech.edu)                                                   *
*                                                                                                  *
***************************************************************************************************/

// $Id$

#include <stdio.h>      // printf() , fopen()  , fclose()
#include <stdlib.h>     // atoi()   , exit()

#include "params.h"
#include "getopt.h"

/***************************************************************************************************
*                                                                                                  *
* Params() -- constructor of the Params class                                                      *
*                                                                                                  *
* Arguments:                                                                                       *
*     int   argc   -- number of command line arguments for the program                             *
*     char* argv[] -- array of command line arguments for the program                              *
*                                                                                                  *
***************************************************************************************************/
Params::Params(const int argc, char* argv[])
{
	this->init();
	
	this->command_line = "";
	for(int i = 0; i < argc; ++i)
	{
		if(i > 0)
			this->command_line += " ";
		this->command_line.append( argv[i] );
	}
	
	GetOpt getopt(argc, argv, "f:m:t:S");
	
	int c;
	while((c = getopt.next_opt()) != EOF) 
	{
		switch (c)
      	{  
         	case 'f':
				this->seq_fn.assign( getopt.optarg );
				break;
         	case 'm':
				this->hmm_fn.assign( getopt.optarg );
				break;
         	case 't':
				this->table_fn.assign( getopt.optarg );
				break;
			case 'S':
				this->seq_starts_with_ATG = true;
				break;
		}
	}
	
	if(this->seq_fn.length() == 0 || this->hmm_fn.length() == 0 )
	{
		this->usage();
		exit(0);
	}
}

void Params::init()
{
	this->display_name        = "GeneTack";
	this->prg_name            = "genetack";
	this->version             = "0.53";
	this->descr               = "";
	this->seq_starts_with_ATG = false;
	this->parse_svn_str("$Id$");
}

void Params::parse_svn_str(const string svn_str)
{
	int start, end;
	
	// Derive svn_revision
	start = (int)svn_str.find(" ", 5) + 1;
	end   = (int)svn_str.find(" ", start);
	this->svn_revision = atoi( svn_str.substr(start, end-start).c_str() );
	
	// Derive svn_date
	start = end+1;
	end   = (int)svn_str.find(" ", start);
	this->svn_date = svn_str.substr(start, end-start);
}

/***************************************************************************************************
*                                                                                                  *
* usage() -- prints usage for FSMark program to stderr                                             *
*                                                                                                  *
***************************************************************************************************/
void Params::usage()
{
	cerr << "\n"+display_name+", version "+version+", revision " << svn_revision << " ("+svn_date+")\n";
	if( !this->descr.empty() )
		cerr << "\n"+this->descr+"\n";
	cerr << "\n";
	cerr << "OPTIONS:\n";
	cerr << "    -t PREFIX  --  save Viterbi tables to files named PREFIX+SEQUENCE_NAME.\n";
	cerr << "                   PREFIX can be directory name, for example 'dir/'.\n";
	cerr << "    -S         --  the input sequence may start with start codon (e.g. ATG)\n";
	cerr << "                   (by default 5' UTR is expected).\n";
	cerr << "\n";
	cerr << "USAGE:\n";
	cerr << "    "+prg_name+" [OPTIONS] -m MODEL.hmm_def -f SEQUENCE.fasta > OUT.txt\n";
	cerr << "\n";
}

/***************************************************************************************************
*                                                                                                  *
* print() -- prints parameters list to stdout. It's useful to put in the output file list of       *
*            parameters with which program had been executed.                                      *
*                                                                                                  *
***************************************************************************************************/
void Params::print()
{
}

