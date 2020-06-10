/***************************************************************************************************
*                                                                                                  *
* Author: Ivan Antonov (ivan.antonov@gatech.edu)                                                   *
*                                                                                                  *
***************************************************************************************************/

// $Id$

#include "sequence.h"

/***************************************************************************************************
*                                                                                                  *
* Sequence() -- constructor of the Sequence class                                                  *
*                                                                                                  *
* Arguments:                                                                                       *
*     string seq_fn -- sequence file name. If seq_fn == "" then stdin is used.                     *
*                      Following formats are currently supported:                                  *
*                          * fasta                                                                 *
*                                                                                                  *
***************************************************************************************************/
Sequence::Sequence(const string fn)
{
}

/***************************************************************************************************
*                                                                                                  *
* length() -- returns length of sequence. Is used because there is a warning when call just        *
*             seq.seq.length(). You should always add (int) at the beginning of this statement.    *
*                                                                                                  *
***************************************************************************************************/
int Sequence::length()
{
	return (int)this->seq.length();
}

/***************************************************************************************************
*                                                                                                  *
* read_fasta() -- reads sequence from fasta file and assigns varialbles of the class               *
*                                                                                                  *
* Returns:                                                                                         *
*     bool -- Has the sequence been successfully read                                              *
*                                                                                                  *
***************************************************************************************************/
bool Sequence::read_next_seq( istream* in )
{
	string tmp;
	
	if( in->eof() )
		return 0;
	
	getline(*in, tmp, '\n');
	trim(tmp, " \t\n\r");
	if( tmp.substr(0,1).compare(">") == 0 )
	{
		// Very first sequence of the file
		this->seqname = tmp.substr(1);
		this->seq     = "";
	}
	else if( this->next_seqname.length() > 0 )
	{
		this->seqname = this->next_seqname;
		this->seq     = tmp;
	}
	else
	{
		cerr << "Something wrong -- Can't determine sequence name" << endl;
		exit(1);
	}
	
	while( !getline(*in, tmp, '\n').eof() )
	{
		trim(tmp, " \t\n\r");
		if( tmp.substr(0,1).compare(">") == 0 )
		{
			this->next_seqname = tmp.substr(1);
			break;
		}
		else
		{
			this->seq += tmp;
		}
	}

	if( tmp.substr(0,1).compare(">") != 0 && tmp.length() > 0 )
		this->seq += tmp;  // Tail of the last sequence
	
	// Check if seqnames are unique
	if( all_seqnames.find(this->seqname) == all_seqnames.end() )
		this->all_seqnames[ this->seqname ] = 1;
	else
	{
		cerr << "ERROR: sequence names must be unique: '" << this->seqname << "' is duplicated!\n";
		exit(1);
	}
	
	return 1;
}

