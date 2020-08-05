/***************************************************************************************************
*                                                                                                  *
* fsmark -- program to search frame shift sequencing errors in procaryotic genes by HMM model.     *
*           FSMark uses GeneMark to preprocess given DNA in order to find genes in it.             *
*                                                                                                  *
* Author: Ivan Antonov (ivan.antonov@gatech.edu)                                                   *
*                                                                                                  *
* TODO:                                                                                            *
*     * Improve output format:                                                                     *
*           - implement Params::print()                                                            *
*     * Bug: read_fasta - doesn't read last file line if the line doesn't have '\n' at the end     *
*     * Use pointers in function arguments because otherwise we always copy the data               *
*     * Add ability to read sequence in GeneBank format                                            *
*     * Is it possible NOT to specify id's for states in the hmm_def file?                         *
*     * Add HMM check: each state has to have unique name and this name shouldn't be -1 (end state)*
*                                                                                                  *
****************************************************************************************************/

// $Id$

// To reduce warnings in MS VS 2005
#define _CRT_SECURE_NO_DEPRECATE

#include <iostream>
#include <iomanip>
#include <fstream>
#include <stdio.h>
#include <math.h>
#include <vector>
#include <ctime>

#include "sequence.h"
#include "hmm.h"
#include "params.h"
#include "path.h"

using namespace std;

void print_fs_array(string seqname, int gene_start, int gene_end, vector <int> fs_arr, Path &path)
{
	int i;

	for(i=0; i < fs_arr.size(); ++i)
	{
		int fs_pos     = fs_arr[i];
		int from_state = path.states_seq[fs_pos-1];
		int to_state   = path.states_seq[fs_pos];
		if( from_state == to_state )
		{
			cerr << "Frameshift position " << fs_pos << " is wrong!!!" << endl;
			continue;
		}
		
		string fs_type = "-1";
		if((from_state==0 && to_state==1) || (from_state==1 && to_state==2) || (from_state==2 && to_state==0))
			fs_type = "+1";
		
		int gene_len = gene_end - gene_start + 1;
		cout << left << setw(13)<< seqname    << "\t" << setw(8) << fs_pos   << "\t" << setw(7) << fs_type  << "\t"
		             << setw(10)<< gene_start << "\t" << setw(8) << gene_end << "\t" << setw(8) << gene_len << endl;
	}
}

void print_header(Params &prm)
{
	// Get date & time. Taken from: http://www.dreamincode.net/code/snippet1102.htm
	time_t curtime = time(0); 
	tm now=*localtime(&curtime); 
	char date[256]={0};
	strftime(date, sizeof(date)-1, "%A, %B %d %Y %X", &now);
	
	cout << "###" << endl
	     << "# "+prm.display_name+" Version "+prm.version+", revision " << prm.svn_revision << " ("+prm.svn_date+")" << endl
	     << "# " << prm.command_line << endl
	     << "# " << date << endl
	     << "#" << endl
	     << "sequence_name\tfs_coord\tfs_type\tgene_start\tgene_end\tgene_len" << endl;
}
	
void output(Params &prm, Sequence &seq, Hmm &hmm, Path &path)
{
	int pos_i, gene_start, overlap_start, prev_state_i, cur_state_i;
	bool inside_gene;
	vector <int> fs_arr;
	
	// We start with max_order because there is no states before 
	inside_gene   = false;
	overlap_start = -1;
	prev_state_i  = -1;
	
	for(pos_i=hmm.get_max_order(); pos_i < path.num_steps; ++pos_i)
	{
		cur_state_i = path.states_seq[pos_i];
		
		// Set prev_state_i if we can
		if( pos_i > hmm.get_max_order() )
			prev_state_i = path.states_seq[pos_i-1];
		
		// Check if it is frameshift
		if(0 <= cur_state_i && cur_state_i <= 2 && 0 <= prev_state_i && prev_state_i <= 2 && cur_state_i!=prev_state_i)
			fs_arr.push_back(pos_i);
		
		if(!inside_gene && (hmm.is_start_state(cur_state_i) || hmm.is_coding_state(cur_state_i)) )
		{
			inside_gene = true;
			gene_start  = pos_i==hmm.get_max_order() ? 0 : (pos_i-2);
		}
		else if(inside_gene && (hmm.is_end_state(cur_state_i) || pos_i == path.num_steps-1) )
		{
			if( fs_arr.size() > 0 )
				print_fs_array(seq.seqname,gene_start,pos_i,fs_arr,path);
			fs_arr.clear();
			
			if(overlap_start > 0)
			{
				// We are inside gene ocerlap
				gene_start    = overlap_start;
				overlap_start = -1;
			}
			else
			{
				inside_gene = false;
			}
		}
		else if(inside_gene && hmm.is_start_state(cur_state_i))
		{
			overlap_start = pos_i;   // Overlap started
		}
	}
}

void output_table(Sequence &seq, Hmm &hmm, Path &path, const string out_fn)
{
	int pos_i, state_i;
	string name;
	
	ofstream out_table;
	out_table.open(out_fn.c_str());
	
	out_table << fixed;   // Not scientific in form 1e+6
	
	out_table << "Position\tLetter\tState\tEmission\tBest_score";
	for(state_i = 0; state_i < hmm.get_num_states(); ++state_i)
	{
		name = hmm.get_state_name( state_i );
		out_table << "\tState " << name << " sum\tState " << name << " ptr";
	}
	out_table << endl;
	
	// Very beginning of sequence -- where we don't have enough letters behind yet
	for(pos_i = 0; pos_i < hmm.get_max_order(); ++pos_i)
		out_table << pos_i << "\t" << seq.seq[pos_i] << endl;
	
	for(; pos_i < path.num_steps; ++pos_i)
	{
		int best_state_i = path.states_seq[pos_i];
		out_table << pos_i << "\t" << seq.seq[pos_i] << "\t" << hmm.get_state_name( best_state_i )
		          << "\t" << hmm.get_emiss_name( path.emiss_seq[pos_i] ) << "\t" << path.all_steps[pos_i][best_state_i].sum;
		for(state_i = 0; state_i < hmm.get_num_states(); ++state_i)
			out_table << "\t" << setprecision(3) << path.all_steps[pos_i][state_i].sum << "\t" << path.all_steps[pos_i][state_i].ptr;
		out_table << endl;
	}
	
	out_table.close();
}

int main(int argc, char* argv[])
{
	Params prm(argc, argv);
	
	Hmm hmm(prm.hmm_fn);
	if( prm.seq_starts_with_ATG )
		hmm.max_seq_offset = 2;
	
	ifstream infile;
	infile.open(prm.seq_fn.c_str(), ios::in);
	if(!infile.is_open())
	{
		cerr << "Unable to open file" << prm.seq_fn << "!\n";
		exit(1);
	}
	
	print_header(prm);
	
	Sequence seq(prm.seq_fn);
	while( seq.read_next_seq( &infile ) )
	{
		Path path(seq.length(), hmm.get_num_states());
		hmm.run_viterbi(seq, path);
		
		if( prm.table_fn.length() > 0 )
			output_table(seq, hmm, path, prm.table_fn+seq.seqname);
		
		output(prm, seq, hmm, path);
	}
	
	return(0);
}

