/***************************************************************************************************
*                                                                                                  *
* Definitions for the Path class. Path class is designed to store pathes HMM goes for a            *
* particular sequence in Viterbi algorithm.                                                        *
*                                                                                                  *
* Author: Ivan Antonov (ivan.antonov@gatech.edu)                                                   *
*                                                                                                  *
***************************************************************************************************/

// $Id$

#ifndef _PATH_H_INCLUDED
#define _PATH_H_INCLUDED

#include "sequence.h"

struct Path_step
{
	double sum;
	int    ptr;   // PoinTeR to the previous state
};

class Path
{
	private:
	public:
		int         num_steps;                // Sequence length
		int         first_step_i;             // Index of the first step (sometimes it may not be 0)
		int         num_states;               // Number of states in each step
		int*        states_seq;               // Best sequence of HMM states determined by Viterbi, for instance
		int*        emiss_seq;                // Sequence of emission ids that corresponds to states_seq
		Path_step** all_steps;                // Represents two tables (ptr and sum) of equal size ( NUM_STEPS x NUM_STATES ).
		                                      // all_steps[row][col]   <=>   all_steps[seq_i][state_i]
		Path(const Path &old_path);           // Copy constructor
		Path(const int num_steps, const int num_states);
		~Path();
};

#endif /* _PATH_H_INCLUDED */
