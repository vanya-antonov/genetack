/***************************************************************************************************
*                                                                                                  *
* Author: Ivan Antonov (ivan.antonov@gatech.edu)                                                   *
*                                                                                                  *
***************************************************************************************************/

// $Id$

#include "path.h"

/***************************************************************************************************
*                                                                                                  *
* Path() -- constructor of the Path class                                                          *
*                                                                                                  *
* Arguments:                                                                                       *
*     int num_steps  -- number of steps in the path. In our case this is the sequence length.      *
*     int num_states -- number of stated in each step. In our case this is the # of states in HMM. *
*                                                                                                  *
* Note:                                                                                            *
*     In fact these two arguments are dimentions of all_steps table: num_steps x num_states        *
*                                                                                                  *
***************************************************************************************************/
Path::Path(const int num_steps, const int num_states)
{
	int i;

	this->num_steps    = num_steps;
	this->first_step_i = 0;
	this->num_states   = num_states;
	this->states_seq   = new int[num_steps];
	this->emiss_seq    = new int[num_steps];
	this->all_steps    = new Path_step*[ num_steps ];
	for(i = 0; i < this->num_steps; ++i)
		this->all_steps[i] = new Path_step[num_states];
}

/***************************************************************************************************
*                                                                                                  *
* Path() -- copy constructor of the Path class. We have to define it because we have pointers      *
*            in this class.                                                                        *
*                                                                                                  *
***************************************************************************************************/
Path::Path(const Path &old_path)
{
	int i, j;
	
	this->num_steps    = old_path.num_steps;
	this->first_step_i = old_path.first_step_i;
	this->num_states   = old_path.num_states;
	this->states_seq   = new int[ this->num_steps ];
	this->all_steps    = new Path_step*[ this->num_steps ];
	for(i = 0; i < this->num_steps; ++i)
		this->all_steps[i] = new Path_step[num_states];
	
	for(i = 0; i < this->num_steps; ++i)
	{
		this->states_seq[i] = old_path.states_seq[i];
		for(j = 0; j < this->num_states; ++j)
			this->all_steps[i][j] = old_path.all_steps[i][j];
	}
}

/***************************************************************************************************
*                                                                                                  *
* ~Path() -- destructor of the Path class                                                          *
*                                                                                                  *
***************************************************************************************************/
Path::~Path()
{
	int i;
	
	delete this->states_seq;
	delete this->emiss_seq;
	
	for(i = 0; i < this->num_steps; i++)
		delete this->all_steps[i];
	delete this->all_steps;
}
