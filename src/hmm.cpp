/***************************************************************************************************
*                                                                                                  *
* Author: Ivan Antonov (ivan.antonov@gatech.edu)                                                   *
*                                                                                                  *
***************************************************************************************************/

// $Id$

#include "hmm.h"

/***************************************************************************************************
*                                                                                                  *
* Hmm() -- constructor of the Hmm class. Reads the hmm_def file with HMM model definition.         *
*                                                                                                  *
* Arguments:                                                                                       *
*     const char* fn -- file name with HMM definition                                              *
*                                                                                                  *
***************************************************************************************************/
Hmm::Hmm(const string fn)
{
	this->init(fn);

	ifstream in( fn.c_str() );
	
	if( in.fail() )
	{
		cerr << "Can't open file '" + fn + "'!";
		exit(1);
	}
	
	this->load_header_info(in);
	this->load_options(in);
	this->load_emission_list(in);
	this->load_state_list(in);
	this->load_transition_list(in);
	in.close();
	this->check_hmm();
	
	this->finalize();
}

/***************************************************************************************************
*                                                                                                  *
* init() -- initialize Hmm class. Set values of the class parameters.                              *
*                                                                                                  *
* Arguments:                                                                                       *
*     const char* fn -- file name with HMM definition                                              *
*                                                                                                  *
***************************************************************************************************/
void Hmm::init(const string fn)
{
	this->src_fn        = fn;
	this->comment_char  = '#';
	this->spaces        = " \t\n\r";
	this->cur_file_line = 0;
	this->max_seq_offset= 1000000;
}

/***************************************************************************************************
*                                                                                                  *
* finalize() -- finalize Hmm object creation                                                       *
*                                                                                                  *
***************************************************************************************************/
void Hmm::finalize()
{
	int i;
	
	// Calculate conditional probabilities
	for(i = 0; i < (int)this->emission_list.size(); ++i)
		this->calculate_cond_emission_probs(i);
}

/***************************************************************************************************
*                                                                                                  *
* get_emission_prob_log() -- wrapper around get_emission_prob() that returns log of the            *
*                            probability and handles the case when probability is 0.               *
*                                                                                                  *
***************************************************************************************************/
double Hmm::get_emission_prob_log(int emission_id, const string &emission_str)
{
	double prob = this->get_emission_prob(emission_id, emission_str);
	if( prob == 0 )
		return ZERO_LOG_VALUE;
	else
		return log(prob);
}

/***************************************************************************************************
*                                                                                                  *
* get_emission_prob() -- returns conditional emission probability for a given state_id,            *
*                        emission_str and period. This function also handles cases where the       *
*                        shorter emission_str is passed.                                           *
*                                                                                                  *
* Returns:                                                                                         *
*      double prob -- conditional emission probability                                             *
*                                                                                                  *
***************************************************************************************************/
double Hmm::get_emission_prob(int emission_id, const string &emission_str)
{
	int len;
	map<string, double>::const_iterator prob;
	
	len = (int)emission_str.length();
	if( len == 0 || len > this->emission_list[emission_id].order+1 )
	{
		cerr << "Wrong emission string = '" + emission_str + "' (HMM order of emission item "
		     << emission_id << " is " << this->emission_list[emission_id].order << ")";
		exit(1);
	}
	else
	{
		prob = this->emission_list[emission_id].cond_probs.find( emission_str );
		if( prob == this->emission_list[emission_id].cond_probs.end() )
		{
			return 0.0;
		}
		else
		{
			return prob->second;
		}
	}
}

/***************************************************************************************************
*                                                                                                  *
* get_transition_prob_log() -- wrapper around get_transition_prob() that returns log of the        *
*                              probability and handles the case when probability is 0.             *
*                                                                                                  *
***************************************************************************************************/
double Hmm::get_transition_prob_log(int from_state_id, int to_state_id, const string &emission_str)
{
	double prob = this->get_transition_prob(from_state_id, to_state_id, emission_str);
	if( prob == 0 )
		return ZERO_LOG_VALUE;
	else
		return log(prob);
}

/***************************************************************************************************
*                                                                                                  *
* get_transition_prob() -- returns transition probability for a given pair                         *
*                          from_state_id - to_state_id. Returns 0 if such a transition is not      *
*                          defined in the HMM.                                                     *
*                                                                                                  *
* Arguments:                                                                                       *
*     int    from_state_id -- from state                                                           *
*     int    to_state_id   -- to state                                                             *
*     string emission_str  -- current emission string                                              *
*                                                                                                  *
* Returns:                                                                                         *
*      double prob -- conditional emission probability                                             *
*                                                                                                  *
***************************************************************************************************/
double Hmm::get_transition_prob(int from_state_id, int to_state_id, const string &emission_str)
{
	map<int, Transition_item>::const_iterator trans;
	map<string, double>::const_iterator exception;
	
	trans = this->state_list[from_state_id].transitions.find( to_state_id );
	if( trans == this->state_list[from_state_id].transitions.end() )
		return 0.0;
	else
	{
		// Try to find transition probability exception
		exception = trans->second.emiss_str_exceptions.find( emission_str );
		if( exception == trans->second.emiss_str_exceptions.end() )
			return trans->second.prob;   // Exception is not found -- return usual transition probability
		else
			return exception->second;    // Exception found -- return EXCEPTIONAL transition probability
	}
}

/***************************************************************************************************
*                                                                                                  *
* get_num_states() -- returns total number of states in the HMM                                    *
*                                                                                                  *
***************************************************************************************************/
int Hmm::get_num_states()
{
	return (int)this->state_list.size();
}

/***************************************************************************************************
*                                                                                                  *
* get_max_order() -- returns maximum order among all emission items of the HMM                     *
*                                                                                                  *
***************************************************************************************************/
int Hmm::get_max_order()
{
	int i, max_order;
	
	max_order = this->emission_list[0].order;
	for(i = 1; i < (int)this->emission_list.size(); ++i)
	{
		if(this->emission_list[i].order > max_order)
		{
			max_order = this->emission_list[i].order;
		}
	}

	if(max_order > this->max_seq_offset)
		max_order = this->max_seq_offset;   // Start with 3rd letter of the sequence in case first 3 letters is a start codon

	return max_order;
}

/***************************************************************************************************
*                                                                                                  *
* get_state_name() -- returns state name by state index                                            *
*                                                                                                  *
***************************************************************************************************/
string Hmm::get_state_name(int state_i)
{
	if(state_i < 0 || state_i >= (int)this->state_list.size())
	{
		cerr << "Wrong state index = " << state_i << endl;
		exit(1);
	}
	
	return this->state_list[state_i].name;
}

/***************************************************************************************************
*                                                                                                  *
* get_emiss_name() -- returns emission name by emission id                                         *
*                                                                                                  *
***************************************************************************************************/
string Hmm::get_emiss_name(int emiss_i)
{
	if(emiss_i < 0 || emiss_i >= (int)this->emission_list.size())
	{
		cerr << "Wrong emission index = " << emiss_i << endl;
		exit(1);
	}
	
	return this->emission_list[emiss_i].name;
}

/***************************************************************************************************
*                                                                                                  *
* get_emission_str() -- returns emission string from a given sequence and pos_i with respect       *
*                       to HMM order.                                                              *
*                                                                                                  *
* Returns:                                                                                         *
*     emission_str  --  emission string in uppercase for a given sequence position                 *
*                                                                                                  *
***************************************************************************************************/
string Hmm::get_emission_str(const Sequence &seq, int pos_i, int emission_id)
{
	string emission_str;
	int order;
	
	if( pos_i < 0 || pos_i >= (int)seq.seq.length())
	{
		cerr << "Hmm::get_emission_str: wrong pos_i = " << pos_i << endl;
		exit(1);
	}
	
	order = this->emission_list[emission_id].order;
	if( pos_i < order )                                  // For the very beginning of sequence
		emission_str = seq.seq.substr(0, pos_i+1);   
	else
		emission_str = seq.seq.substr(pos_i-order, order+1);

	return uc_string(emission_str);
}

/***************************************************************************************************
*                                                                                                  *
* get_period() -- returns period calculated based on pos_i and state periodicity                   *
*                                                                                                  *
***************************************************************************************************/
int Hmm::get_period(int pos_i, int state_i)
{
	return pos_i%this->state_list[state_i].periodicity;
}

/***************************************************************************************************
*                                                                                                  *
* can_start_with() -- can we start HMM path with the state state_i. Namely the method              *
*                     checks whether the state_i is in dont_start_with array                       *
*                                                                                                  *
***************************************************************************************************/
bool Hmm::can_start_with(int state_i)
{
	int i;
	for(i=0; i < (int)this->opts.dont_start_with.size(); ++i)
	{
		if(this->opts.dont_start_with[i] == state_i)
			return false;
	}
	return true;
}

/***************************************************************************************************
*                                                                                                  *
* can_end_with() -- can we end HMM path with the state state_i. Namely the method                  *
*                   checks whether the state_i is in dont_end_with array                           *
*                                                                                                  *
***************************************************************************************************/
bool Hmm::can_end_with(int state_i)
{
	int i;
	for(i=0; i < (int)this->opts.dont_end_with.size(); ++i)
	{
		if(this->opts.dont_end_with[i] == state_i)
			return false;
	}
	return true;
}

/***************************************************************************************************
*                                                                                                  *
* run_viterbi() -- find the most probable sequence of states by Viterbi algorithm                  *
*                                                                                                  *
* Arguments:                                                                                       *
*     Sequence seq    -- sequence to allply Viterbi algorithm to                                   *
*                                                                                                  *
* Returns:                                                                                         *
*     Path path -- full information about all pathes Viterbi considered                            *
*                                                                                                  *
***************************************************************************************************/
void Hmm::run_viterbi(Sequence &seq, Path& path)
{
	int pos_i; // End position of sliding window in sequence coordinates
	
	pos_i = this->do_viterbi_first_steps(seq, this->get_max_order(), path);
	
	for(; pos_i < seq.length(); ++pos_i)
		this->do_viterbi_step(seq, pos_i, path);
	
	this->traceback(path);
}

/***************************************************************************************************
*                                                                                                  *
* do_viterbi_first_step() -- add values to two path's tables for sequence position = 0.            *
*                            For ptr table the same value will be used - const int BEGIN_STATE_ID  *
*                            that corresponds to BEGIN state.                                      *
*                                                                                                  *
* Returns:                                                                                         *
*     int next_pos  --  which step we should proceed with                                          *
*                                                                                                  *
***************************************************************************************************/
int Hmm::do_viterbi_first_steps(Sequence &seq, const int pos_i, Path &path)
{
	int state_i, period, emission_id;
	double sum;
	string emission_str;
	
	for(state_i = 0; state_i < this->get_num_states(); ++state_i)
	{
		path.all_steps[pos_i][state_i].ptr = BEGIN_STATE_ID;  // Previous state is always BEGIN_STATE
		
		if( this->can_start_with(state_i) )                   // Are we allowed to start HMM path with this state?
		{
			period       = this->get_period(pos_i, state_i);
			emission_id  = this->state_list[state_i].emission_ids[period];
			emission_str = this->get_emission_str(seq, pos_i, emission_id);
			sum          = log(1.0/(double)this->get_num_states()) + this->get_emission_prob_log(emission_id, emission_str);
			path.all_steps[pos_i][state_i].sum = sum;
		}
		else
		{
			path.all_steps[pos_i][state_i].sum = ZERO_LOG_VALUE; // Assign very negative value if we are not allowed
		}
	}
	
	return(pos_i + 1);
}

/***************************************************************************************************
*                                                                                                  *
* do_viterbi_step() -- add values to two path's tables for a sequence positon pos_i                *
*                                                                                                  *
***************************************************************************************************/
void Hmm::do_viterbi_step(Sequence &seq, const int pos_i, Path &path)
{
	int    to_state, from_state, period, emission_id;
	double sum;
	string emission_str;
	
	for(to_state = 0; to_state < this->get_num_states(); ++to_state)
	{
		if(pos_i == seq.length()-1 && !this->can_end_with(to_state))
		{
			from_state = 1;                           // Any fake state -- it doesn't matter
			path.all_steps[pos_i][to_state].ptr = from_state;
			path.all_steps[pos_i][to_state].sum = path.all_steps[pos_i-1][from_state].sum + ZERO_LOG_VALUE;
		}
		else
		{
			period       = this->get_period(pos_i, to_state);
			emission_id  = this->state_list[to_state].emission_ids[period];
			emission_str = this->get_emission_str(seq, pos_i, emission_id);
			from_state   = this->get_best_prev_state(path.all_steps[pos_i-1], to_state, emission_str);   // Find the best state to come from
			sum          = this->get_transition_prob_log(from_state, to_state, emission_str) + this->get_emission_prob_log(emission_id, emission_str);
			
			path.all_steps[pos_i][to_state].ptr = from_state;
			path.all_steps[pos_i][to_state].sum = path.all_steps[pos_i-1][from_state].sum + sum;
		}
	}
}

/***************************************************************************************************
*                                                                                                  *
* get_best_prev_state() -- returns the best previous state for the given to_state                  *
*                                                                                                  *
* Arguments:                                                                                       *
*     Path_step* prev_step    -- data about previous step                                          *
*     int        to_state     -- state id where we are going to come                               *
*     string     emission_str -- current emission string                                           *
*                                                                                                  *
* Returns:                                                                                         *
*     int best_prev_state -- id of the previous state to go from which is the best option          *
*                                                                                                  *
***************************************************************************************************/
int Hmm::get_best_prev_state(const Path_step* prev_step, const int to_state, const string &emission_str)
{
	int    prev_state, best_prev_state;
	double sum, best_sum;
	
	best_sum        = prev_step[0].sum + this->get_transition_prob_log(0, to_state, emission_str);
	best_prev_state = 0;
	for(prev_state = 1; prev_state < this->get_num_states(); ++prev_state)
	{
		sum = prev_step[prev_state].sum + this->get_transition_prob_log(prev_state, to_state, emission_str);
		if(sum > best_sum)
		{
			best_sum        = sum;
			best_prev_state = prev_state;
		}
	}
	
	return best_prev_state;
}

/***************************************************************************************************
*                                                                                                  *
* traceback() -- does traceback for a given path and puts the traced sequence of states in         *
*                path.states_seq and path.emiss_seq                                                *
* Arguments:                                                                                       *
*     Path path -- Path object with all_steps array filled completely                              *
*                                                                                                  *
***************************************************************************************************/
void Hmm::traceback(Path &path)
{
	int    cur_state, state_i, pos_i, period;
	double best_sum;
	
	// Find the best state at the end
	cur_state = 0;
	best_sum  = path.all_steps[path.num_steps-1][0].sum;
	for(state_i = 1; state_i < this->get_num_states(); ++state_i)
	{
		if(path.all_steps[path.num_steps-1][state_i].sum > best_sum)
		{
			best_sum   = path.all_steps[path.num_steps-1][state_i].sum;
			cur_state  = state_i;
		}
	}
	
	for(pos_i = path.num_steps-1; pos_i >= 0; --pos_i)
	{
		period = this->get_period(pos_i, cur_state);
		path.states_seq[pos_i] = cur_state;
		path.emiss_seq [pos_i] = this->state_list[cur_state].emission_ids[period];
		cur_state = path.all_steps[pos_i][cur_state].ptr;
		
		if(cur_state == BEGIN_STATE_ID)
		{
			path.first_step_i = pos_i;
			break;
		}
	}
}

/***************************************************************************************************
*                                                                                                  *
* load_header_info() -- loads information from FILE_INFO and HMM_MODEL sections                    *
*                                                                                                  *
* Arguments:                                                                                       *
*     istream in -- istream object for the opened hmm_def file                                     *
*                                                                                                  *
***************************************************************************************************/
void Hmm::load_header_info(ifstream& in)
{
	string line, tmp;
	
	//>>>>>>>>>> FILE_INFO section <<<<<<<<<<
	this->next_line(in, line);
	this->check_line(line, "<FILE_INFO>");

	this->next_line(in, line);
	this->parse_name_value(line, tmp, "version");
	this->file_version = strtod(tmp.c_str(), NULL);
	
	this->next_line(in, line);
	this->parse_name_value(line, this->file_svn_str, "svn_str");
	
	this->next_line(in, line);
	this->check_line(line, "</FILE_INFO>");

	//>>>>>>>>>> HMM_MODEL section <<<<<<<<<<
	this->next_line(in, line);
	this->check_line(line, "<HMM_MODEL>");
	
	this->next_line(in, line);
	this->parse_name_value(line, this->model_name, "name");

	this->next_line(in, line);
	this->parse_name_value(line, this->model_name, "descr");
	
	this->next_line(in, line);
	this->check_line(line, "</HMM_MODEL>");
}

/***************************************************************************************************
*                                                                                                  *
* load_options() -- loads information from OPTIONS section                                         *
*                                                                                                  *
* Arguments:                                                                                       *
*     istream in -- istream object for the opened hmm_def file                                     *
*                                                                                                  *
***************************************************************************************************/
void Hmm::load_options(ifstream& in)
{
	string line, name, tmp;
	
	this->next_line(in, line);
	this->check_line(line, "<OPTIONS>");

	while( this->next_line(in,line) )
	{
		if(line.compare("</OPTIONS>") == 0)
			break;
		else                                // We still inside options section and we have name value pair now
		{
			name = this->parse_name_value(line, tmp);
			if( name.compare("dont_start_with") == 0 )
				this->parse_i_string(tmp, this->opts.dont_start_with);
			else if( name.compare("dont_end_with") == 0 )
				this->parse_i_string(tmp, this->opts.dont_end_with);
			else
				this->wrong_fmt_die("Unknown option '"+name+"'");
		}
	}
}

/***************************************************************************************************
*                                                                                                  *
* load_emission_list() -- loads information from EMISSION_LIST section                             *
*                                                                                                  *
* Arguments:                                                                                       *
*     istream in -- istream object for the opened hmm_def file                                     *
*                                                                                                  *
***************************************************************************************************/
void Hmm::load_emission_list(ifstream& in)
{
	int id, expected_id;
	string line;
	stringstream ss_id, ss_expected_id;
	
	this->next_line(in, line);
	this->check_line(line, "<EMISSION_LIST>");
	
	expected_id = 0;
	while( this->next_line(in,line) )
	{
		if(line.compare("</EMISSION_LIST>") == 0)
			break;
		else if( line.compare("<ITEM>") == 0 )
		{
			id = this->load_emission_item( in );
			if(id != expected_id)
			{
				ss_id << id;
				ss_expected_id << expected_id;
				this->wrong_fmt_die("unexpected emission item id = " + ss_id.str() + ". I expect " + ss_expected_id.str() + " here");
			}
			expected_id++;
		}
		else
			this->wrong_fmt_die("unexpected line '"+line+"', I expect '</EMISSION_LIST>' or '<ITEM>' here");
	}
}

/***************************************************************************************************
*                                                                                                  *
* load_emission_item() -- loads information from a single ITEM section of EMISSION_LIST            *
*                                                                                                  *
* Arguments:                                                                                       *
*     istream in -- istream object for the opened hmm_def file with just read <ITEM> line          *
*                                                                                                  *
* Returns:                                                                                         *
*     int id -- id of the just loaded item                                                         *
*                                                                                                  *
***************************************************************************************************/
int Hmm::load_emission_item(ifstream& in)
{
	double        prob;
	string        line, key, tmp;
	Emission_item item;
	
	this->next_line(in, line);
	this->parse_name_value(line, tmp, "id");
	item.id = atoi( tmp.c_str() );
	
	this->next_line(in, line);
	this->parse_name_value(line, item.name, "name");
	
	this->next_line(in, line);
	this->check_line(line, "<PROBABILITIES>");
	
	while( this->next_line(in, line) )
	{
		if( line.compare("</PROBABILITIES>") == 0 )
			break;
		key  = this->parse_name_value(line, tmp);
		prob = strtod(tmp.c_str(), NULL);           // change type to double
		uc_string(key);                             // translate emission string in uppercase
		
		// Check if we have already saw such key
		if( item.init_freqs.find(key) == item.init_freqs.end() )
			item.init_freqs.insert( pair<string, double>(key, prob) );   // Everything is OK: no such key
		else
			this->wrong_fmt_die("key '" + key + "' is duplicated in <PROBABILITIES> section");
	}
	
	this->emission_list.push_back(item);
	
	this->next_line(in, line);
	this->check_line(line, "</ITEM>");
	
	return item.id;
}

/***************************************************************************************************
*                                                                                                  *
* load_emission_list() -- loads information from STATE_LIST section                                *
*                                                                                                  *
* Arguments:                                                                                       *
*     istream in -- istream object for the opened hmm_def file                                     *
*                                                                                                  *
***************************************************************************************************/
void Hmm::load_state_list(ifstream& in)
{
	int id, expected_id;
	string line;
	stringstream ss_id, ss_expected_id;
	
	this->next_line(in, line);
	this->check_line(line, "<STATE_LIST>");
	
	expected_id = 0;
	while( this->next_line(in,line) )
	{
		if(line.compare("</STATE_LIST>") == 0)
			break;
		else if( line.compare("<ITEM>") == 0 )
		{
			id = this->load_state_item( in );
			if(id != expected_id)
			{
				ss_id << id;
				ss_expected_id << expected_id;
				this->wrong_fmt_die("unexpected state item id = " + ss_id.str() + ". I expect " + ss_expected_id.str() + " here");
			}
			expected_id++;
		}
		else
			this->wrong_fmt_die("unexpected line '"+line+"', I expect '</STATE_LIST>' or '<ITEM>' here");
	}
}

/***************************************************************************************************
*                                                                                                  *
* load_state_item() -- loads information from a single ITEM section of STATE_LIST                  *
*                                                                                                  *
* Arguments:                                                                                       *
*     istream in -- istream object for the opened hmm_def file with just read <ITEM> line          *
*                                                                                                  *
* Returns:                                                                                         *
*     int id -- id of the just loaded item                                                         *
*                                                                                                  *
***************************************************************************************************/
int Hmm::load_state_item(ifstream& in)
{
	string     line, tmp, id;
	State_item item;
	
	this->next_line(in, line);
	this->parse_name_value(line, tmp, "id");
	item.id = atoi( tmp.c_str() );
	
	this->next_line(in, line);
	this->parse_name_value(line, item.name, "name");
	
	this->next_line(in, line);
	this->parse_name_value(line, tmp, "periodicity");
	item.periodicity = atoi( tmp.c_str() );
	
	this->next_line(in, line);
	this->parse_name_value(line, tmp, "emission_set");
	this->parse_i_string(tmp, item.emission_ids);
	
	this->state_list.push_back( item );
	
	this->next_line(in, line);
	this->check_line(line, "</ITEM>");
	
	return item.id;
}


/***************************************************************************************************
*                                                                                                  *
* load_transition_list() -- loads information from TRANSITION_LIST section                         *
*                                                                                                  *
* Arguments:                                                                                       *
*     istream in -- istream object for the opened hmm_def file                                     *
*                                                                                                  *
***************************************************************************************************/
void Hmm::load_transition_list(ifstream& in)
{
	string line;
	
	this->next_line(in, line);
	this->check_line(line, "<TRANSITION_LIST>");
	
	while( this->next_line(in,line) )
	{
		if(line.compare("</TRANSITION_LIST>") == 0)
			break;
		else if( line.compare("<ITEM>") == 0 )
			this->load_transition_item( in );
		else
			this->wrong_fmt_die("unexpected line '"+line+"', I expect '</TRANSITION_LIST>' or '<ITEM>' here");
	}
}

/***************************************************************************************************
*                                                                                                  *
* load_transition_item() -- loads information from a single ITEM section of TRANSITION_LIST        *
*                                                                                                  *
* Arguments:                                                                                       *
*     istream in -- istream object for the opened hmm_def file with just read <ITEM> line          *
*                                                                                                  *
***************************************************************************************************/
void Hmm::load_transition_item(ifstream& in)
{
	string          line, tmp;
	Transition_item item;

	this->next_line(in, line);
	this->parse_name_value(line, tmp, "from_state");
	item.from_state = atoi( tmp.c_str() );
	if(item.from_state < 0 || item.from_state >= this->get_num_states())
		this->wrong_fmt_die("wrong from_state id = '" + tmp + "'");
	
	this->next_line(in, line);
	this->parse_name_value(line, tmp, "to_state");
	item.to_state = atoi( tmp.c_str() );
	if(item.to_state < 0 || item.to_state >= this->get_num_states())
		this->wrong_fmt_die("wrong to_state id = '" + tmp + "'");
	
	this->next_line(in, line);
	this->parse_name_value(line, tmp, "probability");
	item.prob = strtod(tmp.c_str(), NULL);
	
	this->next_line(in, line);
	if( line.compare("<EMISSION_STR_EXCEPTIONS>") == 0 )
	{
		this->load_transition_item_emiss_str_exceptions(in, item);
		this->next_line(in, line);
	}
	
	this->check_line(line, "</ITEM>");
	this->state_list[ item.from_state ].transitions.insert(pair <int, Transition_item> (item.to_state, item));
}

/***************************************************************************************************
*                                                                                                  *
* load_transition_item_emiss_str_exceptions() -- loads EMISSION_STR_EXCEPTIONS for a particular    *
*                                                ITEM from TRANSITION_LIST                         *
*                                                                                                  *
* Arguments:                                                                                       *
*     istream         in   -- istream object for the opened hmm_def file with just read <ITEM> line*
*     Transition_item item -- Transition_item objet to add exceptions to                           *
*                                                                                                  *
***************************************************************************************************/
void Hmm::load_transition_item_emiss_str_exceptions(std::ifstream &in, Transition_item &item)
{
	string key, tmp, line;
	double prob;

	while( this->next_line(in, line) )
	{
		if( line.compare("</EMISSION_STR_EXCEPTIONS>") == 0 )
			break;
		key  = this->parse_name_value(line, tmp);
		prob = strtod(tmp.c_str(), NULL);           // change type to double
		uc_string(key);                             // translate emission string in uppercase
		
		// Check if we have already saw such key
		if( item.emiss_str_exceptions.find(key) == item.emiss_str_exceptions.end() )
			item.emiss_str_exceptions.insert( pair<string, double>(key, prob) );   // Everything is OK: no such key
		else
			this->wrong_fmt_die("key '" + key + "' is duplicated in <PROBABILITIES> section");
	}
}

/***************************************************************************************************
*                                                                                                  *
* next_line() -- reads next line from the hmm_def file. The line is trimmed on both ends.          *
*                Comment and empty lines are ignored.                                              *
*                                                                                                  *
* Arguments:                                                                                       *
*     istream in -- istream object for the opened hmm_def file                                     *
*     string str -- string to store the line                                                       *
*                                                                                                  *
***************************************************************************************************/
int Hmm::next_line(istream& in, string& str)
{
	while(1)
	{
		getline(in, str);
		this->cur_file_line++;
		
		if( in.eof() )
			this->wrong_fmt_die( "unexpected end of file" );
		
		trim(str, this->spaces);
		
		if(str.length() > 0 && str[0] != this->comment_char)
			break;
	}
	return 1;
}

/***************************************************************************************************
*                                                                                                  *
* parse_name_value() -- parses simple name-value string and stores value in the variable value     *
*                                                                                                  *
* Arguments:                                                                                       *
*     string line       -- line to parse. Format: '<name> <spaces> = <spaces> <value>'.            *
*                          <spaces> is 0 or more spaces that are specified in this->spaces         *
*                          variable. <name> must NOT contain '=' character.                        *
*     string value      -- string to store <value>                                                 *
*     string valid_name -- (optional) string to validate <name>. If valid_name is "-" validation   *
*                          is not performed. Default value is "-".                                 *
*                                                                                                  *
* Returns:                                                                                         *
*     string name       -- the <name>                                                              *
*                                                                                                  *
***************************************************************************************************/
string Hmm::parse_name_value(string& line, string& value, const string valid_name)
{
	int pos;
	string name;
	
	trim(line, this->spaces); // Remove spaces from the beginnign and end - just in case...
	
	pos = (int)line.find_first_of("=");
	if(pos == -1)
		this->wrong_fmt_die("line '"+line+"' doesn't contain '=' character");
	
	name = line.substr(0, pos);
	trim(name, this->spaces);
	
	value = line.substr(pos+1);
	trim(value, this->spaces);
	
	if( valid_name.compare("-") != 0 && valid_name.compare( name ) != 0 )   // Validate name
		this->wrong_fmt_die("unexpected name in the line '"+line+"', I expect name '"+valid_name+"' here");
	
	return name;
}

/***************************************************************************************************
*                                                                                                  *
* parse_i_string() -- parses string of integers and puts them into the given vector                *
*                                                                                                  *
* Arguments:                                                                                       *
*     string str -- string of integers separated by commas and (optional) spaces, ex: "0, 1,2"     *
*                                                                                                  *
* Returns:                                                                                         *
*     int t_num  --  total num integers derived from the string                                      *
*                                                                                                  *
***************************************************************************************************/
int Hmm::parse_i_string(string str, vector <int> &v)
{
	int pos, t_num;
	string str_i;
	
	t_num = 0;
	while(1)
	{
		pos = (int)str.find_first_of(",");
		if(pos == -1)
			break;
		
		str_i = str.substr(0, pos);
		trim(str_i, this->spaces);
		v.push_back( atoi(str_i.c_str()) );
		
		str = str.substr(pos+1);
		trim(str, this->spaces);
		
		t_num++;
	}
	v.push_back( atoi(str.c_str()) );
	
	return(t_num+1);
}

/***************************************************************************************************
*                                                                                                  *
* check_line() -- compares observed_str with expected_str. If they are not equal terminates        *
*                 program by calling wrong_fmt_die() function.                                     *
*                                                                                                  *
* Arguments:                                                                                       *
*     string observed_str -- string to check                                                       *
*     string expected_str -- expected right string                                                 *
*                                                                                                  *
***************************************************************************************************/
void Hmm::check_line(const string observed_str, const string expected_str)
{
	if(observed_str.compare( expected_str ) != 0)
		this->wrong_fmt_die("unexpected line '"+observed_str+"', I expect '"+expected_str+"' here");
}

/***************************************************************************************************
*                                                                                                  *
* wrong_fmt_die() -- prints given message to STDERR and  terminates the program                    *
*                                                                                                  *
* Arguments:                                                                                       *
*     string message -- what's wrong with file format                                              *
*                                                                                                  *
***************************************************************************************************/
void Hmm::wrong_fmt_die(const string message)
{
	cerr << "[" + this->src_fn + " : line " << this->cur_file_line << "] Wrong format: " + message;
	exit(1);
}

/***************************************************************************************************
*                                                                                                  *
* check_hmm() -- validates just loaded HMM and terminates the program if HMM is wrong.             *
*                This function also set some values to this.                                       *
*                                                                                                  *
***************************************************************************************************/
void Hmm::check_hmm()
{
	int i;
	
	
	
	for(i = 0; i < (int)this->emission_list.size(); ++i)
		this->check_emission_item(i);
	
	for(i = 0; i < (int)this->state_list.size(); ++i)
		this->check_state_item(i);
}

/***************************************************************************************************
*                                                                                                  *
* check_emission_item() -- validates this->options. Namely it checks that:                         *
*                             * all dont_start_with states must exist                              *
*                             * all dont_end_with states must exist                                *
*                                                                                                  *
***************************************************************************************************/
void Hmm::check_opts()
{
	int i;
	
	for(i=0; i < (int)this->opts.dont_start_with.size(); ++i)
	{
		if( this->opts.dont_start_with[i] >= (int)this->state_list.size() )
		{
			cerr << "State with id "<< this->opts.dont_start_with[i] <<" doesn't exist (see dont_start_with option)";
			exit(1);
		}
	}

	for(i=0; i < (int)this->opts.dont_end_with.size(); ++i)
	{
		if( this->opts.dont_end_with[i] >= (int)this->state_list.size() )
		{
			cerr << "State with id "<< this->opts.dont_end_with[i] <<" doesn't exist (see dont_end_with option)";
			exit(1);
		}
	}
}

/***************************************************************************************************
*                                                                                                  *
* check_emission_item() -- validates given emission item and terminates the program if the item    *
*                          is not valid. Also sets value to this->emission_list[id].order.         *
*                          Properties of a valid item:                                             *
*                             * all emission strings have the same order                           *
*                             * item id is a non negative integer                                  *
*                             * sum of all probabilities of the item is equal to 1                 *
*                                                                                                  *
* Arguments:                                                                                       *
*     int item -- emission item id to check                                                        *
*                                                                                                  *
***************************************************************************************************/
void Hmm::check_emission_item(int id)
{
	double sum;
	bool first_loop;
	map <string, double>::const_iterator freq_i;
	
	sum = 0;
	first_loop = true;
	for(freq_i = this->emission_list[id].init_freqs.begin(); freq_i != this->emission_list[id].init_freqs.end(); ++freq_i)
	{
		// Check emission order
		if( first_loop )
		{
			this->emission_list[id].order = (int)freq_i->first.length()-1;
			first_loop = false;
		}
		else if( (int)freq_i->first.length() != this->emission_list[id].order+1 )
		{
			cerr << "Emission item with id=" << id << " has wrong order of emission string: '"
			     << freq_i->first << "' (right order is " << this->emission_list[id].order << ")" << endl;
			exit(1);
		}
		
		// Check frequency value
		if( freq_i->second < 0 )
		{
			cerr << "In emission item " << id << " emission probability for '" << freq_i->first << "' = "
			     << freq_i->second << " (negative!)" << endl;
			exit(1);
		}
		sum += freq_i->second;
	}
	
	if( fabs(sum-1.0) > 0.05 )
	{
		cerr << "Probability sum for the item with id=" << id << " = " << sum << " (!= 1)" << endl;
		exit(1);
	}
}

/***************************************************************************************************
*                                                                                                  *
* check_state_item() -- validates given state item and terminates the program if the item is not   *
*                       valid. Properties of the valid item:                                       *
*                           * Sum of all outcoming TRANSITION probabilities = 1                    *
*                           * Number of emission_ids must be equal to state periodicity            *
*                                                                                                  *
* Arguments:                                                                                       *
*     int item -- emission item id to check                                                        *
*                                                                                                  *
***************************************************************************************************/
void Hmm::check_state_item(int id)
{
	double sum;
	map <int, Transition_item>::const_iterator trans;
	
	// Sum of all outcoming TRANSITION probabilities for each state = 1
	sum = 0;
	for(trans = this->state_list[id].transitions.begin(); trans != this->state_list[id].transitions.end(); ++trans)
	{
		if( trans->second.prob < 0 )
		{
			cerr << "In state item " << id << " transition probability " << trans->second.from_state << " -> "
			     << trans->second.to_state << " is " << trans->second.prob << " (negative!)" << endl;
			exit(1);
		}
		sum += trans->second.prob;
	}
	
	if( fabs(sum-1.0) > 0.0001 )
	{
		cerr << "Sum of all transition probabilities for state = " << id << " = " << sum << " (!= 1)" << endl;
		exit(1);
	}
	
	// Number of emission_ids must be equal to state periodicity
	if( this->state_list[id].emission_ids.size() != this->state_list[id].periodicity )
	{
		cerr << "emission_set size != periodicity for state = " << id << endl;
		exit(1);
	}
}

/***************************************************************************************************
*                                                                                                  *
* calculate_cond_emission_probs() -- calculates conditional emission probabilitis from the initial *
*                                    frequencies read from the file                                *
*                                                                                                  *
* Description:                                                                                     *
*     In the hmm_def file user specified just frequencies for emission string. Sum of all          *
*     all frequencies is equal to 1. There frequencies are stored in the init_freqs map. Here we   *
*     we calculate CONDITONAL PROBABILITIES from the frequencies and store probabilites obtained   *
*     cond_probs map.                                                                              *
*                                                                                                  *
* Arguments:                                                                                       *
*     int emiss_id -- emission item id to calculate conditional probabilities for                  *
*                                                                                                  *
***************************************************************************************************/
void Hmm::calculate_cond_emission_probs(int emiss_id)
{
	int prefix_len;
	double prob, prefix_freq;
	string prefix;
	map<string, double> prefix_freqs;
	map<string, double>::iterator prefix_item;
	map<string, double>::const_iterator freq;
	
	if( this->emission_list[emiss_id].order == 0 )   // For zero order HMM conditional probabilities equal to frequencies
	{
		this->emission_list[emiss_id].cond_probs = this->emission_list[emiss_id].init_freqs;
		return;
	}
	
	// Calculate all frequencies for the prefixes of all length
	prefix_freqs = this->emission_list[emiss_id].init_freqs;       // Prefix Length = this->order + 1
	for(prefix_len = this->emission_list[emiss_id].order; prefix_len > 0; --prefix_len)
	{
		for(freq = prefix_freqs.begin(); freq != prefix_freqs.end(); ++freq)
		{
			if(freq->first.length() != prefix_len+1)   // Here we need other prefixes that longer on one letter only
				continue;
			
			prefix      = freq->first.substr(0, prefix_len);
			prefix_item = prefix_freqs.find(prefix);
			if( prefix_item == prefix_freqs.end() )            // We didn't see such prefix yet
				prefix_freqs.insert( pair <string, double> (prefix, freq->second) ); 
			else
				prefix_item->second += freq->second;
		}
	}
	prefix_freqs.insert( pair <string, double> ("", 1) ); 
	
	// Calculate conditional probabilities itself: cond_prob = FREQ(prefix+letter) / FREQ(prefix)
	for(prefix_len = this->emission_list[emiss_id].order; prefix_len >= 0; --prefix_len)
	{
		for(freq = prefix_freqs.begin(); freq != prefix_freqs.end(); ++freq)
		{
			if(freq->first.length() != prefix_len+1)   // prefix + letter => prefix_len+1
				continue;
			
			prefix      = freq->first.substr(0, prefix_len);
			prefix_freq = prefix_freqs.find(prefix)->second;
			prob        = prefix_freq == 0 ? 0 : freq->second / prefix_freq;
			this->emission_list[emiss_id].cond_probs.insert( pair <string, double> (freq->first, prob) );
		}
	}
}

bool Hmm::is_coding_state(int state_i)
{
	if( this->get_state_name(state_i).compare("0")   == 0 ||
	    this->get_state_name(state_i).compare("1")   == 0 ||
	    this->get_state_name(state_i).compare("2")   == 0 ||
	    this->get_state_name(state_i).compare("0+1") == 0 ||
	    this->get_state_name(state_i).compare("0+2") == 0 ||
	    this->get_state_name(state_i).compare("1+1") == 0 ||
	    this->get_state_name(state_i).compare("1+2") == 0 ||
	    this->get_state_name(state_i).compare("2+1") == 0 ||
	    this->get_state_name(state_i).compare("2+2") == 0
	)
		return true;
	else
		return false;
}

bool Hmm::is_end_state(int state_i)
{
	if( this->get_state_name(state_i).find("stop") == -1 )
		return false;
	else
		return true;
}

bool Hmm::is_start_state(int state_i)
{
	if( this->get_state_name(state_i).find("start") == -1 )
		return false;
	else
		return true;
}

