/***************************************************************************************************
*                                                                                                  *
* Definitions for the Hmm class and other classes it uses.                                         *
* Hmm class is designed to store HMM and to provide functions connected with it.                   *
*                                                                                                  *
* Author: Ivan Antonov (ivan.antonov@gatech.edu)                                                   *
*                                                                                                  *
***************************************************************************************************/

// $Id$

#ifndef _HMM_H_INCLUDED
#define _HMM_H_INCLUDED

#include <math.h>
#include <string>
#include <sstream>
#include <fstream>
#include <vector>
#include <map>
#include <stdlib.h>

#include "base_util.h"
#include "path.h"

using namespace std;

const int BEGIN_STATE_ID = -1;         // Id of the begin state in the HMM. Will be used in path for the first sequence position.
const int ZERO_LOG_VALUE = -1000000;   // Is used instead of log(0)

struct Transition_item
{
	int    from_state;
	int    to_state;
	double prob;
	map    <string, double> emiss_str_exceptions;
};

struct State_item
{
	int    id;
	string name;
	int    periodicity;
	map    <int, Transition_item> transitions;    // List of states TO which current state can go: to_state_id => Transition_item
	vector <int>                  emission_ids;   // List of emission probabilities with respect to periodicity, i.e. 
	                                              // index in emission_ids vector 
};

struct Emission_item
{
	int id;
	int order;                       // Order of the emission item: order = length(emission_string) - 1
	string name;
	map <string, double> init_freqs; // Initial emission frequencies table read from file: emission_string => frequency
	                                 // Sum of all frequencies is equal to 1
	map <string, double> cond_probs; // Conditional probabilities calculated from init_freqs
	                                 // Sum of all probabilities
	                                 //     * for nucleotide 0-order HMM is equal to 1
	                                 //     * for nucleotide 1-order HMM is equal to 4
	                                 //     * for nucleotide 2-order HMM is equal to 16
	                                 //     * for nucleotide n-order HMM is equal to 4^n
};

/* Corresponds to <OPTIONS> section in hmm_def file */
struct Hmm_opts
{
	vector <int> dont_start_with;     // A list of state IDs that can NOT be at the beginning of the predicted FSMark path
	vector <int> dont_end_with;       // A list of state IDs that can NOT be at the end of the predicted FSMark path
};

class Hmm
{
	private:
		//>>>>>>>>>> Technical variables needed while reading HMM from file  <<<<<<<<<<
		int    cur_file_line;  // Current line of the file we read HMM from - is used in error messages
		char   comment_char;   // Character used to start a comment line
		string spaces;         // All spaces that should be trimmed while reading lines from file
		string src_fn;         // Source file name from which HMM model has been read
		
		//>>>>>>>>>> Variables to store HMM <<<<<<<<<<
		Hmm_opts               opts;            // HMM options
		vector <Emission_item> emission_list;   // List of all emission probabilities: emission_item_id => Emission_item
		vector <State_item>    state_list;      // List of all states: state_list[id] = State_item
		
		double file_version;
		string file_svn_str;
		string model_name;
		string model_descr;
		
		//>>>>>>>>>> Methods to read HMM from file <<<<<<<<<<
		void   load_header_info(ifstream& in);
		void   load_options(ifstream& in);
		void   load_emission_list(ifstream& in);
		int    load_emission_item(ifstream& in);
		void   load_state_list(ifstream& in);
		int    load_state_item(ifstream& in);
		void   load_transition_list(ifstream& in);
		void   load_transition_item(ifstream& in);
		void   load_transition_item_emiss_str_exceptions(ifstream& in, Transition_item& item);
		void   check_line(const string observed_str, const string expected_str);
		void   check_hmm();
		void   check_opts();
		void   check_emission_item(int id);
		void   check_state_item(int id);
		int    next_line(istream& in, string& str);
		void   wrong_fmt_die(const string message);
		string parse_name_value(string& line, string& name, const string value = "-");
		int    parse_i_string(string str, vector <int> &v);
		
		//>>>>>>>>>> Viterbi methods <<<<<<<<<<
		bool   can_start_with(int state_i);
		bool   can_end_with(int state_i);
		int    get_period(int pos_i, int state_i);
		int    do_viterbi_first_steps(Sequence &seq, const int pos_i, Path &path);
		void   do_viterbi_step(Sequence &seq, const int pos_i, Path &path);
		int    get_best_prev_state(const Path_step* prev, const int from_state, const string &emission_str);
		void   traceback(Path &path);
		string get_emission_str(const Sequence &seq, int pos_i, int emission_id);
		
		//>>>>>>>>>> Initialization methods <<<<<<<<<<
		void init(const string fn);                       // This member function is called at the very beginning of constructor
		void finalize();                                  // This member function is called at the very end of constructor
		void calculate_cond_emission_probs(int emiss_id);
		
		//>>>>>>>>>> Other methods <<<<<<<<<<
		
	public:
		int max_seq_offset;    // Set this to 2 if input sequence starts with a start codon (ATG)
		                       // and the HMM order is higher than 2nd order
		
		Hmm(const string fn);
		int    get_num_states();
		int    get_max_order();
		string get_state_name(int state_i);
		string get_emiss_name(int emiss_i);
		double get_emission_prob(int emission_id, const string &emission_str);
		double get_emission_prob_log(int emission_id, const string &emission_str);
		double get_transition_prob(int from_state_id, int to_state_id, const string &emission_str);
		double get_transition_prob_log(int from_state_id, int to_state_id, const string &emission_str);
		void   run_viterbi(Sequence &seq, Path &path);
		bool   is_start_state(int state_i);
		bool   is_end_state(int state_i);
		bool   is_coding_state(int state_i);
};

#endif /* _HMM_H_INCLUDED */

