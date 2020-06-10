package GeneTack::FSMark;

use strict;
use warnings;

# $Id$

###
# Antonov Ivan
#
###
#
# $self
#   |
#   |--->{tbl}    --  reference to array of hashes
#   |     tbl[i]
#   |      |
#   |      |--->{pos}
#   |      |--->{let}
#   |      |--->{state}
#   |      |--->{emiss_name}
#   |      |--->{best_score}
#   |
#   |--->{fs}     --  predicted frame shifts
#   |     fs[i]
#   |      |
#   |      |--->{from_state}                  --  
#   |      |--->{to_state}                    --  
#   |      |--->{pos}                         --  last position of sequence in state 'from_state'
#   |      |--->{left_cod_shoulder_len}       --  length of coding region in the left direction from FS
#   |      |--->{right_cod_shoulder_len}      --  length of coding region in the right direction from FS
#   |      |--->{left_gene_nc_shoulder_len}   --  length from FS to the gene left border + length of n/c region to the next gene to the left
#   |      |--->{right_gene_nc_shoulder_len}  --  length from FS to the gene right border + length of n/c region to the next gene to the right
#   |      |--->{type}                        --  '+1' or '-1'
#   |      |
#   |  =============================
#   |  calculate_fs_scores function
#   |  =============================
#   |      |
#   |      |--->{score}
#   |      |--->{lg_fs_score}
#   |      |--->{fs_path_score}
#   |      |--->{wo_fs_path_score}
#   |      |--->{gene_nc_len}              --  length of the frameshifted gene + full length flanking nc shoulders on both sides
#   |
#   |--->{CACHE}
#   |       |
#   |       |
#

use Data::Dumper;
use Carp qw(confess); 
use File::Temp qw(tempfile);
use HTML::Template;

use GeneTack::Lib qw(write_fasta translate_dna_to_prot);

###
# CONSTANTS
our $PRG_PATH = $ENV{HOME}.'/bin/FSMark.exe';
my %FS_TYPES = (
	'01' => '+1',
	'12' => '+1',
	'20' => '+1',
	'10' => '-1',
	'21' => '-1',
	'02' => '-1',
);
our $ORDER_COD_FREQS = {           # Different order HMMs have different order of coding freqs
#	0 => [1,2,3],
#	1 => [3,1,2],
	2 => [2,3,1],
	3 => [1,2,3],
	4 => [3,1,2],
	5 => [2,3,1],
};
my %CODING_STATES = map {$_ => 1} (0, 1, 2);
my $STOP_CODONS   = [qw(TAG TAA TGA)];

###
# CONSTRUCTOR
sub new
{
	my $class = shift;
	my($fsmark_fn, %opts)   = @_;
	$opts{prg_path} ||= $PRG_PATH;
	
	my $self = bless {
		tbl => [],
		fs  => [],
	}, $class;
	
	open(my $fh, '<', $fsmark_fn) or die "Can't open file '$fsmark_fn': $!";
	$self->_load($fh);
	close $fh;
	
	$self->_find_fs();
	
	return $self;
}

sub new_from_seq
{
	my $class = shift;
	my($hmm_def_fn, $seq_fn, %opts) = @_;
	my $fsmark_fn = run_fsmark($hmm_def_fn, $seq_fn, %opts);
	
	my $obj = $class->new($fsmark_fn, %opts);
	
	unlink($fsmark_fn);
	
	return $obj;
}

###
# PRIVATE METHODS
sub _load
{
	my $self = shift;
	my($fh)  = @_;
	
	# Skip header
	$_ = <$fh>;
	
	while( <$fh> )
	{
		s/[\n\r]$//g;
		next if /^\s*$/;
		my($pos, $let, $state, $emiss, $score) = split /\t/;
		push @{$self->{tbl}}, {
			'pos'        => $pos,
			'let'        => $let,
			'state'      => $state,
			'emiss_name' => $emiss,
			'best_score' => $score,
		};
	}
}

sub _find_fs
{
	my $self = shift;
	
	for(my $i = 0; $i < @{$self->{tbl}}-1; $i++)
	{
		my($this_state, $next_state) = ($self->{tbl}[$i]{state}, $self->{tbl}[$i+1]{state});
		next if !defined $this_state || !defined $next_state;
		if($this_state ne $next_state && $CODING_STATES{$this_state} && $CODING_STATES{$next_state})
		{
			push @{$self->{fs}}, {
				'from_state'           => $self->{tbl}[$i]{state},
				'to_state'             => $self->{tbl}[$i+1]{state},
				'pos'                  => $i,
				type                   => $FS_TYPES{$self->{tbl}[$i]{state}.$self->{tbl}[$i+1]{state}},
				left_cod_shoulder_len  => $self->_get_left_shoulder_len( $i ),
				right_cod_shoulder_len => $self->_get_right_shoulder_len( $i ),
				left_gene_nc_shoulder_len  => $self->_get_left_gene_nc_shoulder_len( $i ),
				right_gene_nc_shoulder_len => $self->_get_right_gene_nc_shoulder_len( $i ),
			};
		}
	}
}

sub _get_left_gene_nc_shoulder_len
{
	my $self = shift;
	my($fs_coord) = @_;
	
	my $passed_gene_boundary = 0;
	my $shoulder_len   = 0;
	foreach(my $i = $fs_coord; $i >= 0; $i--)
	{
		my $state = $self->{tbl}[$i]{state};
		last if !defined $state;   # end of sequence
		last if $passed_gene_boundary && $state ne 'nc';
		$passed_gene_boundary = 1 if $state =~ /^start.+nc$/ || $state =~ /^stop.+nc$/;
		$shoulder_len++;
	}
	
	return $shoulder_len;
}

sub _get_right_gene_nc_shoulder_len
{
	my $self = shift;
	my($fs_coord) = @_;

	my $passed_gene_boundary = 0;
	my $shoulder_len   = 0;
	foreach(my $i = $fs_coord+1; $i < @{$self->{tbl}}; $i++)
	{
		my $state = $self->{tbl}[$i]{state};
		last if !defined $state;   # end of sequence
		last if $passed_gene_boundary && $state ne 'nc';
		$passed_gene_boundary = 1 if $state =~ /^start.+nc$/ || $state =~ /^stop.+nc$/;
		$shoulder_len++;
	}
	
	return $shoulder_len;
}

sub _get_left_shoulder_len
{
	my $self = shift;
	my($fs_coord) = @_;
	
	my $shoulder_state = $self->{tbl}[$fs_coord]{state};
	my $shoulder_len   = 0;
	foreach(my $i = $fs_coord; $i >= 0; $i--)
	{
		last if !defined $self->{tbl}[$i]{state} || $shoulder_state ne $self->{tbl}[$i]{state};
		$shoulder_len++;
	}
	
	return $shoulder_len;
}


sub _get_right_shoulder_len
{
	my $self = shift;
	my($fs_coord) = @_;

	my $shoulder_state = $self->{tbl}[$fs_coord+1]{state};
	my $shoulder_len   = 0;
	foreach(my $i = $fs_coord+1; $i < @{$self->{tbl}}; $i++)
	{
		last if $shoulder_state ne $self->{tbl}[$i]{state};
		$shoulder_len++;
	}

	return $shoulder_len;
}

###
# PUBLIC METHODS

###
# Arguments:
# 
# $tse -- TRANS_START_EXCEPT
# 
###
# %opts
#   |
#   |--->{fn_body}
# 
sub calculate_fs_scores
{
	my($self) = shift;
	my($mod_fn, $hmmdef_tmpl, $fs_prob, $tse, %opts) = @_;
	
	my($fn_body) = $opts{fn_body};
	die "Please specify fn_body or implement tmp_fn_body!!!" if !defined $fn_body;
	
	# Generate .hmm_def with and without possible transition between coding states (possibility of FS)
	my $mod = GeneTack::GeneMark::Mod->new($mod_fn);
	my $usual_hmmdef_fn = "$fn_body.usual.hmm_def";
	my $wo_fs_hmmdef_fn = "$fn_body.wo_fs.hmmdef";
	GeneTack::FSMark::create_hmm_def_file($hmmdef_tmpl, $mod, $usual_hmmdef_fn, fs_prob => $fs_prob, tse => $tse);
	GeneTack::FSMark::create_hmm_def_file($hmmdef_tmpl, $mod, $wo_fs_hmmdef_fn,
		fs_prob           => 0,
		trans_cod_deadend => 2*$fs_prob,           # All fs probability we redirect to deadend
		tse               => $tse,
	);
	
	my $seq = $self->get_seq;
	foreach my $fs ( @{$self->get_fs_info} )
	{
#		warn "Calculating score for frameshift at relative position $fs->{pos}...\n";
		$fs->{gene_nc_len} = $fs->{left_gene_nc_shoulder_len}+$fs->{right_gene_nc_shoulder_len};
		my $gene_nc_seq = substr($seq, $fs->{'pos'}-$fs->{left_gene_nc_shoulder_len}+1, $fs->{gene_nc_len});
		
		my $fasta_seq_fn    = "$fn_body.fs_score.$fs->{pos}.fasta";
		my $usual_fsmark_fn = "$fn_body.fs_score.$fs->{pos}.usual_hmmdef.fsmark";
		my $wo_fs_fsmark_fn = "$fn_body.fs_score.$fs->{pos}.wo_fs_hmmdef.fsmark";
		write_fasta($fasta_seq_fn, [{seq=>$gene_nc_seq, seqname=>"Seq for FS at position $fs->{pos}"}], width => 100);
		
		GeneTack::FSMark::run_fsmark($usual_hmmdef_fn, $fasta_seq_fn, out_fn => $usual_fsmark_fn);
		my $usual_fsm = GeneTack::FSMark->new($usual_fsmark_fn);
		
		if( $usual_fsm->num_fs == 0 )
		{
			# If a frameshift didn't appear in a shorter sequence we assume that this is a weak frameshift and assign score -1
			$fs->{fs_path_score}    = -1;
			$fs->{wo_fs_path_score} = -1;
			$fs->{score}            = -1;
			$fs->{lg_fs_score}      = -100;
		}
		else
		{
			# Apply .wo_fs.hmm_def that doesn't allow frameshifts
			GeneTack::FSMark::run_fsmark($wo_fs_hmmdef_fn, $fasta_seq_fn, out_fn => $wo_fs_fsmark_fn);
			my $wo_fs_fsm = GeneTack::FSMark->new($wo_fs_fsmark_fn);
			
			$fs->{fs_path_score}    = $usual_fsm->get_best_path_score;
			$fs->{wo_fs_path_score} = $wo_fs_fsm->get_best_path_score;
			$fs->{score}            = $fs->{wo_fs_path_score}/$fs->{fs_path_score};  # Both number are negative this is why we reverese the division
			$fs->{lg_fs_score}      = $fs->{score} > 1 ? log($fs->{score}-1.0)/log(10) : -100;
		}
		unlink($fasta_seq_fn, $usual_fsmark_fn, $wo_fs_fsmark_fn);
	}
	
	unlink($usual_hmmdef_fn, $wo_fs_hmmdef_fn);
}

###
# DESCRIPTION
#     Counts number of states in $self->{tbl} that are listed in $states
#
###
# Arguments:
#     $states -- reference to array with state names
#
sub cout_states
{
	my $self = shift;
	my($states) = @_;
	my $hash = array2hash( @$states );
	my $count = scalar( grep { $hash->{$_->{state}} } @{$self->{tbl}} );
	return $count;
}

sub get_state_at_pos
{
	my $self = shift;
	my($pos_i) = @_;
	return $pos_i < @{$self->{tbl}} ? $self->{tbl}[$pos_i]{state} : undef;
}

sub get_emiss_name_at_pos
{
	my $self = shift;
	my($pos_i) = @_;
	return $pos_i < @{$self->{tbl}} ? $self->{tbl}[$pos_i]{emiss_name} : undef;
}

sub get_best_path_score
{
	my $self = shift;
	return $self->{tbl}[-1]{best_score};
}

sub get_seq
{
	my $self = shift;
	return join('', map { $_->{let} } @{$self->{tbl}});
}

sub get_seq_letter_at_pos
{
	my $self = shift;
	my($pos_i) = @_;
	return $pos_i < @{$self->{tbl}} ? $self->{tbl}[$pos_i]{let} : undef;
}

###
# For a given fs returns sequence from both sides of the fs, i.e. sequence which corresponds to 
# states 111111111111111122222222222222222222 for instance.
# The sequence always starts with the zero-nucleotide, i.e. in right protein-coding frame
sub get_flanked_seqs_for_fs
{
	my $self = shift;
	my($fs) = @_;
	
	my $left_seq  = '';
	for(my $i = $fs->{'pos'}; ; --$i)
	{
		if(!defined $self->{tbl}[$i]{state} || $self->{tbl}[$i]{state} ne $fs->{from_state} || $i<0)  # We have reached the end of left sequence
		{
			$i++;   # One step back
			
			# Refine left coordinate -- find 'zero-nucleotide'
			$i++ while $self->{tbl}[$i]{emiss_name} ne '0';
			$left_seq = join '', map { $self->{tbl}[$_]{let} } $i..$fs->{'pos'};
			last;
		}
	}
	
	# Find start for right sequence
	my $r_start = $fs->{'pos'}+1;
	$r_start++ while $r_start < $self->seq_len && $self->{tbl}[$r_start]{emiss_name} ne '0' ;
	
	my $right_seq = '';
	for(my $i = $r_start; ; ++$i)
	{
		if($i>=$self->seq_len || $self->{tbl}[$i]{state} ne $fs->{to_state})  # We have reached the end of right sequence
		{
			$i--;   # One step back
			
			$i-- while $self->{tbl}[$i]{emiss_name} ne '0';
			$right_seq = join '', map { $self->{tbl}[$_]{let} } $r_start..($i-1);
			last;
		}
	}
	
	return($left_seq, $right_seq);
}

sub translate_fs
{
    my $self = shift;
	my($fs_coord,%opts) = @_;
	
	$opts{min_aa_shoulder} = 10 if !defined $opts{min_aa_shoulder};
	
	my $report = {};
	$self->get_protein_seq_for_fs($fs_coord, $opts{min_aa_shoulder}, $report);
	
	return $report;
}

sub get_protein_seq_for_fs
{
	my $self = shift;
	my($fs_coord, $min_aa_shoulder, $report) = @_;
	$min_aa_shoulder = 0 if !defined($min_aa_shoulder);
	
	my $info = $self->get_seq_info_for_fs_at_pos($fs_coord);
	
	if( defined $report )
	{
		foreach my $k (keys(%$info))
		{
			$report->{$k} = $info->{$k};
		}
	}
	
	my($left_len, $right_len) = (length($info->{left_prot}), length($info->{right_prot}));
	if( $left_len < $min_aa_shoulder || $left_len < $min_aa_shoulder )
	{
		if( defined $report )
		{
			push(@{$report->{notes_arr}},
				"Left/right ($left_len/$right_len) prot shoulder < threshold len");
			$report->{notes} = join "; ", @{$report->{notes_arr}};
		}
		return undef;
	}
	
	return $info->{prot_seq};
}

###
# $report  --  reference to hash. Following keys will be added/modified:
#    |
#    |--->{prot_seq}
#    |--->{gene_seq}     --  initial (not corrected) fsgene seq
#    |--->{left_dna}     --  left part of the gene
#    |--->{right_dna}    --  right part of the gene
#    |--->{left_prot}
#    |--->{right_prot}
#    |--->{prot_fs_coord}
#    |--->{notes_arr}
#    |--->{notes}
#
sub get_seq_info_for_fs_at_pos
{
	my $self = shift;
	my($fs_coord) = @_;
	
	my $report = {notes_arr => []};
	my $fs_info = $self->get_info_for_fs_at_pos($fs_coord);
	if( !$fs_info )
	{
		push @{$report->{notes_arr}}, "Can't get fs_info";
		$report->{notes} = join "; ", @{$report->{notes_arr}};
		return undef;
	}

	my($left_dna, $right_dna) = $self->get_flanked_seqs_for_fs( $fs_info );
	if( !$left_dna && !$right_dna )
	{
		push @{$report->{notes_arr}}, "Can't get left/right dna";
		$report->{notes} = join "; ", @{$report->{notes_arr}};
		return undef;
	}
	
	$report->{left_dna}  = $left_dna;
	$report->{right_dna} = $right_dna;
	$report->{left_prot}     = translate_dna_to_prot($left_dna);
	$report->{right_prot}    = translate_dna_to_prot($right_dna);
	$report->{prot_fs_coord} = length($report->{left_prot});
	$report->{prot_seq}      = lc($report->{left_prot}).uc($report->{right_prot});
	$report->{prot_seq_full} = $report->{prot_seq};
	$report->{gene_seq}      = lc($left_dna).uc($right_dna);
	$report->{notes}         = join "; ", @{$report->{notes_arr}};
	
	return $report;
}

sub get_info_for_fs_at_pos
{
       my $self = shift;
       my($fs_pos) = @_;

       foreach my $fs ( @{$self->{fs}} )
       {
               return $fs if $fs->{'pos'} == $fs_pos;
       }

       warn "No frameshift found at position $fs_pos. Existing fs:".join(',',map{$_->{'pos'}}@{$self->{fs}});
       return undef;
}

###
#    $fs     --  reference to array 
#     fs[i]
#      |
#      |--->{from_state}  --  
#      |--->{to_state}    --  
#      |--->{pos}         --  last position of sequence in state 'from_state'
# 
sub get_fs_info
{
	my $self = shift;
	return $self->{fs};
}


# First letter of the init_gene_seq corresponds to 1st codon position in the init_frame
# and the last letter corresponds to the 3rd codon position in the fs-frame
sub get_init_gene_seq_for_fs
{
	my $self = shift;
	my($fs_coord) = @_;
	
	my $fs = $self->get_info_for_fs_at_pos($fs_coord);	
	warn("get_init_gene_seq_for_fs: can't get info for fs_coord=$fs_coord!!!") and return undef if !$fs;
	
	my $left_seq  = '';
	for(my $i = $fs->{'pos'}; ; --$i)
	{
		if(!defined $self->{tbl}[$i]{state} || $self->{tbl}[$i]{state} ne $fs->{from_state} || $i<0)  # We have reached the end of left sequence
		{
			# We want start codon as well
			if($self->{tbl}[$i]{state} =~ /start/i){
				$i-=2
			} else {
				$i++;  # In case it is 0 emiss name in the other frame
				$i++ while $self->{tbl}[$i]{emiss_name} !~ /^0/; # start of the gene should be start of the codon
			}
			
			# Refine left coordinate -- find 'zero-nucleotide'
			$left_seq = join '', map { $self->{tbl}[$_]{let} } $i..$fs->{'pos'};
			last;
		}
	}
	
	# Find start for right sequence
	my $r_start   = $fs->{'pos'}+1;
	my $right_seq = '';
	for(my $i = $r_start; ; ++$i)
	{
		if($i>=$self->seq_len || $self->{tbl}[$i]{state} ne $fs->{to_state})  # We have reached the end of right sequence
		{
			if($self->{tbl}[$i]{state} !~ /stop/)
			{
				$i--;  # In case it is 2 emiss name in the other frame
				$i-- while $self->{tbl}[$i]{emiss_name} !~ /^2/; # last letter should be last position of codon
			}
			
			$right_seq = join '', map { $self->{tbl}[$_]{let} } $r_start..($i);
			last;
		}
	}
	
	return lc($left_seq).uc($right_seq);
}

###
# The function searches for a start state in the upstream region from the given frameshift
# Returns start state coordinate (-2 nt in upstream direction) or start of the fragment if no start state found
# 
sub get_gene_start_for_fs
{
	my $self = shift;
	my($fs_coord) = @_;
	
	foreach(my $start_coord = $fs_coord; $start_coord >= 0; $start_coord--)
	{
		last if !defined $self->{tbl}[$start_coord]{state};  # Start of fragment
		return($start_coord-2) if $self->{tbl}[$start_coord]{state} =~ /start/;
	}
	
	return 0;  # Return start of the framgment coordinate if no start state found

}

###
# The function searches for a stop state in the downstream region from the given frameshift
# Returns stop state coordinate or end of the fragment if no stop state found
# 
sub get_gene_end_for_fs
{
	my $self = shift;
	my($fs_coord) = @_;
	
	foreach(my $stop_coord = $fs_coord+1; $stop_coord < @{$self->{tbl}}; $stop_coord++)
	{
		return $stop_coord if $self->{tbl}[$stop_coord]{state} =~ /stop/;
	}
	
	return scalar(@{$self->{tbl}})-1;    # Return end of the framgment coordinate if no stop state found
}

sub num_fs
{
	my $self = shift;
	return scalar( @{$self->{fs}} );
}

sub has_fs
{
	my $self = shift;
	return scalar( @{$self->{fs}} );
}

sub seq_len
{
	my $self = shift;
	return scalar( @{$self->{tbl}} );
}

###
# SUBROUTINES

sub run_fsmark
{
	my($hmm_def_fn, $seq_fn, %opts) = @_;
	$opts{prg_path} ||= $PRG_PATH;
	
	if( !$opts{out_fn} )   # create out_file name if it's not provided
	{
		my($tmp_fh, $out_fn) = tempfile('XXXXXXXX', SUFFIX => '.fsmark');
		close $tmp_fh;
		$opts{out_fn} = $out_fn;
	}
	
	system($opts{prg_path}.' '.$hmm_def_fn.' '.$seq_fn.' > '.$opts{out_fn});
	
	return $opts{out_fn};
}

###
# Arguments:
#     $tp_fn  --  template file name
#     $mod    --  GeneTack::GeneMark::Mod object
#     $out_fn --  name of file that will be created (usually with .hmm_def extention)
#     %opts
#       |
#       |--->{fs_prob}            --  corresponds to X (see below)
#       |--->{tse}                --  (trans_start_except) exception transition probability to start codon
#       |--->{trans_cod_deadend}  --  transition probability from coding state to deadend
# 
# Notation:
#     v -- Usual transition probability from coding state to stop codon (and from stop state only one way to n/c state exists)
#     w -- Transition probability from N/C to itself
#     x -- [Nov 21, 2008] transition to another coding state:   0   -> 1  -- probability of frame shift, actually
#     y -- [Nov 21, 2008] transition to the same coding state:  0   -> 0
#     z -- [Nov 21, 2008] transition to the same overlap state: 0+1 -> 0+1
# 
sub create_hmm_def_file
{
	my($tp_fn, $mod, $out_fn, %opts) = @_;
	my $tse     = $opts{tse};
	my $x       = $opts{fs_prob};
	my $deadend = $opts{trans_cod_deadend} || 0;
	confess('Some required options are not specified: '.Data::Dumper->Dump([\%opts], ['$opts'])) if !defined $x || !defined $tse;
	
	my $v = 0;
	my $y = 1 - 2*$x - $v - $deadend;
	my $w = 1;
	my $z = $y;
	
	my $tp = HTML::Template->new(
		filename          => $tp_fn,
		shared_cache      => 0,
		die_on_bad_params => 0,
		loop_context_vars => 1,
	);
	
	$tp->param(
		head_comment_line => 'This file was automatically generated by '.File::Spec->splitpath($0).' from '.$mod->init_fn().' on '.localtime(),
		file_info_ver     => 0.01,
		hmm_model_name    => 'FSMark_HMM',
		hmm_model_descr   => $mod->get_order().' order model of coding states',
		emission_list_arr => get_emission_list_arr($mod, $opts{cod_model}),
		                                                            ### Transition probabilities: ###
		trans_cod1_cod1   => to_num($y),                            # From coding state to itself
		trans_cod1_cod2   => to_num($x),                            # Between coding states
		trans_cod_ovlp    => to_num((1-2*$x-$y-$v-$deadend)/2),     # From coding state to overlap state
		trans_ovlp_ovlp   => to_num($z),                            # From overlap state to itself
		trans_ovlp_stop   => to_num(1-$z),                          # From overlap state to stop state
		trans_cod_nc      => to_num($v),                            # From coding state to stop codon
		trans_nc_nc       => to_num($w),                            # From N/C to itself
		trans_nc_cod      => to_num((1-$w)/3),                      # From N/C to start codon
		trans_start_except=> $tse,                                  # Exception transition probability to start codon
		trans_cod_deadend => $deadend,
	);
	
	open(my $out_fh, '>', $out_fn) or confess("Can't open file '$out_fn' to write: $!");
	print $out_fh $tp->output();
	close $out_fh;
}

sub get_emission_list_arr
{
	my($mod, $cod_model) = @_;
	$cod_model = 'COD1' if !$cod_model;
	
	my $freqs_order = $ORDER_COD_FREQS->{ $mod->get_order() };
	confess('Order '.$mod->get_order().' is not supported') if !defined $freqs_order;
	
	my %freqs_0     = %{$mod->{freqs}{$cod_model}{$freqs_order->[0]}};  # Freqs_0  --  correspond to emission of nucleotide in the 1st codon position
	my %freqs_1     = %{$mod->{freqs}{$cod_model}{$freqs_order->[1]}};  # Freqs_1  --  correspond to emission of nucleotide in the 2nd codon position
	my %freqs_2     = %{$mod->{freqs}{$cod_model}{$freqs_order->[2]}};  # Freqs_2  --  correspond to emission of nucleotide in the 3rd codon position
	my %f_0_wo_stop = remove_freqs_for_cods($STOP_CODONS, \%freqs_0);
	my %f_1_wo_stop = remove_freqs_for_cods($STOP_CODONS, \%freqs_1);
	my %f_2_wo_stop = remove_freqs_for_cods($STOP_CODONS, \%freqs_2);
	my %freqs_nc    = %{$mod->get_noncoding_freqs()};
	my %stop_freqs  = (TAG => 0.34, TAA => 0.33, TGA => 0.33);
	my %start_freqs = (ATG => 0.80, GTG => 0.10, TTG => 0.10);
	my %zero_freqs  = (XXX => 1);
	
	my $tmpl = "\t\t\t%s = %f";
	return [
		{
			id        => 0,
			name      => 0,
			comment   => "Frequencies of triplets 1 nt shifted from coding frame ($cod_model), but in condition probabilities they actually emit nucleotide in 1st codon position",
			probs_str => join("\n", map { sprintf($tmpl, $_, $freqs_0{$_}) } sort keys %freqs_0),
		},
		{
			id        => 1,
			name      => 1,
			comment   => "Frequencies of triplets 2 nt shifted from coding frame ($cod_model), but in condition probabilities they actually emit nucleotide in 2nd codon position",
			probs_str => join("\n", map { sprintf($tmpl, $_, $freqs_1{$_}) } sort keys %freqs_1),
		},
		{
			id        => 2,
			name      => 2,
			comment   => "Frequencies of triplets in coding frame ($cod_model), but in condition probabilities they emit nucleotide actually in 3rd codon position",
			probs_str => join("\n", map { sprintf($tmpl, $_, $freqs_2{$_}) } sort keys %freqs_2),
		},
		{
			id        => 3,
			name      => '0_1',
			probs_str => join("\n", map { sprintf($tmpl, $_, 0.5*$freqs_0{$_}+0.5*$freqs_1{$_}) } sort keys %freqs_0),
		},
		{
			id        => 4,
			name      => 'stop',
			probs_str => join("\n", map { sprintf($tmpl, $_, $stop_freqs{$_}) } sort keys %stop_freqs),
		},
		{
			id        => 5,
			name      => 'start',
			probs_str => join("\n", map { sprintf($tmpl, $_, $start_freqs{$_}) } sort keys %start_freqs),
		},
		{
			id        => 6,
			name      => 'nc',
			probs_str => join("\n", map { sprintf($tmpl, $_, $freqs_nc{$_}) } sort keys %freqs_nc),
		},
		{
			id        => 7,
			name      => 'zero',
			probs_str => join("\n", map { sprintf($tmpl, $_, $zero_freqs{$_}) } sort keys %zero_freqs),
		},
		{
			id        => 8,
			name      => '0_wo_stop',
			probs_str => join("\n", map { sprintf($tmpl, $_, $f_0_wo_stop{$_}) } sort keys %f_0_wo_stop),
		},
		{
			id        => 9,
			name      => '1_wo_stop',
			probs_str => join("\n", map { sprintf($tmpl, $_, $f_1_wo_stop{$_}) } sort keys %f_1_wo_stop),
		},
		{
			id        => 10,
			name      => '2_wo_stop',
			probs_str => join("\n", map { sprintf($tmpl, $_, $f_2_wo_stop{$_}) } sort keys %f_2_wo_stop),
		},
	];
}

sub remove_freqs_for_cods
{
	my($codons_to_remove, $freqs) = @_;
	
	# Copy only codons that are not in @$codons_to_remove
	my %res_freqs = ();
	my $total_removed_freq = 0;
	foreach my $cod ( keys %$freqs )
	{
		if( grep { $cod =~ /$_$/i } @$codons_to_remove )
		{
			$total_removed_freq += $freqs->{$cod}; # We are removing this freq. So let's remember how much freq we removed
		}
		else
		{
			$res_freqs{$cod} = $freqs->{$cod};
		}
	}
	
	# Assign all removed freq to non-existing codon -- so the sum is still 1
	my $empty_codon = 'X' x length( each %$freqs );
	$res_freqs{$empty_codon} = $total_removed_freq;
	
	return %res_freqs;
}

sub to_num
{
	return 0 if $_[0] == 0;
	
	# 5e-005 => 0.0000500000 => 0.00005
	(my $num = sprintf('%.10f', $_[0])) =~ s/\.?0*$//;
	
	# Remove side effect: "-0"  =>  "0"
	$num =~ s/^-[0\.]+$/0/;
	
	return $num;
}

sub sum
{
	my $res = 0;
	$res += $_ foreach @_;
	return $res;
}

sub array2hash
{
	my %hash = ();
	$hash{$_}++ foreach @_;
	return \%hash;
}

1;

