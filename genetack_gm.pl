#!/usr/bin/perl --

use strict;
use warnings;

my $VERSION = '1.30';

# $Id$

###
# Ivan Antonov (vanya.antonov@gmail.com)
#

$| = 1; # Turn off buffering

use FindBin;
use lib "$FindBin::Bin/lib";
$ENV{PATH} .= ":$FindBin::Bin";

use Data::Dumper;
$Data::Dumper::Deepcopy = 1;
use Getopt::Long;
use Cwd qw(abs_path);
use File::Spec;
use File::Path qw(remove_tree);
use File::Copy;

use GeneTack::Lib qw(create_tmp_dir write_fasta read_fasta get_gc_content);
use GeneTack::FSMark;
use GeneTack::GeneMark;
use GeneTack::GeneMark::Lst;
use GeneTack::GeneMark::Mod;

###
# CONSTANTS

# $FS_PROB -- transition probability between coding states.
# From meeting_log.doc on November 21, 2008: "Value of variable X depends on the sequencing technique used.
# This is one of the FSMark parameters. Based on information from experts from a sequencing lab we decided
# to set default value to 10-4 (see meeting log on September 18, 2008)."
my $NUM_THREADS         = 1;
my $FS_PROB             = 0.0006;
my $TSE                 = 0.0001;     # TRANS_START_EXCEPT -- Based of optimization done on April 21, 2009
my $TSE_HIGH_GC         = 0.001;      # TRANS_START_EXCEPT for high GC content genomes
my $COD_ORDER           = 3;          # Default value for --cod_order option
my $GM_ORDER            = 5;          # Default value for --gm_order option
my $GM_CDEC             = 300;        # Default value for --gm_cdec option
my $SAVE_SEQS           = '';         # Default value for --save_seqs option
my $SAVE_FILTER_REGS    = '';
my $SAVE_SEQS_INFO      = '';         # Default value for --save_seqs_info option
my $SAVE_FSGENE_SEQS    = '';         # Default value for --save_fsgene_seqs option
my $SAVE_FSPROT_SEQS    = '';         # Default value for --save_fsprot_seqs option
my $SAVE_GM_MOD         = '';
my $SAVE_FS_MOD         = '';
my $FS_COD_SHOULDER_LEN = 50;         # Default value for --fs_cod_shoulder_len option
my $MIN_OVERLAP_LEN     = 8;          # Default value for --min_overlap_len option
my $MIN_ADJ_LEN         = 1;          # Default value for --min_adj_len option
my $MAX_OVERLAP_RBS     = 2;          # Default value for --max_overlap_rbs option
my $MIN_LG_FS_SCORE     = -4.5;       # Default value for --min_lg_fs_score option
my $OUTPUT_FILTERED_FS  = 0;          # Default value for --output_filtered_fs option
my $IGNORE_VICINITY     = 50;         # Default value for --ignore_vicinity option
my $MAX_ADJACENT_GAP    = 100;        # Default value for --max_adjacent_gap option
my $MAX_ADJACENT_RBS    = 1;          # Default value for --max_adjacent_rbs option
my $ADJ_STOP_EXTRA_LEN  = 20;         # Default value for --adj_stop_extra_len option 
my $SAVE_FSMARK_FILES   = '';         # Default value for --save_fsmark_files option
my $CALC_FS_SCORES      = 0;          # Default value for --calc_fs_scores option
my $HIGH_GC_THR         = 0.65;       # Threshold for high GC content genomes
my $GM_MOD_FN           = '';
my $FS_MOD_FN           = '';
my $TMP_DIR             = '.';

###
# Parse input data
GetOptions(
	'num_threads=i'        => \$NUM_THREADS,
	'fs_prob=f'            => \$FS_PROB,
	'gm_order=i'           => \$GM_ORDER,
	'gm_cdec=i'            => \$GM_CDEC,
	'cod_order=i'          => \$COD_ORDER,
	'save_seqs=s'          => \$SAVE_SEQS,
	'save_filter_regs=s'   => \$SAVE_FILTER_REGS,
	'save_seqs_info=s'     => \$SAVE_SEQS_INFO,
	'save_fsgene_seqs=s'   => \$SAVE_FSGENE_SEQS,
	'save_fsprot_seqs=s'   => \$SAVE_FSPROT_SEQS,
	'save_gm_mod=s'        => \$SAVE_GM_MOD,
	'save_fs_mod=s'        => \$SAVE_FS_MOD,
	'fs_cod_shoulder_len=i'=> \$FS_COD_SHOULDER_LEN,
	'min_overlap_len=i'    => \$MIN_OVERLAP_LEN,
	'output_filtered_fs'   => \$OUTPUT_FILTERED_FS,
	'ignore_vicinity=i'    => \$IGNORE_VICINITY,
	'max_overlap_rbs=f'    => \$MAX_OVERLAP_RBS,
	'min_lg_fs_score=f'    => \$MIN_LG_FS_SCORE,
	'max_adjacent_gap=i'   => \$MAX_ADJACENT_GAP,
	'max_adjacent_rbs=f'   => \$MAX_ADJACENT_RBS,
	'adj_stop_extra_len=i' => \$ADJ_STOP_EXTRA_LEN,
	'save_fsmark_files=s'  => \$SAVE_FSMARK_FILES,
	'calc_fs_scores'       => \$CALC_FS_SCORES,
	'high_gc_thr=f'        => \$HIGH_GC_THR,
	'gm_mod_fn=s'          => \$GM_MOD_FN,
	'fs_mod_fn=s'          => \$FS_MOD_FN,
	'tmp_dir=s'            => \$TMP_DIR,
) || die usage();

die usage() if @ARGV!=1;
die usage("Value of --cod_order $COD_ORDER is not supported") if !exists $GeneTack::FSMark::ORDER_COD_FREQS->{$COD_ORDER};

###
my $START_TIME = time;

run(
	seq_fn             => $ARGV[0],
	num_threads        => $NUM_THREADS,
	fs_prob            => $FS_PROB,
	tse                => $TSE,
	tse_high_gc        => $TSE_HIGH_GC,
	cod_order          => $COD_ORDER,
	gm_order           => $GM_ORDER,
	gm_cdec            => $GM_CDEC,
	save_seqs          => $SAVE_SEQS,
	save_filter_regs   => $SAVE_FILTER_REGS,
	save_seqs_info     => $SAVE_SEQS_INFO,
	save_fsgene_seqs   => $SAVE_FSGENE_SEQS,
	save_fsprot_seqs   => $SAVE_FSPROT_SEQS,
	save_gm_mod        => $SAVE_GM_MOD,
	save_fs_mod        => $SAVE_FS_MOD,
	output_filtered_fs => $OUTPUT_FILTERED_FS,
	ignore_vicinity    => $IGNORE_VICINITY,
	save_fsmark_files  => $SAVE_FSMARK_FILES,
	calc_fs_scores     => $CALC_FS_SCORES,
	high_gc_thr        => $HIGH_GC_THR,
	gm_mod_fn          => $GM_MOD_FN ? abs_path($GM_MOD_FN) : undef,
	fs_mod_fn          => $FS_MOD_FN ? abs_path($FS_MOD_FN) : undef,
	fn_body            => 'GENOME',
	genemark_path      => '',
	hmmdef_tmpl        => "$FindBin::Bin/hmm_def_files/genetack_gm.hmm_def",
	tmp_root           => $TMP_DIR,
	filter_params      => {
		min_overlap_len     => $MIN_OVERLAP_LEN,
		min_adj_len         => $MIN_ADJ_LEN,
		fs_cod_shoulder_len => $FS_COD_SHOULDER_LEN,
		max_overlap_rbs     => $MAX_OVERLAP_RBS,
		max_adjacent_gap    => $MAX_ADJACENT_GAP,
		max_adjacent_rbs    => $MAX_ADJACENT_RBS,
		adj_stop_extra_len  => $ADJ_STOP_EXTRA_LEN,
		min_lg_fs_score     => $MIN_LG_FS_SCORE,
	},
);

warn "\nElapsed time: ".(time-$START_TIME)." sec\n";
###

###
# SUBROUTINES
sub run
{
	my %opts = @_;

	$opts{tmp_dir} = create_tmp_dir(root_dir => $opts{tmp_root}, prefix => '__genetack_');

	# Copy genome to the tmp dir, so that all the files will be created inside the same folder
	copy($opts{seq_fn}, "$opts{tmp_dir}/$opts{fn_body}.fasta");
	$opts{seq_fn} = abs_path("$opts{tmp_dir}/$opts{fn_body}.fasta");

	# Read genome sequence and determine whether it is high GC genome
	my $Seqs     = read_fasta( $opts{seq_fn} );
	my $totSeq = join '', values %{ $Seqs };
	my $high_gc = get_gc_content( \$totSeq ) > $opts{high_gc_thr} ? 1 : 0;	#???

	# Run GeneMark. Create .gdata file for the high GC genomes only.
	my($gm_mod_fn, $fs_mod_fn, $lst_fn, $gdata_fn) = GeneTack::GeneMark::run_gm($opts{seq_fn}, $opts{cod_order}, $opts{gm_order}, $opts{gm_cdec},
		run_dir       => $opts{tmp_dir},
		genemark_path => $opts{genemark_path},
		gdata         => $high_gc,
		fs_mod_fn     => $opts{fs_mod_fn},
		gm_mod_fn     => $opts{gm_mod_fn},
	);

	my $tse_prob = $high_gc ? $opts{tse_high_gc} : $opts{tse};

	my $ifs = 1; # FrameShift id

	for my $refID ( keys %{ $Seqs } )
	{
		my $seq = $Seqs->{ $refID };

		# Prepare TODO for FSMark -- genome chunks which will be used as input to run FSMark
		my $lst  = GeneTack::GeneMark::Lst->new($lst_fn, $refID, gdata_fn => $gdata_fn);
		my $todo = get_fsmark_todo( $lst, $seq, defined( $gdata_fn ));
		my $res  = run_fsmark($todo, $seq, $opts{hmmdef_tmpl}, $fs_mod_fn, $opts{fn_body}, $opts{fs_prob}, $tse_prob, %opts);

		my $all_filters = apply_filters($lst, $seq, $res, $opts{filter_params}, $opts{ignore_vicinity});

		# Output
		output($res, $opts{output_filtered_fs}, $opts{calc_fs_scores}, $refID, \$ifs);

		save_seqs($todo, $opts{save_seqs})                                 if $opts{save_seqs};
		save_seqs_info($todo, $opts{save_seqs_info})                       if $opts{save_seqs_info};
		save_seqs_fasta($res, 'prot_seq', $opts{save_fsprot_seqs}, %opts)  if $opts{save_fsprot_seqs};
		save_seqs_fasta($res, 'gene_seq', $opts{save_fsgene_seqs}, %opts)  if $opts{save_fsgene_seqs};
		save_filter_regs($opts{save_filter_regs}, $all_filters)            if $opts{save_filter_regs};
		copy($gm_mod_fn, $opts{save_gm_mod})                               if $opts{save_gm_mod};
		copy($fs_mod_fn, $opts{save_fs_mod})                               if $opts{save_fs_mod};
	}

	remove_tree($opts{tmp_dir}) or warn "Couldn't remove tmp dir $opts{tmp_dir}";
}

###
# Input Parameters:
#    $data      --   reference to array
#     data[i]
#       |
#       |--->{strand}
#       |--->{left}
#       |--->{right}
#       |--->{fsmark_fn}
#       |--->{fs}
#             fs[i]
#             |
#         =====================================
#          All keys from
#          GeneTack::FSMark::get_fs_info
#         =====================================
#             |
#             |--->{_coord}
#             |--->{_g_left}
#             |--->{_g_right}
#             |--->{_filter}
#
sub output
{
	my( $data, $filtered_fs, $fs_scores, $refID, $ifs ) = @_;

	my $tmpl = "%-6s\t%-15s\t%-15s\t%-15s\t%-7s\t%-6s\t%-15s\t%-15s\t%-15s\t%-15s";
	$tmpl .= "\t%-12s\t%s\t%s\t%s\t%s" if $fs_scores;
	$tmpl .= "\t%-15s" if $filtered_fs;
	$tmpl .= "\n";

	if( $$ifs < 2 ){
		my @head = qw(Num Seq_ID FS_coord FS_coord_adj FS_type Strand Gene_Left Gene_Right Fragment_Left Fragment_Right);
		push @head, qw(Gene_NC_Len FS_Path_Score Wo_FS_Path_score FS_Score LG_FS_Score) if $fs_scores;
		push @head, 'Filter' if $filtered_fs;
		printf($tmpl, @head);
	}

	@$data = sort { $a->{left} <=> $b->{left} } @$data;

	foreach my $d ( @$data )
	{
		foreach my $fs ( @{$d->{fs}} )
		{
			my @values = ( $$ifs, $refID, $fs->{_coord}, $fs->{_coord_adj}, $fs->{type}, $d->{strand}, $fs->{_g_left}, $fs->{_g_right}, $d->{left}, $d->{right});
			push(@values, $fs->{gene_nc_len}, $fs->{fs_path_score}, $fs->{wo_fs_path_score}, $fs->{score}, $fs->{lg_fs_score}) if $fs_scores;
			if( $filtered_fs )
			{
				push @values, $fs->{_filter} || '-';
			}
			else
			{
				next if $fs->{_filter};
			}
			printf($tmpl, @values);
			++$$ifs;
		}
	}
}

sub apply_filters
{
	my($lst, $genome, $data, $filter_params, $ignore_vicinity) = @_;

	my $all_filters = get_all_filters($lst, $genome, $filter_params, $ignore_vicinity);

	# hash of references to frame shifts by they coordinate
	my %fs_coords = ();
	foreach my $item ( @$data )
	{
		$fs_coords{$_->{_coord}} = $_ foreach @{$item->{fs}};
	}

	foreach my $filter ( @$all_filters )
	{
		foreach my $reg ( @{$filter->{ignore_regs}} )
		{
			foreach my $coord ( grep { $reg->{left} <= $_ && $_ <= $reg->{right} } keys %fs_coords )   # Find FS inside the ignore region
			{
				$fs_coords{$coord}{_filter} = $filter->{name} if !$fs_coords{$coord}{_filter};         # Filter the FS if it is not filtered yet
			}
		}
	}

	my $cod_shoulder = $filter_params->{fs_cod_shoulder_len};
	foreach my $item ( @$data )
	{
		foreach my $fs ( @{$item->{fs}} )
		{
			if( $fs->{left_cod_shoulder_len}<$cod_shoulder || $fs->{right_cod_shoulder_len}<$cod_shoulder)
			{
				$fs->{_filter} = 'short_fs_shoulder' if !$fs->{_filter};
			}
			elsif ( defined $fs->{lg_fs_score} && $fs->{lg_fs_score} < $filter_params->{min_lg_fs_score} )
			{
				$fs->{_filter} = 'small_lg_fs_score' if !$fs->{_filter};
			}
		}
	}

	return $all_filters;
}

###
# Returns:
#     $all_filters      --   array of hashes
#      all_filters[0]
#            |
#            |--->{name}             --  filter name
#            |--->{ignore_regs}      --  reference to array of hashes
#                  ignore_regs[i]
#                     |
#                     |--->{left}
#                     |--->{right}
# 
sub get_all_filters
{
	my($lst, $genome, $filter_params, $ignore_vicinity) = @_;

	return [
		{
			name        => 'ovlp_len',
			ignore_regs => get_ovlp_len_ignore_regs($lst, $filter_params->{min_overlap_len}, $ignore_vicinity),
		},
		{
			name        => 'adj_len',
			ignore_regs => get_adj_len_ignore_regs($lst, $filter_params->{min_adj_len}, $ignore_vicinity),
		},
		{
			name        => 'ovlp_rbs',
			ignore_regs => get_ovlp_rbs_ignore_regs($lst, $filter_params->{max_overlap_rbs}, $ignore_vicinity),
		},
		{
			name        => 'adj_rbs',
			ignore_regs => get_adj_rbs_ignore_regs($lst, $filter_params->{max_adjacent_gap}, $filter_params->{max_adjacent_rbs}, $ignore_vicinity),
		},
		{
			name        => 'adj_stop',
			ignore_regs => get_adj_stop_ignore_regs($lst, $genome, $filter_params->{max_adjacent_gap}, $filter_params->{adj_stop_extra_len}, $ignore_vicinity),
		},
	];
}

sub get_adj_stop_ignore_regs
{
	my($lst, $genome, $gap_limit, $extra_len, $ignore_vicinity) = @_;

	my $all_adjs    = $lst->get_all_adjacent_gene_pairs($gap_limit);
	my $ignore_regs = [];
	foreach my $adj ( @$all_adjs )
	{
		my $gap_seq = get_adj_gap_seq($genome, $adj, $extra_len);
		next if !$gap_seq;
		
		if( $gap_seq =~ /(TAG|TGA|TAA)(...)*(ATG|GTG|TTG)$/i ) 
		{
			my $gap_left  = $adj->{strand} eq '+' ? $adj->{up_gene}{right}  :  $adj->{down_gene}{right};
			my $gap_right = $adj->{strand} eq '+' ? $adj->{down_gene}{left} :  $adj->{up_gene}{left};
			warn "Something is wrong with coordinates: ".Data::Dumper->Dump([$adj], ['adj']) and next if $gap_left > $gap_right;   # Double check coordinates
			push @$ignore_regs, {
				left  => $gap_left  - $ignore_vicinity,
				right => $gap_right + $ignore_vicinity,
			};	
		}
	}

	return $ignore_regs;
}

sub get_adj_gap_seq
{
	my($genome, $adj, $extra_len) = @_;

	my $gap_seq;
	if( $adj->{strand} eq '+' )
	{
		my $start = $adj->{down_gene}{left}-$adj->{gap_len}-$extra_len;
		$gap_seq  = substr($genome, $start, $adj->{gap_len}+$extra_len+2);
	}
	else
	{
		my $start = $adj->{down_gene}{right}-3;
		$gap_seq  = substr($genome, $start, $adj->{gap_len}+$extra_len+3);
		$gap_seq  = revcomp($gap_seq);
	}

	if($gap_seq !~ /([AGT]TG)$/i)
	{
		warn "Something is wrong with seq = '$gap_seq': ".Data::Dumper->Dump([$adj], ['$adj']);
		return undef;
	}

	return $gap_seq;
}

sub get_adj_rbs_ignore_regs
{
	my($lst, $gap_limit, $max_rbs, $ignore_vicinity) = @_;

	my $all_adjs    = $lst->get_all_adjacent_gene_pairs($gap_limit);
	my $ignore_regs = [];
	foreach my $adj ( @$all_adjs )
	{
		if( $adj->{down_gene}{rbs_score} > $max_rbs )
		{
			my $gap_left  = $adj->{strand} eq '+' ? $adj->{up_gene}{right}  :  $adj->{down_gene}{right};
			my $gap_right = $adj->{strand} eq '+' ? $adj->{down_gene}{left} :  $adj->{up_gene}{left};
			die "Something wrong: ".Data::Dumper->Dump([$adj], ['adj']) if $gap_left > $gap_right;   # Double check coordinates
			push @$ignore_regs, {
				left  => $gap_left  - $ignore_vicinity,
				right => $gap_right + $ignore_vicinity,
			};
		}
	}
	return $ignore_regs;
}

sub get_ovlp_rbs_ignore_regs
{
	my($lst, $max_overlap_rbs, $ignore_vicinity) = @_;

	my $ignore_regs = [];
	foreach my $o ( @{$lst->get_pair_overlaps_info} )
	{
		my($g1, $g2) = ($lst->get_gene($o->{genes_i}[0]), $lst->get_gene($o->{genes_i}[1]));
		my $down_rbs_score = $o->{strand} eq '+' ? $g2->{rbs_score} : $g1->{rbs_score};
		if( $down_rbs_score > $max_overlap_rbs )
		{
			push @$ignore_regs, {
				left  => $g2->{left}  - $ignore_vicinity,
				right => $g1->{right} + $ignore_vicinity,
			};		
		}
	}

	return $ignore_regs;
}

sub get_ovlp_len_ignore_regs
{
	my($lst, $min_overlap_len, $ignore_vicinity) = @_;

	my $ignore_regs = [];
	foreach my $o ( @{$lst->get_pair_overlaps_info} )
	{
		if ( $o->{len} < $min_overlap_len )
		{
			my($g1, $g2) = ($lst->get_gene($o->{genes_i}[0]), $lst->get_gene($o->{genes_i}[1]));
			die "Something wrong: ".Data::Dumper->Dump([$o], ['ovlp']) if $g1->{right} < $g2->{left};   # Check whether it is really overlap
			push @$ignore_regs, {
				left  => $g2->{left}  - $ignore_vicinity,
				right => $g1->{right} + $ignore_vicinity,
			};
		}
	}

	return $ignore_regs;
}

sub get_adj_len_ignore_regs
{
	my($lst, $min_adj_len, $ignore_vicinity) = @_;

	my $ignore_regs = [];
	foreach my $adj ( @{$lst->get_all_adjacent_gene_pairs($min_adj_len+1)} )
	{
		my $gap_left  = $adj->{strand} eq '+' ? $adj->{up_gene}{right}  : $adj->{down_gene}{right};
		my $gap_right = $adj->{strand} eq '+' ? $adj->{down_gene}{left} : $adj->{up_gene}{left};
		push @$ignore_regs, {
			left  => $gap_left  - $ignore_vicinity,
			right => $gap_right + $ignore_vicinity,
		};
	}

	return $ignore_regs;
}

###
# Arguments:
# tse -- TRANS_START_EXCEPT
#
###
# %opts
#   |
#   |--->{fsmark_dir}
# 
sub run_fsmark
{
	my($todo, $genome, $hmmdef_tmpl, $mod_fn, $fn_body, $fs_prob, $tse, %opts) = @_;

	# Generate hmm_def file for FSMark
	my $mod       = GeneTack::GeneMark::Mod->new($mod_fn);
	my $hmmdef_fn = "$opts{tmp_dir}/$fn_body.hmm_def";
	GeneTack::FSMark::create_hmm_def_file($hmmdef_tmpl, $mod, $hmmdef_fn, fs_prob => $fs_prob, tse => $tse);

	my $chunk_dir = "$opts{tmp_dir}/seq_chunk";
	mkdir $chunk_dir if !-d $chunk_dir;

	my $fsmark_dir = "$opts{tmp_dir}/fsmark_files";
	mkdir $fsmark_dir if !-d $fsmark_dir;

	mkdir $opts{save_fsmark_files} if !-d $opts{save_fsmark_files};

	my $fs_score_dir = undef;
	if( $opts{calc_fs_scores} )
	{
		### Temporary
		$fs_score_dir = "$opts{tmp_dir}/fs_score/";
		mkdir $fs_score_dir if !-d $fs_score_dir;
	}

	my($i, $res) = (1, []);
	if( $opts{num_threads} == 1 )
	{
		foreach my $item ( @$todo )
		{
			print STDERR "\rGeneTack is running for fragment ".($i++)."/".@$todo.": $item->{left} .. $item->{right}           ";
			
			my $fsmark_fn = _generate_fsmark_file($item->{_seq}, $item->{left}, $hmmdef_fn, $chunk_dir, $fsmark_dir);
			my $fs_hash   = _get_all_fs_from_fsmark_file($fsmark_fn, $item, $mod_fn, $hmmdef_tmpl, $fs_prob, $tse, $fs_score_dir, %opts);
			push(@$res, $fs_hash) if $fs_hash;
			
			# Save the first 4 columns from the fs-mark file if requested
			`cut -f1-4 $fsmark_fn  >  $opts{save_fsmark_files}/$item->{left}.fsmark` if $opts{save_fsmark_files};
		}
	}
	else
	{
		# https://perlmaven.com/speed-up-calculation-by-running-in-parallel
		require Parallel::ForkManager;

		my $pm = Parallel::ForkManager->new( $opts{num_threads} );
		$pm->run_on_finish( sub {
			my ($pid, $exit_code, $ident, $exit_signal, $core_dump, $data) = @_;

			my $fs_hash = _get_all_fs_from_fsmark_file($data->{fsmark_fn}, $data->{item}, $mod_fn, $hmmdef_tmpl, $fs_prob, $tse, $fs_score_dir, %opts);
			push(@$res, $fs_hash) if $fs_hash;

			# Save the first 4 columns from the fs-mark file if requested
			`cut -f1-4 $data->{fsmark_fn}  >  $opts{save_fsmark_files}/$data->{item}->{left}.fsmark` if $opts{save_fsmark_files};
		});

		foreach my $item ( @$todo )
		{
			print STDERR "\rGeneTack is running for fragment ".($i++)."/".@$todo.": $item->{left} .. $item->{right}           ";
			my $pid = $pm->start and next;
			my $fsmark_fn = _generate_fsmark_file($item->{_seq}, $item->{left}, $hmmdef_fn, $chunk_dir, $fsmark_dir);
			$pm->finish(0, {fsmark_fn => $fsmark_fn, item => $item});
		}
		$pm->wait_all_children;
	}
	print STDERR "\n";

	return $res;
}

sub _get_all_fs_from_fsmark_file
{
	my($fsmark_fn, $item, $mod_fn, $hmmdef_tmpl, $fs_prob, $tse, $fs_score_dir, %opts) = @_;

	my $fsm = GeneTack::FSMark->new( $fsmark_fn );

	$fsm->calculate_fs_scores($mod_fn, $hmmdef_tmpl, $fs_prob, $tse, fn_body => $fs_score_dir.$item->{left}) if $opts{calc_fs_scores};

	my $fs_hash = undef;
	if( $fsm->has_fs )
	{
		$fs_hash = {
			fs     => $fsm->get_fs_info,
			left   => $item->{left},
			right  => $item->{right},
			strand => $item->{strand},
			origin => $item->{origin},
			fsmark_fn => $fsmark_fn,
		};

		foreach my $fs ( @{$fs_hash->{fs}} )   # Calculate absolute coordintate
		{
			my $g_start       = $fsm->get_gene_start_for_fs($fs->{'pos'});
			my $g_end         = $fsm->get_gene_end_for_fs($fs->{'pos'})+1;
			$fs->{_coord}     = $item->{strand} eq '+' ? ($item->{left} + $fs->{'pos'})   : ($item->{right} - $fs->{'pos'});
			$fs->{_coord_adj} = $item->{strand} eq '+' ? ($item->{left} + $fs->{pos_adj}) : ($item->{right} - $fs->{pos_adj});
			$fs->{_g_left}    = $item->{strand} eq '+' ? ($item->{left} + $g_start)       : ($item->{right} - $g_end      );
			$fs->{_g_right}   = $item->{strand} eq '+' ? ($item->{left} + $g_end)         : ($item->{right} - $g_start    );
		}
	}

	return $fs_hash;
}

sub _generate_fsmark_file
{
	my($seq, $seqname, $hmmdef_fn, $chunk_dir, $fsmark_dir) = @_;
	
	my $tmp_seq_fn = "$chunk_dir/$seqname.fasta";
	write_fasta($tmp_seq_fn, [{seq=>$seq, seqname=>$seqname}]);
	system("genetack  -t $fsmark_dir/  -m $hmmdef_fn  -f $tmp_seq_fn  >  /dev/null");

	return abs_path("$fsmark_dir/$seqname");
}

###
# Returns reference to array:
#     $todo
#     $todo[i]
#        |
#        |--->{strand}
#        |--->{num_genes}
#        |--->{left}
#        |--->{right}
#        |--->{origin}
#        |--->{_seq}
#        |--->{_seqname}
# 
sub get_fsmark_todo
{
	my($lst, $seq, $high_gc) = @_;

	my @all_genes   = @{$lst->get_all_genes_zero_based()};
	modify_gene_strand(\@all_genes) if $high_gc;
	my $all_contigs = get_contigs_info(\@all_genes);
	
	my($i, $todo) = (0, []);
	foreach my $contig ( @$all_contigs )
	{
		die "Something is wrong: $i > ".@all_genes if $i > @all_genes;

		my @genes = @all_genes[ $i .. $i+$contig->{num_genes}-1];
		my($left_gene, $right_gene) = ($genes[0], $genes[-1]);
		die "Something is wrong: ".Data::Dumper->Dump([$left_gene, $right_gene, $contig],['$left_gene', '$right_gene', '$contig']) if $left_gene->{strand} ne $right_gene->{strand};
		$i += $contig->{num_genes};
		
		my($s_left, $s_right, $seqname, $origin);
		if($contig->{type} eq '+' || $contig->{type} eq '-')      # Usual genes -- cut them out with shoulders
		{
			my $tmp_g = {_group=>'both_shoulders', left=>$left_gene->{left}, right=>$right_gene->{right}};
			($s_left, $s_right) = calculate_shoulder_coords(\@all_genes, $tmp_g, 500, length($seq));
			$seqname = "$contig->{num_genes}-genes chunk ($right_gene->{strand}): $s_left -- $left_gene->{left} .. $right_gene->{right} -- $s_right";
			$origin  = 'normal';
		}
		else
		{
			die "Unknown contig type = '$contig->{type}'";
		}

		my $seq_chunk = get_seq_chunk_with_marked_genes($seq, $s_left, $s_right, \@genes);
		warn "Sequence fragment $s_left-$s_right contains non-ACGT character(s) and will be omitted!\n" and next if $seq_chunk =~ /[^ACGT]/i;

		push @$todo, {
			strand   => $right_gene->{strand},
			left     => $s_left,
			right    => $s_right,
			origin   => $origin,
			num_genes=> $contig->{num_genes},
			_seq     => $seq_chunk,
			_seqname => $seqname,
		};
		
	}

	return $todo;
}

###
# Fighting with FN in high GC genomes: If a gene has coding potential on any frame of opposite strand > 0.6 ,
# then we change strand of the gene.
# 
sub modify_gene_strand
{
	my($all_genes) = @_;

	foreach my $g ( @$all_genes )
	{
		if( $g->{strand} eq '-' )
		{
			if( $g->{cod_p}{f1} > 0.6 || $g->{cod_p}{f2} > 0.6 || $g->{cod_p}{f3} > 0.6 )
			{
				$g->{_modified_strand} = 1;
				$g->{strand} = '+';
			}
		}
		else
		{
			if( $g->{cod_p}{f4} > 0.6 || $g->{cod_p}{f5} > 0.6 || $g->{cod_p}{f6} > 0.6 )
			{
				$g->{_modified_strand} = 1;
				$g->{strand} = '-';
			}
		}
	}
}

###
# Returns:
#    $contigs       --  reference to array of hashes
#     contigs[i]
#        |
#        |--->{num_genes}  
#        |--->{type}             --  '+', '-'
# 
sub get_contigs_info
{
	my($all_genes) = @_;

	my $contig_str = join '', map { $_->{strand} } @$all_genes;

	# Insert delimiters 
	$contig_str =~ s/(.)(?!\1)/$1|/g;   # Separate different gene types: '++-+xxxx+-----xx'  =>  '++|-|+|xxxx|+|-----|xx|'

	my $all_contigs = [];
	foreach my $contig ( split /\|/, $contig_str )
	{
		my $type = substr($contig, 0, 1);

		push @$all_contigs, {
			num_genes => length($contig),
			type      => $type,
			shoulders => 'both_shoulders',
		};
	}

	return $all_contigs;
}

sub save_filter_regs
{
	my($out_fn, $all_filters) = @_;
	
	open(OUT,'>',$out_fn) or die "Can't open file '$out_fn' to write: $!";
	foreach my $filter ( @$all_filters )
	{
		print OUT "$filter->{name}\t$_->{left}\t$_->{right}\n" foreach @{$filter->{ignore_regs}};
	}
	close OUT;
}

sub save_seqs_fasta
{
	my($res, $info_key, $out_fn, %opts) = @_;
	
	print STDERR "Saving '$info_key' sequences to $out_fn...\n";
	
	my @all_seqs = ();
	foreach my $contig ( @$res )
	{
		my $fsm = GeneTack::FSMark->new( $contig->{fsmark_fn} );
		foreach my $fs ( @{$contig->{fs}} )
		{
			next if $fs->{_filter} && !$opts{output_filtered_fs};
			my $info = $fsm->get_seq_info_for_fs_at_pos($fs->{'pos'});
			if($info && $info->{$info_key})
			{
				push(@all_seqs, {seq => $info->{$info_key}, seqname => $fs->{_coord}});
			}
		}
	}
	write_fasta($out_fn, \@all_seqs);
}

sub save_seqs
{
	my($all_items, $out_fn) = @_;
	my @all_seqs = map +{ seq => $_->{_seq}, seqname => $_->{_seqname} }, @$all_items;
	write_fasta($out_fn, \@all_seqs, width => 100);
}

sub save_seqs_info
{
	my($all_items, $out_fn) = @_;

	open(OUT, '>'.$out_fn) or confess("Can't open file '$out_fn' to write: $!");

	print OUT "Num\tStrand\tNum_Genes\tLeft\tRight\tLength\tOrigin\n";

	my $i = 1;
	foreach my $item ( @$all_items )
	{
		print OUT join("\t",
			$i++,
			$item->{strand},
			$item->{origin},
			$item->{num_genes},
			$item->{left},
			$item->{right},
			$item->{right}-$item->{left}+1
		)."\n";
	}

	close OUT;
}

sub get_seq_chunk_with_marked_genes
{
	my($genome, $left, $right, $genes) = @_;

	my $seq_len = $right - $left;
	my $seq     = uc substr($genome, $left, $seq_len);

	# Lower case shoulders
	my $left_s_len = $genes->[0]{left} - $left;
	$seq =~ s/^(.{$left_s_len})/lc($1)/e;
	my $right_s_len = $right - $genes->[-1]{right};
	$seq =~ s/(.{$right_s_len})$/lc($1)/e;

	# Lower case intergenic regions 
	for(my $i=0; $i < scalar(@$genes)-1; $i++)
	{
		my($offset, $len);
		if( $genes->[$i+1]{left} <= $genes->[$i]{right} )   # Gene overlap
		{
			$offset = $genes->[$i+1]{left} - $left;
			$len    = $genes->[$i]{right} - $genes->[$i+1]{left} + 1;
		}
		else                                                # Adjacent genes
		{
			$offset = $genes->[$i]{right} - $left + 1;
			$len    = $genes->[$i+1]{left} - $genes->[$i]{right} - 1;
		}
		$seq = substr($seq, 0, $offset).lc(substr($seq, $offset, $len)).substr($seq, $offset+$len);
		die "Something wrong for $offset/$len different seq length: ".length($seq)." != $seq_len" if length($seq) != $seq_len;
	}

	return $genes->[0]{strand} eq '-' ? revcomp($seq) : $seq;
}

###
# Arguments:
#     $g
#      |
#      |--->{_group}  --  'both_shoulders' or 'down_shoulder'
#      |--->{right}   --  current right coordinate to add shoulder to
#      |--->{left}    --  current left coordinate to add shoulder to
#      |--->{strand}
# 
sub calculate_shoulder_coords
{
	my($all_genes, $g, $shoulder_len, $genome_len) = @_;

	# $o_gene -- overlap gene
	my($s_left, $s_right);
	if($g->{_group} eq 'both_shoulders')                    # Shoulders from both sides
	{
		$s_right = ($g->{right}+$shoulder_len) > $genome_len ? $genome_len : ($g->{right}+$shoulder_len);
		$s_right = get_g_left_coord_in_region($all_genes, $g->{right}, $s_right); # Adjust shoulder coordinate wrt other genes in the shoulder
		
		$s_left  = ($g->{left}-$shoulder_len) <= 0 ? 0 : ($g->{left}-$shoulder_len);
		$s_left  = get_g_right_coord_in_region($all_genes, $s_left, $g->{left});  # Adjust shoulder coordinate wrt other genes in the shoulder
	}
	elsif($g->{_group} eq 'down_shoulder')                  # Downstream shoulder only
	{
		if( $g->{strand} eq '+' )
		{
			$s_left  = $g->{left};
			
			$s_right = ($g->{right}+$shoulder_len) > $genome_len ? $genome_len : ($g->{right}+$shoulder_len);
			$s_right = get_g_left_coord_in_region($all_genes, $g->{right}, $s_right); # Adjust shoulder coordinate wrt other genes in the shoulder
		}
		else
		{
			$s_right = $g->{right};
			
			$s_left  = ($g->{left}-$shoulder_len) <= 0 ? 0 : ($g->{left}-$shoulder_len);
			$s_left  = get_g_right_coord_in_region($all_genes, $s_left, $g->{left});  # Adjust shoulder coordinate wrt other genes in the shoulder
		}
	}
	elsif( $g->{_group} eq 'no_shoulders' )                  # In order not to make if in the where we call this function
	{
		$s_left  = $g->{left};
		$s_right = $g->{right};
	}
	else
	{
		die "calculate_shoulder_coords: unknown group: '$g->{_group}'";
	}

	# Final adjustment of coordinates
	$s_left  = $g->{left}  if $g->{left}  < $s_left;
	$s_right = $g->{right} if $g->{right} > $s_right;

	return($s_left, $s_right);
}

###
# If no gene in the region found $reg_left is returned
sub get_g_right_coord_in_region
{
	my($all_genes, $reg_left, $reg_right) = @_;

	my $g_right_coord = undef;
	foreach my $g ( @$all_genes )
	{
		# If either left or right end of gene is inside the region -- this is what we need!
		if( $reg_left < $g->{right} && $g->{right} < $reg_right )
		{
			$g_right_coord = $g->{right} if !$g_right_coord || $g->{right} > $g_right_coord;
		}
		elsif( $g->{left} < $reg_right && $reg_right < $g->{right} )
		{
			$g_right_coord = $reg_right;  # Gene overlap -- no shoulder
			last;
		}
	}

	return $g_right_coord || $reg_left;
}

###
# If no gene in the region found $reg_right is returned
sub get_g_left_coord_in_region
{
	my($all_genes, $reg_left, $reg_right) = @_;

	my $g_left_coord = undef;
	foreach my $g ( @$all_genes )
	{
		# If either left or right end of gene is inside the region -- this is what we need!
		if( $reg_left < $g->{left}  && $g->{left}  < $reg_right )
		{
			$g_left_coord = $g->{left} if !$g_left_coord || $g->{left} < $g_left_coord;
		}
		elsif( $g->{left} < $reg_left && $reg_left < $g->{right} )
		{
			$g_left_coord = $reg_left;  # Gene overlap -- no shoulder
			last;
		}
	}

	return $g_left_coord || $reg_right;
}

sub revcomp
{
	my($seq) = @_;
	$seq = reverse $seq;
	$seq =~ tr/ACGTacgt/TGCAtgca/;
	return $seq;
}

sub usage
{
	my($msg) = @_;
	$msg = $msg ? $msg."\n" : '';

	my $script = File::Spec->splitpath($0);
	return"$msg
GeneTack-GM, version $VERSION

USAGE:
    $script   [OPTIONS]   <sequence.fasta>   >   <out.fsmarkgm>

OPTIONS:
    --num_threads         <N>           --  number of threads to use. Default is $NUM_THREADS
    --output_filtered_fs                --  output ALL predicted FS + add Filter column in the output table
    --save_seqs           <FILE.fasta>  --  save sequences of all regions analyzed by FSMark
    --save_filter_regs    <FILE.txt>    --  save ingnored regions where frameshifts are filtered out
    --save_seqs_info      <FILE.txt>    --  save information about sequences analyzed by FSMark
    --save_fsgene_seqs    <FILE.fna>    --  save nt sequences of all frameshifted CDSs 
    --save_fsprot_seqs    <FILE.faa>    --  save sequences of all proteins conceptually translated from predicted FS-genes
    --save_gm_mod         <FILE.mod>    --  save model that was used for gene prediction
    --save_fs_mod         <FILE.mod>    --  save model that was used to generate the .hmm_def file
    --save_fsmark_files   <DIR>         --  save FSMark output files in the directory. Name of a file is 
                                            the left coordinate of contig.

    --fs_prob             <NUM>         --  probability of a frameshift -- probability of direct transition from one 
                                            coding state to another. Default value is $FS_PROB
    --gm_order            <NUM>         --  order to run GeneMark and predict genes in genome. Default value is $GM_ORDER
    --gm_cdec             <NUM>         --  value of CDEC option for GeneMark HMM. Default value is $GM_CDEC
    --cod_order           <NUM>         --  order of emission strings in coding states in FSMark HMM.
                                            Supported values: ".join(',',sort keys %{$GeneTack::FSMark::ORDER_COD_FREQS}).". Default value is $COD_ORDER
    --high_gc_thr         <NUM>         --  threshold to consider genome as high GC. Default value $HIGH_GC_THR
    --calc_fs_scores                    --  calculate score for every predicted FS (turned off by default).
    --gm_mod_fn           <FILE.mod>    --  use this model for gene prediction
    --fs_mod_fn           <FILE.mod>    --  use this model to generate .hmm_def file and predict frameshifts
    --tmp_dir             <DIR>         --  where the temporary folder will be created. Default value: '$TMP_DIR'

FILTER OPTIONS:
    --ignore_vicinity     <NUM>         --  vicinity around special places (such as short gene overlap region)
                                            where we ignore predicted FS. Default value is $IGNORE_VICINITY
    --fs_cod_shoulder_len <NUM>         --  min length of coding region shoulders which flank predicted frame shift.
                                            All other frame shifts will be filtered out. Default value is $FS_COD_SHOULDER_LEN
    --min_overlap_len     <NUM>         --  frame shifts caused by gene overlaps with overlap region smaller than this
                                            value will be ignored. Default value is $MIN_OVERLAP_LEN
    --min_adj_len         <NUM>         --  frame shifts caused by adjacent genes with gap length less than this
                                            value will be ignored. Default value is $MIN_ADJ_LEN
    --max_overlap_rbs     <NUM>         --  maximum value of RBS site score of downstream gene in gene overlap.
                                            Default value is $MAX_OVERLAP_RBS
    --max_adjacent_gap    <NUM>         --  maximum length of the gap between two genes to say that they are adjacent.
                                            Default value is $MAX_ADJACENT_GAP
    --max_adjacent_rbs    <NUM>         --  maximum value of RBS site score of downstream gene in adjacent genes.
                                            Default value is $MAX_ADJACENT_RBS
    --min_lg_fs_score     <NUM>         --  frameshifts with the LF_FS_SCORE less than this value will be ignored. Filter is
                                            enabled if --calc_fs_scores option is used only. Default value is $MIN_LG_FS_SCORE
    --adj_stop_extra_len  <NUM>         --  extra length (+ gap length) to the upstream region of downstream gene in
                                            adjacent genes to search for STOP codon in. Default value is $ADJ_STOP_EXTRA_LEN

";
}

