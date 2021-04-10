package GeneTack::GeneMark;

use strict;
use warnings;

# $Id$

###
# Ivan Antonov (antonov1986@gmail.com)
#
###
#
# $self
#   |
#   |--->{init_fn}  -- file name that was read and its data were used to initialize the object
#

use Data::Dumper;
use Carp qw(confess); 
use Cwd  qw(abs_path);

###
# CONSTANTS
my $GM_BIN = $ENV{HOME}.'/bin/';

###
# CONSTRUCTOR
sub new
{
	my $class = shift;

	my $self = bless {}, $class;
	
	return $self;
}

###
# PRIVATE METHODS

###
# Description:
#     This function will be called inside new() method of all child classes
# 
# Arguments:
#     %opts
#        |
#        |--->{init_fn}  -- to set in $self
# 
sub _init
{
	my $self = shift;
	my(%opts) = @_;
	$self->{init_fn} = $opts{init_fn};
}

###
# PUBLIC METHODS

sub init_fn
{
	my($self) = @_;
	return $self->{init_fn};
}

###
# SUBROUTINES

###
# Run GeneMark and returns names of created files:
#   * $fs_mod_fn
#   * $lst_fn
#   * $gdata_fn   --  [optional] file created only if $opts{gdata} flag is received.
#                     This value is undef in absence of $opts{gdata}.
#   * gm_mod      --  model to run GeneMark to predict genes
#   * fs_mod      --  model to be returned that will be used to build .hmm_def
#
sub run_gm
{
	my($seq_fn, $fs_order, $gm_order, $gm_cdec, %opts) = @_;

	my $orig_cwd = undef;
	if( $opts{run_dir} )    # Becase GeneMark creates many files in CWD
	{
		$orig_cwd = abs_path('.');
		chdir $opts{run_dir};
		$opts{run_dir} .= '/' if $opts{run_dir} !~ /[\\\/]$/;
	}
	else
	{
		$opts{run_dir} = '';
	}

	$opts{genemark_path} = $GM_BIN if !defined $opts{genemark_path};
	# add '/' to the path, if it is provided
	$opts{genemark_path} .= '/' if $opts{genemark_path} =~ /\S/ && $opts{genemark_path} !~ /[\\\/]$/;
	
	my($base_name, $gm_mod_fn, $gm_mat_fn) = ('MODEL', undef, undef);
	if( $opts{gm_mod_fn} )
	{
		$gm_mod_fn = $opts{gm_mod_fn};
	}
	else
	{
		$gm_mod_fn = $opts{run_dir}.$base_name.'.gm_hmm.mod'; # Generate model file name for GeneMark.hmm
		$gm_mat_fn = $opts{run_dir}.$base_name.'.gm_gm.mat';  # Generate mattrix file name for GeneMark

		my $run = "$opts{genemark_path}gmsn.pl --order $gm_order --gm --clean --name ${base_name}.gm $seq_fn";
		warn "Creating model file: $run\n";
		system("$run > /dev/null");
=comment
# After in run_dir will be create:
	GENOME.fasta.lst
	gms.log
	MODEL.gm_gm.mat
	MODEL.gm_hmm_heuristic.mod
		(sometimes MODEL.gm_gm_heuristic.mat)
	MODEL.gm_hmm.mod
=cut
		&modify_mod_file( $gm_mod_fn, [{re=>qr/^CDEC\s+\d+/, subst=>"CDEC $gm_cdec"}] ) || return undef;
	}

	my $lst_fn = $opts{run_dir}.$base_name.'.lst';
	warn "Locating genes: $opts{genemark_path}gmhmmp -r -k -m $gm_mod_fn -o $lst_fn $seq_fn \n";
	system("$opts{genemark_path}gmhmmp -r -k -m $gm_mod_fn -o $lst_fn $seq_fn");

	my $gdata_fn = undef;
	if( $opts{gdata} && $gm_mat_fn )
	{
		# /home/ivan/_my/tmp/1/1Mb.fasta   =>   1Mb.fasta
		my($name_wo_path) = $seq_fn =~ /([^\\\/]+)$/;
		$gdata_fn = $opts{run_dir}.$name_wo_path.'.gdata';
		warn "Creating gdata file: $opts{genemark_path}gm -D -m $gm_mat_fn $seq_fn \n";
		system("$opts{genemark_path}gm -D -m $gm_mat_fn $seq_fn");
	}

	my $fs_mod_fn = undef;
	if( $opts{fs_mod_fn} )
	{
		$fs_mod_fn = $opts{fs_mod_fn};
	}
	else
	{
		$fs_mod_fn = $opts{run_dir}.$base_name.'.fs_hmm.mod';           # Generate model file name for FSMark
		my $run = "$opts{genemark_path}gmsn.pl --order $fs_order --clean --name ${base_name}.fs $seq_fn";
		warn "Creating model file: $run\n";
		system("$run > /dev/null");
=comment
# After in run_dir will be create:
	gms.log
	GENOME.fasta.lst
	GENOME.fasta.gdata
	GENOME.fasta.ldata
	GENOME.fasta.lst
	MODEL.fs_hmm_heuristic.mod
	MODEL.fs_hmm.mod
	MODEL.lst
=cut

	}

	chdir($orig_cwd) if $orig_cwd;

	return( $gm_mod_fn, $fs_mod_fn, $lst_fn, $gdata_fn );
}

###
# $todo
# $todo[i]
#    |
#    |--->{re}     --  RE to confirm that we need to modify the string
#    |--->{subst}  --  substitution string
# 
sub modify_mod_file
{
	my($mod_fn, $todo) = @_;
	my $content = '';

	if( open(my $in_fh, "<", $mod_fn) ){
		while( my $str = <$in_fh> )
		{
			for my $item ( @$todo )
			{
				if( $str =~ $item->{re} )
				{
					$item->{_applied}++ if $str =~ s/$item->{re}/$item->{subst}/;
				}
			}
			$content .= $str;
		}
		close $in_fh;
	}else{
		warn "Can't open file '$mod_fn': $!";
		return undef;
	}

	warn "modify_mod_file: some TODOs didn't work!" if grep { !$_->{_applied} } @$todo;

	open( my $out_fh, ">", $mod_fn) || die "Can't open file '$mod_fn' to write: $!";
	print $out_fh $content;
	close $out_fh;

	return 1;
}

1;
