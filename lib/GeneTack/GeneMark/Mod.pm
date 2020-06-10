package GeneTack::GeneMark::Mod;

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
#   |--->{_alphabet}     --  [aref]
#   |
#   |--->{order}         --  [number] order of a model
#   |
#   |--->{freqs}
#           |
#           |--->{COD1}  --  [reference to hash of hashes: {frame}->{codon}] typical genes
#           |       |
#           |       |--->{1}
#           |       |     |
#           |       |     |--->{AAA}
#           |       |     |--->{AAC}
#           |       |    ...
#           |       |
#           |       |--->{2}
#           |       |     |
#           |       |    ...
#           |       |
#           |       |--->{3}
#           |             |
#           |            ...
#           |
#           |--->{COD2}  --  (optional) atypical genes -- if we read _hmm_combine.mod file only
#           |       |
#           |      ...
#           |
#           |--->{NONC}  --  non-coding region
#                   |
#                   |--->{AAA}
#                   |--->{AAC}
#                  ...
#

use Data::Dumper;
use Carp qw(confess); 

use base qw(GeneTack::GeneMark);

###
# CONSTANTS

###
# CONSTRUCTOR
sub new
{
	my $class = shift;
	my($fn) = @_;
	
	my $self = bless {
		_alphabet => ['A', 'C', 'G', 'T'],
		freqs     => {
			COD1 => {1=>{}, 2=>{}, 3=>{}},
			COD2 => {1=>{}, 2=>{}, 3=>{}},
			NONC => {},
		},
	}, $class;
	
	$self->_init($fn);
	
	$self->SUPER::_init(init_fn => $fn);
	
	return $self;
}

###
# PRIVATE METHODS

sub _init
{
	my $self = shift;
	my($fn) = @_;
	
	open(my $fh, '<', $fn) or die "Can't open file '$fn': $!";

	# Load order
	while(<$fh>){ /^ORDM\s+(\d+)/ and $self->{order} = $1 and last };
	confess("Someting wrong: can't find ORDM in the file '$fn'") if $_ !~ /^ORDM/;
	
	# Skip lines until COD1
	while(<$fh>){ /^COD1/ and last };
	confess("Someting wrong: can't find COD1 here") if $_ !~ /^COD1/;
	
	# Load COD1
	$self->_load_cod_freqs($fh, 'COD1');
	
	$_ = <$fh>;
	$_ = <$fh> while $_ =~ /^\s*$/;
	if( /^COD2/ ) {
		$self->_load_cod_freqs($fh, 'COD2');
	}
	elsif( /^NONC/ ) {
		$self->_load_nonc_freqs($fh);
	}
	else {
		confess("Someting wrong: I see '$_' but expect 'COD2' or 'NONC' here");
	}

	$_ = <$fh>;
	$_ = <$fh> while $_ =~ /^\s*$/;
	if( /^COD2/ ) {
		$self->_load_cod_freqs($fh, 'COD2');
	}
	elsif( /^NONC/ ) {
		$self->_load_nonc_freqs($fh);
	}
	else {
		confess("Someting wrong: I see '$_' but expect 'COD2' or 'NONC' here");
	}
	
	close $fh;
}

sub _load_nonc_freqs
{
	my $self = shift;
	my($fh)  = @_;
	
	foreach my $cod ( @{$self->_get_all_codons} )
	{
		# Here $line will be = '0.03313'
		(my $line  = <$fh>) =~ s/[\n\r]//g;
		confess("Someting wrong: line = '$line'") if $line !~ /^\d+\.\d+$/;
		confess("Codon duplication found '$cod'") if exists $self->{freqs}{NONC}{$cod};
		$self->{freqs}{NONC}{$cod} = $line;
	}
}

sub _load_cod_freqs
{
	my $self = shift;
	my($fh, $name)  = @_;

	foreach my $cod ( @{$self->_get_all_codons} )
	{
		# Here $line will be = '0.03313 0.02450 0.01561'
		my $line  = <$fh>;
		my @freqs = split /\s+/, $line;
		confess("Someting wrong for $name (codon '$cod'): ".Data::Dumper->Dump([\@freqs], ['freqs'])) if @freqs != 3;
		confess("Codon duplication found '$cod'") if exists $self->{freqs}{$name}{1}{$cod};
		$self->{freqs}{$name}{$_+1}{$cod} = $freqs[$_] foreach (0..2);
	}
}

sub _get_all_codons
{
	my $self = shift;

	my $alpha_size  = scalar(@{$self->{_alphabet}});
	my $cur_indexes = [ map { 0 } 0..$self->{order} ];
	my @all_codons  = ( join('', map { $self->{_alphabet}[ $cur_indexes->[$_] ] } 0..$self->{order}) );
	while( $cur_indexes = get_next_index_set($cur_indexes, $alpha_size-1) )
	{
		push @all_codons, join('', map { $self->{_alphabet}[ $cur_indexes->[$_] ] } 0..$self->{order});
	}
	
	return \@all_codons;
}

sub get_next_index_set
{
	my($cur_indexes, $last_index) = @_;
	
	my @next_indexes = @$cur_indexes;
	for(my $i = $#$cur_indexes; $i >= 0; --$i)
	{
		if( $cur_indexes->[$i] == $last_index )
		{
			$next_indexes[$i] = 0;
		}
		else
		{
			$next_indexes[$i]++;
			last;
		}
	}
	
	# $cur_indexes consists of zeros only at the very first call of function or when we completed whole cicle.
	my $to_return = join('', @next_indexes) =~ /^0*$/ ? undef() : \@next_indexes;
	return $to_return;
}

###
# PUBLIC METHODS

###
# Description:
#     Typical genes are assumed
# 
# Arguments:
#     $frame  --  [number] 1, 2 or 3
# 
sub get_coding_freqs
{
	my $self = shift;
	my($frame) = @_;
	confess("Wrong frame = '$frame'") if !exists $self->{freqs}{COD1}{$frame};
	return $self->{freqs}{COD1}{$frame};
}

sub get_noncoding_freqs
{
	my $self = shift;
	return $self->{freqs}{NONC};
}

sub get_order
{
	my $self = shift;
	return $self->{order};
}

###
# SUBROUTINES

1;

