package GeneTack::GeneMark::Meta;

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
# 

use Data::Dumper;
use Carp qw(confess); 

###
# CONSTANTS

###
# CONSTRUCTOR
# 
sub new
{
	my $class = shift;
	my($fn) = @_;
	
	my $self = bless {
		genes => [],
	}, $class;
	
	$self->_read_fn($fn);
	
	return $self;
}

###
# PRIVATE METHODS

sub _read_fn
{
	my $self = shift;
	my($fn) = @_;
	
	open(my $fh, '<', $fn) or die "Can't open file '$fn': $!";
	
	my $seqname = undef;
	while( <$fh> )
	{
		if( /^FASTA definition line\: (.+)/ )
		{
			$seqname = $1;
		}
		elsif( /^\s+\d/ )
		{
			s/^\s*//;
			my @vals = split /\s+/;
			push @{$self->{genes}}, {
				seqname => $seqname,
				id      => $vals[0],
				strand  => $vals[1],
				left    => $vals[2],
				right   => $vals[3],
				len     => $vals[4],
				class   => $vals[5],
				spacer  => $vals[6],
				rbs     => $vals[7],
			};
		}
	}
	
	close $fh;
}

###
# PUBLIC METHODS

sub num_genes
{
	my $self = shift;
	return scalar( @{$self->{genes}} );
}

sub get_all_genes
{
	my $self = shift;
	return $self->{genes};
}

sub get_gene
{
	my $self = shift;
	my($i) = @_;
	return $self->{genes}[$i];
}

1;

