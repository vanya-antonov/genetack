package GeneTack::GeneMark::Lst;

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
#   |--->{genes}    --  reference to array of hashes
#   |     genes[i]
#   |       |
#   |       |--->{id}
#   |       |--->{strand}
#   |       |--->{left}
#   |       |--->{right}
#   |       |--->{len}
#   |       |--->{class}
#   |       |--->{frame}      --  [number] 1-6
#   |       |--->{spacer}     --  optional
#   |       |--->{rbs_score}  --  optional
#   |       |--->{overlap}    --  [boolean]
#   |       |--->{cod_p}      --  [optional] coding potential in 6 different frames (if gdata_fn file is provided)
#   |               |
#   |               |
#   |               |--->{f1}
#   |               |--->{f2}
#   |               |--->{f3}
#   |               |--->{f4}
#   |               |--->{f5}
#   |               |--->{f6}
#   |
#   |--->{overlaps}  --  reference to array of arrays: overlaps[i] = [gene1_i, gene2_i, ...]
#
###
#
# For '+' strand gene:
#  ATG ... TGA
#  ^         ^
#  |         |
#  left    right
# 
# !!! First genome position is 1 !!!
# 

use Data::Dumper;
use Carp qw(confess cluck);
use Storable qw(dclone);

use base qw(GeneTack::GeneMark);

###
# CONSTANTS
my %OPT_COLS = (
	Spacer => {
		key => 'spacer',
	},
	RBS_score => {
		key => 'rbs_score',
	},
	RBS => {
		key => 'rbs_score',
	},
);

###
# CONSTRUCTOR
# 
# %opts
#   |
#   |--->{gdata_fn}  --  .gdata file with information about coding potential in 6 frames
# 
sub new
{
	my $class = shift;
	my($fn, $refID, %opts) = @_;

	my $self = bless {
		genes    => [],
		overlaps => [],
	}, $class;

	$self->_init( $fn, $refID ) || return undef;

	# add $self->{genes}->[i]->{cod_p} (coding potential) to each gene
	$self->_calculate_coding_potential($opts{gdata_fn}) if $opts{gdata_fn};

	$self->SUPER::_init( init_fn => $fn ); # GeneTack::GeneMark::_init

	return $self;
}

###
# PRIVATE METHODS
 
sub _calculate_coding_potential
{
	my $self = shift;
	my($fn) = @_;

	my $gdata = _read_gdata_fn($fn);

	my $last_left_i = 0; # for optimization purposes
	foreach my $gene ( @{$self->{genes}} )
	{
		my @cod_vals = ();
		
		foreach my $val ( @$gdata[$last_left_i..$#$gdata] )
		{
			push(@cod_vals, $val) if $gene->{left} < $val->{seq_i} && $val->{seq_i} < $gene->{right};
			if( $val->{seq_i} > $gene->{right} )   # No reason to look further
			{
				$last_left_i = $val->{i} - 10; # -10 back in case of gene overlap
				last;
			}
		}

		foreach my $frame ( qw(f1 f2 f3 f4 f5 f6) )
		{
			$gene->{cod_p}{$frame} = @cod_vals == 0 ? 0 : aver( map { $_->{$frame} } @cod_vals);
		}
	}
}

sub _read_gdata_fn
{
	my($fn) = @_;
	
	open(my $fh, '<', $fn) or die "Can't open file '$fn': $!";
	my($gdata, $i) = ([], 0);
	while( <$fh> )
	{
		next if /^#/ || /^\s*$/;
		# here $_ = 828	0.000755	0	0	0	0	0
		s/[\n\r]//g;
		my @vals = split /\s+/;
		die "Something wrong with the string '$_': $!" if @vals != 7;
		push @$gdata, {
			i     => $i++,
			seq_i => $vals[0],
			f1    => $vals[1],
			f2    => $vals[2],
			f3    => $vals[3],
			f4    => $vals[4],
			f5    => $vals[5],
			f6    => $vals[6],
		};
	}
	close $fh;
	
	return $gdata;
}

sub _init
{
	my $self = shift;
	my( $fn, $refID ) = @_;

	$self->_read_fn( $fn, $refID ) || return undef;

	$self->_calculate_frames();
	$self->_find_overlaps_in_strand('+');
	$self->_find_overlaps_in_strand('-');

	return 1;
}

###
# Text below taken from the .ldata file generated by gm:
# 
#  Frames are counted 1-6, Direct 1-3 and Complement 4-6
#  The first base of the input sequence is numbered 1 not 0
#  The frame is computed by the direct strand gene start
#    (LEnd for direct strand, REnd for complement strand)
#     Frame=(LEnd-1)%3+1  or (REnd-1)%3+4
# 
sub _calculate_frames
{
	my $self = shift;
	
	foreach my $gene ( @{$self->{genes}} )
	{
		if( $gene->{strand} eq '+' )
		{
			$gene->{frame} = ($gene->{left}-1)%3 + 1;
		}
		else
		{
			$gene->{frame} = ($gene->{right}-1)%3 + 4;
		}
	}
}

sub _find_overlaps_in_strand
{
	my $self    = shift;
	my($strand) = @_;
	
	my @all_strand_i = grep { $self->{genes}[$_]{strand} eq $strand } 0..@{$self->{genes}}-1;
	return if @all_strand_i < 2;
	
	my $this_i = shift @all_strand_i;
	my $next_i = shift @all_strand_i;
	
	while ( @all_strand_i )
	{
		my $item       = [$this_i];                       # list of consecuite overlapping genes
		my $this_right = $self->{genes}[$this_i]{right};  # right coordinate of the gene $i
		my $next_left  = $self->{genes}[$next_i]{left};   # left coordinate of next gene: $i+1
		
		while( $next_left <= $this_right )
		{
			# If we come here overlap is found!
			push @$item, $next_i;
			$self->{genes}[$this_i]{overlap} = 1;
			$self->{genes}[$next_i]{overlap} = 1;
			
			# Do not process last gene of the strand
			last if @all_strand_i == 0;
			
			# Step forward
			$this_i = $next_i;
			$next_i = shift @all_strand_i;
			
			$this_right = $self->{genes}[$this_i]{right};
			$next_left  = $self->{genes}[$next_i]{left};
		}
		push @{$self->{overlaps}}, $item if @$item > 1;
		
		# Do not process last gene of the strand
		last if @all_strand_i == 0;
		
		# Step forward
		$this_i = $next_i;
		$next_i = shift @all_strand_i;
	}
}


sub _read_fn
{
	my $self = shift;
	my( $fn, $refID ) = @_;	# MODEL.lst file, and sequence ID

	my @opt_keys;
	my $i;
	my $flg = 0;
	my $emsg = "Sequence '$refID' in file '$fn'";

	open(my $fh, '<', $fn) or die "Can't open file '$fn': $!";
	while(<$fh>){
=comment
...
 FASTA definition line: NZ_KI535340.1 Abiotrophia defectiva ATCC 49176 Scfld0, whole genome shotgun sequence
Predicted genes
   Gene    Strand    LeftEnd    RightEnd       Gene     Class   Spacer     RBS
    #                                         Length            Length    score
    1        -          46         834          789        1        7   3.3786
=cut
		if( /FASTA\s+definition\s+line:\s+(\S+)/i ){

			if( $refID eq $1 ){
				$flg = 1; # Found a result with the required ID
			}elsif( $flg > 1 ){
				last;	# Finish the search
			}else{
				$flg = 0; # Skip lines
			}
			next;
		}

		if( $flg == 1 and /^Predicted\s+genes/i ){ # and last;

#   Gene    Strand    LeftEnd    RightEnd       Gene     Class   Spacer     RBS
			my $head = <$fh>;

			@opt_keys = ();

			for my $name ( $head =~ /\s+Class\s+(\S+)\s+(\S+)\s*$/i ) # From optional column names to keys: Spacer     RBS
			{
				confess("File '$fn' has unknown optional column '$name'") if !exists $OPT_COLS{ $name };
				push @opt_keys, $OPT_COLS{ $name }{'key'};
			}
			confess("$emsg does NOT have optional columns") unless @opt_keys;

			++$flg;
			next;
		}
		next if $flg < 2;

		next if /^\s+$/ || /^\s*#/; # Skip line: #                                         Length            Length    score
		s/^\s+//;
=comment
1439        +          <3         151          150        1       -1   0.0333
1591        -      178725     >179279          555        1       13  -1.6296
=cut
		my @vals = split /\s+/;
		confess("$emsg has WRONG FORMAT: ".Data::Dumper->Dump([\@vals], ['vals'])) if @vals != 6+@opt_keys;

		# Substitute coordinates: '<3' --> '3' and >179279 --> 179279
		$vals[$_] =~ s/\D+//g for 2, 3;

		my %opt_vals = map { $opt_keys[$_] => $vals[6+$_] } 0 .. $#opt_keys;
		push @{ $self->{'genes'} }, {
			id      => ++$i,
			gmID    => $vals[0], # geneMark ID
			strand  => $vals[1],
			left    => $vals[2],
			right   => $vals[3],
			len     => $vals[4],
			class   => $vals[5],
			overlap => 0,
#			refID   => $refID,
			%opt_vals,
		};

	}
	close $fh;

	cluck("$emsg does NOT have \x1b[31mpredicted genes\x1b[0m") unless $i; # was confess()
	return $i;
}


###
# Description:
#     Returns hash with information about the overlap.
#     Warning: SPECIFIED MUST BE INDECES OF OVERLAPPING GENES. METHOD DOESN'T DOUBLE CHECK.
# 
###
# Arguments:
#     $genes_i  -  reference to array with indeces of CONSECUTIVE OVERLAPPING genes. May contains > 2 indeces.
# 
###
# Returns:
#     $info   --  reference to hash
#        |
#        |--->{genes_i}    --  reference to array of gene indeces
#        |--->{num_genes}  --  [number] number of genes in the overlap
#        |--->{len}        --  [number] sum of lengths of all overlap regions
#        |--->{total_len}  --  [number] total length including non-overlapping regions -- from first gene start to last gene end
#        |--->{left}       --  [number] left end of the first gene in the soverlap
#        |--->{right}      --  [number] rigth end of the last gene in the overlap
#        |--->{strand}     --  '+' or '-'
#        |
#    ======================
#    For pair overlaps only
#    ======================
#        |
#        |--->{len_before_overlap}  --  legth of the upstream gene chunk before overlap
#        |--->{ovlp_type}           --  how many nucleotides overlapping genes are shifted wrt to each other.
#        |                              Namely, how many nucleotides we need to shift any gene in downstream direction so both genes 
#        |                              are in the same frame. So overlaps ATGA and ATGxTGA has type 1 and ATGxxTGA has type 2.
#        |--->{ovlp_seq}            --  overlap sequence (if option seq of provided)
###
#     %opts
#       |
#       |--->{seq}  --  genome sequence to cut overlap sequences from
# 
sub _get_overlap_info
{
	my $self = shift;
	my($genes_i, %opts) = @_;

	my @genes = map { $self->{genes}[$_] } @$genes_i;
	
	my $info = {
		genes_i   => $genes_i,
		num_genes => scalar(@genes),
		len       => sum( map { $genes[$_]{right}-$genes[$_+1]{left}+1 } 0..$#genes-1 ),
		total_len => $genes[-1]{right} - $genes[0]{left} + 1,
		left      => $genes[0]{left},
		right     => $genes[-1]{right},
		strand    => $genes[0]{strand},
	};

	if ( @genes == 2 )
	{
		$info->{len_before_overlap} = $genes[1]{left} - $genes[0]{left};
		$info->{ovlp_type} = $info->{len}%3;
		warn "Overlapping genes in the same frame found: ".Data::Dumper->Dump([\@genes,$info], ['genes','info']) if $info->{ovlp_type} == 0;
		if( $opts{seq} )
		{
			my $ovlp_seq = substr($opts{seq}, $info->{left}+$info->{len_before_overlap}-1, $info->{len});
			$info->{ovlp_seq} = $info->{strand} eq '-' ? revcomp($ovlp_seq) : $ovlp_seq;
			confess("Wrong overlap sequence '$ovlp_seq'".Data::Dumper->Dump([\@genes, $info], ['genes', 'info'])) if $info->{ovlp_seq} !~ /^(ATG|GTG|TTG)/i || $info->{ovlp_seq} !~ /(TAG|TGA|TAA)$/i;
		}
	}
	
	if( $info->{len} <= 0 )
	{
		confess("Something wrong with overlap: ".Data::Dumper->Dump([\@genes, $info], ['genes', 'info']));
	}
	
	return $info;
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
	return dclone($self->{genes});
}

sub get_all_genes_zero_based
{
	my $self = shift;

	my @genes = ();
	for my $g ( @{$self->{genes}} )
	{
		my %h = %$g;
		$h{left}--;
		push(@genes, \%h);
	}

	return \@genes;
}

sub get_gene
{
	my $self = shift;
	my($i) = @_;
	return $self->{genes}[$i];
}

sub num_overlaps
{
	my $self = shift;
	return scalar( @{$self->{overlaps}} );
}

###
# Description:
#     This function considers overlap of 3 or more consecutive genes as a single overlap.
#     To get pair overlaps only, use get_pair_overlaps() method
# 
###
# Returns:
#     $overlaps  --  reference to array of hashes. See _get_overlap_info() function for more details.
# 
sub get_all_overlaps_info
{
	my $self = shift;
	my @info = map { $self->_get_overlap_info($_) } @{$self->{overlaps}};
	return [ sort { $a->{left} <=> $b->{left} } @info ];
}

###
# Description:
#     This function considers pair overlaps only, i.e. it splits overlap of 3 consecutive 
#     genes into two pair overlaps.
#     So, scalar(@pair_ovelaps) >= scalar(@{$self->{overlaps}})
# 
sub get_pair_overlaps_info
{
	my $self = shift;
	my(%opts) = @_;
	my @info = map { $self->_get_overlap_info($_, %opts) } @{$self->get_pair_overlaps};
	return [ sort { $a->{left} <=> $b->{left} } @info ];
}

sub get_pair_overlaps
{
	my $self = shift;
	
	my @pair_ovelaps = ();
	foreach my $overlap ( @{$self->{overlaps}} )
	{
		for my $i ( 0 .. $#$overlap-1 )
		{
			push @pair_ovelaps, [$overlap->[$i], $overlap->[$i+1]];
		}
	}
	return \@pair_ovelaps;
}

sub aver_pair_overlap_len
{
	my $self = shift;
	return aver( map { $_->{len} } @{$self->get_pair_overlaps_info} );
}

sub aver_len_from_gene_start_to_overlap
{
	my $self = shift;
	return aver( map { $_->{len_before_overlap} } @{$self->get_pair_overlaps_info} );
}

###
# $gap_len_limit -- maximum length of the gap between two genes for the pair to be considered as adjacent
sub get_all_adjacent_gene_pairs
{
	my $self = shift;
	my ($gap_len_limit) = @_;
	
	my @res;
	for (my $up_i = 0; $up_i < $self->num_genes(); $up_i++)
	{
		my $up_gene = $self->{genes}[$up_i];
		my($down_gene, $gap_len) = $self->get_closest_downstream_gene( $up_gene );
		next if !defined $down_gene;
		if( $gap_len < $gap_len_limit )
		{
			my $strand = $down_gene->{strand};
			my $left   = $strand eq '+' ? $up_gene->{left}    : $down_gene->{left};
			my $right  = $strand eq '+' ? $down_gene->{right} : $up_gene->{right};
			push @res, {
				up_gene   => $up_gene,
				down_gene => $down_gene,
				gap_len   => $gap_len,
				strand    => $strand,
				left      => $left,
				right     => $right,
			};
		}
	}
	
	return \@res;
}

###
# Arguments:
#    %opts
#       |
#       |--->{overlap}   --  [boolean] 
#       |--->{adjacnet}  --  [integer] gap length limit for adjacent genes
# 
sub get_all_genes_that_are_NOT
{
	my $self = shift;
	my(%opts) = @_;
	
	# Deep copy of $self->{genes}
	my $genes = dclone($self->{genes});
	
	# Filter out adjacent genes if needed
	if( $opts{adjacent} )
	{
		my $all_adjs = $self->get_all_adjacent_gene_pairs($opts{adjacent});
		my %adj_gene_lefts = map {$_->{up_gene}{left}=>1, $_->{down_gene}{left}=>1} @$all_adjs;
		@$genes = grep { !$adj_gene_lefts{ $_->{left} } } @$genes;
	}
	
	# Filter out overlaps if needed
	if( $opts{overlap} )
	{
		@$genes = grep { !$_->{overlap} } @$genes;
	}
	
	return $genes;
}

###
# If the given coord located inside a gene -- the gene is returned, otherwise undef
sub is_the_coord_inside_gene
{
	my $self = shift;
	my($coord) = @_;
	
	foreach my $g ( @{$self->{genes}} )
	{
		return $g if $g->{left} < $coord && $coord < $g->{right};
	}
	
	return undef;
}

sub get_closest_downstream_gene_for_coord
{
	my $self = shift;
	my($coord, $strand) = @_;
	
	my $start_name = $strand eq '+' ? 'left' : 'right';
	my($min_dist,$gene) = (10000000, undef);
	foreach my $g ( @{$self->{genes}} )
	{
		next if $g->{strand} ne $strand;
		my $dist = abs( $coord - $g->{$start_name} );
		if(!$gene || $dist < $min_dist)
		{
			$min_dist = $dist;
			$gene     = $g;
		}
	}
	
	return {gene => $gene, dist => $min_dist};
}

sub get_closest_downstream_gene
{
	my $self = shift;
	my($up_gene) = @_;
	
	# We assume that $self->{genes} are sorted by their left coordinate
	if( $up_gene->{strand} eq '+' )      # Look foreard
	{
		for(my $i = $up_gene->{id}; $i < $self->num_genes(); $i++)
		{
			if( $self->{genes}[$i]{strand} eq $up_gene->{strand} )
			{
				my $gap_len = $self->{genes}[$i]{left} - $up_gene->{right} - 1;
				next if $gap_len <= 0;   # We don't consider gene overlaps
				return($self->{genes}[$i], $gap_len);
			}
		}
	}
	elsif( $up_gene->{strand} eq '-' )   # Look backward
	{
		for(my $i = $up_gene->{id}-2; $i >= 0; $i--)
		{
			if( $self->{genes}[$i]{strand} eq $up_gene->{strand} )
			{
				my $gap_len = $up_gene->{left} - $self->{genes}[$i]{right} - 1;
				next if $gap_len <= 0;   # We don't consider gene overlaps
				return($self->{genes}[$i], $gap_len);
			}
		}
	}
	else
	{
		confess("Unknown strand: ".Data::Dumper->Dump([$up_gene], ['$up_gene']));
	}

	return (undef, undef);
}

sub print_genes
{
	my $self = shift;
	print "left\tright\tf1\tf2\tf3\tf4\tf5\tf6\n";
	foreach my $g ( @{$self->get_all_genes} )
	{
		print join("\t",
			$g->{left},
			$g->{right},
			$g->{cod_p}{f1},
			$g->{cod_p}{f2},
			$g->{cod_p}{f3},
			$g->{cod_p}{f4},
			$g->{cod_p}{f5},
			$g->{cod_p}{f6}
		)."\n";
	}
}

###
# SUBROUTINES
sub sum
{
	my $res = 0;
	$res += $_ foreach @_;
	return $res;
}

sub aver
{
	return @_ == 0 ? undef : sum(@_)/@_;
}

sub revcomp
{
	my($seq) = @_;
	$seq = reverse $seq;
	$seq =~ tr/ACGTacgt/TGCAtgca/;
	return $seq;
}

1;

