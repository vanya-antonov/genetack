package GeneTack::Lib;

use strict;
use warnings;

# $Id$

use Data::Dumper;
use Cwd qw(abs_path);

use Exporter;
use vars qw(@ISA @EXPORT_OK);
@ISA       = qw(Exporter);
@EXPORT_OK = qw(
	create_tmp_dir
	get_gc_content
	read_fasta
	translate_dna_to_prot
	write_fasta
);


###
# CONSTANTS
our %GENETIC_CODE = (
	'TCA'=>'S','TCC'=>'S','TCG'=>'S','TCT'=>'S',
	'TTC'=>'F','TTT'=>'F','TTA'=>'L','TTG'=>'L',
	'TAC'=>'Y','TAT'=>'Y','TAA'=>'*','TAG'=>'*',
	'TGC'=>'C','TGT'=>'C','TGA'=>'*','TGG'=>'W',
	'CTA'=>'L','CTC'=>'L','CTG'=>'L','CTT'=>'L',
	'CCA'=>'P','CCC'=>'P','CCG'=>'P','CCT'=>'P',
	'CAC'=>'H','CAT'=>'H','CAA'=>'Q','CAG'=>'Q',
	'CGA'=>'R','CGC'=>'R','CGG'=>'R','CGT'=>'R',
	'ATA'=>'I','ATC'=>'I','ATT'=>'I','ATG'=>'M',
	'ACA'=>'T','ACC'=>'T','ACG'=>'T','ACT'=>'T',
	'AAC'=>'N','AAT'=>'N','AAA'=>'K','AAG'=>'K',
	'AGC'=>'S','AGT'=>'S','AGA'=>'R','AGG'=>'R',
	'GTA'=>'V','GTC'=>'V','GTG'=>'V','GTT'=>'V',
	'GCA'=>'A','GCC'=>'A','GCG'=>'A','GCT'=>'A',
	'GAC'=>'D','GAT'=>'D','GAA'=>'E','GAG'=>'E',
	'GGA'=>'G','GGC'=>'G','GGG'=>'G','GGT'=>'G'
);


###
# SUBROUTINES

sub translate_dna_to_prot
{
	my($dna_seq) = @_;
	
	my $protein = '';
	for(my $i=0; $i<(length($dna_seq)-2); $i+=3)
	{
		my $codon   = substr($dna_seq,$i,3);
		my $u_codon = uc $codon;
		if( $GENETIC_CODE{$u_codon} )
		{
			my $aa = $GENETIC_CODE{$u_codon};
			$protein .= $codon =~ /[a-z]/ ? lc($aa) : $aa;
		}
		else
		{
#			confess("Unknown codon: $codon");
			warn("Unknown codon: $codon");
			return undef;
		}
	}
	
	return $protein;
}

sub create_tmp_dir
{
	my(%opts) = @_;
	$opts{root_dir} ||= '.';
	$opts{prefix}   ||= '';
	
	my $dir;
	while(1)
	{
		$dir = $opts{root_dir}.'/'.$opts{prefix}.int(rand(1000000));
		last unless -d $dir;
	}
	
	mkdir($dir);
	system("chmod 0777 $dir");
	
	return abs_path($dir);
}

###
# Returns array of hashes where:
#    seqs[i]
#       |
#       |--->{seqname}
#       |--->{fullname}
#       |--->{seq}
#
sub read_fasta
{
	my( $fn ) = @_;
	my( $refID, %Seqs );

	open(F, '<', $fn) or die "Can't open file '$fn'";
	while( <F> )
	{
		if( /^>(\S+)/ ){
			$refID = $1;
			next;

		}elsif( $refID ){
			s/\s+//g;
			next if /^$/;

			$Seqs{ $refID } .= $_;
		}
	}

	close F;

	return ( $refID && scalar( keys %Seqs )) ? \%Seqs : undef;
}


# INPUT: \$seq --- ref to nucleotide sequence
sub get_gc_content
{
	my( $seq ) = @_;
	my $num_gc = $$seq =~ tr/GCgc/GCgc/;
	return $num_gc / length( $$seq );
}

#################################
# Arguments:
#    $fn   - output file name ("test.fasta", for example)
#    $seqs - (optional argument) reference to array of hashes
#    seqs[i]
#       |
#       |--->{fullname}  -  (optional)
#       |--->{seqname}
#       |--->{seq}
#
sub write_fasta
{
	my($fn, $seqs, %opts) = @_;
	
	open(FASTA, ">", $fn) or die "Can't open file $fn to write: $!";
	
	my $seq_str = '';
	foreach my $h ( @$seqs )
	{
		my $head = $h->{fullname} || $h->{seqname};
		my $seq  = $h->{seq};
		
		if( $opts{width} )
		{
			$seq =~ s/(\w{$opts{width}})/$1\n/g ;
			$seq =~ s/[\n\r]+$//g;
		}
		
		print FASTA ">$head\n$seq\n";
	}
	
	close FASTA;
}

1;

