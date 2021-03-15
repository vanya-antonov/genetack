#!/usr/bin/perl -w

use strict;
use warnings;

my $VERSION = '1.0';

# $Id$

###
# Alessandro Gorohovski (an.gorohovski@gmail.com) & Ivan Antonov (vanya.antonov@gmail.com)
#

$|++; # Turn off buffering

use Test::More 'no_plan';
use Data::Dumper;

my $input_fna        = 'two_contigs.fna';
my $save_fsgene_seqs = 'fsgene_seqs.fna';
my $save_fsprot_seqs = 'fsprot_seqs.faa';

(my $log_file = $input_fna ) =~s/fna$/log/;

`rm -f $save_fsgene_seqs $save_fsprot_seqs $log_file`;
`../genetack_gm.pl --save_fsgene_seqs $save_fsgene_seqs --save_fsprot_seqs $save_fsprot_seqs ../examples/$input_fna > $log_file`;

# Read test file
my $tdata = &read_file("data/$log_file");

# Read log file (results)
my $rdata = &read_file( $log_file );

is_deeply( $tdata, $rdata, 'Check multiple contigs FASTA');

`rm -f $save_fsgene_seqs $save_fsprot_seqs $log_file`;
exit;


sub read_file {
	my( $in_file ) = @_;
	my %dd;

	open LOGF, $in_file or die "ERROR: '$in_file' file is missing: $!";
	while(<LOGF>){
		s/^\s+|\s+$//g;
		next if /^$/;
=comment
Num	Seq_ID	FS_coord	FS_coord_adj	FS_type	Strand	Gene_Left	Gene_Right	Fragment_Left	Fragment_Right
1	NZ_JOBF01000001.1	61593	61591	-1	+	60547	61746	60277	63023
=cut
		next unless /^\d+\t\D/;

		my( $Num, $Seq_ID, $FS_coord, $FS_coord_adj, @d ) = split /\t/;
		my $k = "$Seq_ID,$FS_coord,$FS_coord_adj";
		push @{ $dd{ $k } }, @d;
	}
	close LOGF;

	return \%dd;
}

