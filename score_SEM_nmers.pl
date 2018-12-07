#!/usr/bin/perl

use strict;

my @char = ('A','C','G','T');

my $motif_file = shift(@ARGV);
my $seq = shift(@ARGV);

open(INF, "$motif_file") or die("Could not open file $motif_file\n");

my $header = <INF>;
my @mat;
while(my $line = <INF>) {
        chomp($line);
        my @wt = split("\t", $line);
        shift(@wt);
        push(@mat, \@wt);
}
close(INF);

if(length($seq) <= $#mat) {
	print stderr "Sequence shorter than motif!\n";
	exit(0);
}

print get_sem_sum_score($seq, \@mat, \@char);

#####################

sub get_sem_sum_score {
	my($seq, $mat_ref, $char_ref) = @_;
	my @NUM_SEQ = char_to_num($seq, $char_ref);
	my @mat = @{$mat_ref};

	my $max_score = -10000000;

	for(my $k = 0; $k < 2; $k++) {
		for(my $i=0;$i <= $#NUM_SEQ - $#mat + 1; $i++) {
		 	my $bit_score = 0;
			for(my $j = 0; $j <= $#mat; $j++) {
				$bit_score += $mat[$j]->[$NUM_SEQ[$i+$j]]; 
  			}
 			#print join("", @NUM_SEQ) . " $bit_score\n";
			if($bit_score > $max_score) {
				$max_score = $bit_score;
			}
		}
		@NUM_SEQ = char_to_num(revdnacomp($seq), $char_ref);
	}

	return $max_score;
}

sub char_to_num{
 my $sequence = shift;
 my $array_ref = shift;
 my @seq_array = split('', $sequence);
 my @numeric_sequence = ();
 foreach(@seq_array) {
  if($_ eq $array_ref->[0]) {
	push(@numeric_sequence, 0);
  } elsif($_ eq $array_ref->[1]) {
	push(@numeric_sequence, 1);
  }elsif($_ eq $array_ref->[2]) {
	push(@numeric_sequence, 2);
  }elsif($_ eq $array_ref->[3]) {
	push(@numeric_sequence, 3);
  }else {
	print stderr "Unkown character!\n";
	exit(0);
  }
 }
  return(@numeric_sequence);
}

sub revdnacomp {
  # my $dna = @_;  
  # the above means $dna gets the number of 
  # arguments in @_, since it's a scalar context!

  my $dna = shift; # or   my $dna = shift @_;
  # ah, scalar context of scalar gives expected results.
  # my ($dna) = @_; # would work, too

  my $revcomp = reverse($dna);
  $revcomp =~ tr/ACGTacgt/TGCAtgca/;
  return $revcomp;
}

