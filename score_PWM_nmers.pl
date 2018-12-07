#!/usr/bin/perl

use strict;

my @char = ('A','C','G','T');

my $motif_file = shift(@ARGV);
my $seq = shift(@ARGV);

open(INF, "$motif_file") or die("Could not open file $motif_file\n");

my $header = <INF>;
my $motif_name;
my $logflag = 0;

($motif_name) = $header=~/DE\t(.+)$/;
chomp($motif_name);
$motif_name =~ s/\t.*//g;
#print stderr $motif_name . "\n";

unless($motif_name){die "Couldn't extract motif name from PWM header.\n"}


my @mat;
while(my $line = <INF>) {
	chomp($line);
	next if($line eq "XX");
	my @wt = split /\t/, $line;
	shift(@wt);
	pop(@wt);
	push(@mat, \@wt);
}
close(INF);

if(length($seq) <= $#mat) {
	print stderr "Sequence shorter than motif!\n";
	exit(0);
}

print get_max_bit_score($seq, \@mat, \@char, $logflag);


#####################

sub get_max_bit_score {
	my($seq, $mat_ref, $char_ref, $logflag) = @_;
	my @NUM_SEQ = char_to_num($seq, $char_ref);
	my @mat = @{$mat_ref};

	my $max_score = -10000000;

	for(my $k = 0; $k < 2; $k++) {
		for(my $i=0;$i <= $#NUM_SEQ - $#mat + 1; $i++) {
			my $bit_score = 0;
			for(my $j = 0; $j <= $#mat; $j++) {
				my $sum = $mat[$j]->[0] + $mat[$j]->[1] + $mat[$j]->[2] + $mat[$j]->[3]; 
				$bit_score += log2(($mat[$j]->[$NUM_SEQ[$i+$j]] + 0.25) / ($sum + 1)) - log2(0.25); 
  			}
# 			print join("", @NUM_SEQ) . " $bit_score\n";
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

sub get_max_add{ 
my @mat = @_;

my $len = scalar(@mat);

my @maxs; 
my @mins; 

for(my $i=0;$i <= $#mat;$i++) { 
	my @list = sort {$b <=> $a} @{$mat[$i]}; 
	# max is element 0 
        my $sum = $mat[$i]->[0] + $mat[$i]->[1] + $mat[$i]->[2] + $mat[$i]->[3]; 
      	$maxs[$i] = &log2(($list[0] + 0.25) / ($sum + 1)) - &log2(0.25); 
      	$mins[$i] = &log2(($list[3] + 0.25) / ($sum + 1)) - &log2(0.25); 
}

return(\@mins,\@maxs); 
}

sub log2{ 
	my $n = shift;
	return log($n)/log(2.0);
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

