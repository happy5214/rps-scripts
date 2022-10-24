#!/usr/bin/perl

use 5.006;
use strict;
use warnings;

my $decimal_regex = qr/[+-]?(\d+\.\d+|\d+\.|\.\d+)/;

my $results_filename = shift // 'lresults.txt';
my $maxlen_filename = shift // 'maxlen.txt';

my @fftlens = ();

open(my $maxlen_file, '<', $maxlen_filename) or die "Could not open max length file: $!";
while (<$maxlen_file>) {
	chomp;
	my @pair = split;
	push @fftlens, 0 + $pair[0];
}
close($maxlen_file) or warn "Could not close max length file: $!";

my %times = ();

my $index = 0;

open(my $results_file, '<', $results_filename) or die "Could not open LLR results file: $!";
while (<$results_file>) {
	chomp;
	my $result = $_;
	my ($n) = ($result =~ /\d+\*2\^(\d+)-1 is not prime\./);
	my ($ms) = ($result =~ /Time : $decimal_regex ms\./);
	unless ($ms) {
		my ($seconds) = ($result =~ /Time : $decimal_regex sec\./);
		$ms = $seconds * 1000;
	}
	$times{$fftlens[$index]} = $ms / $n;
	$index++;
}
close($results_file) or warn "Could not close LLR results file: $!";

printf "%8s%8s\n", 'FFT len', 'Time';
for my $fftlen (@fftlens) {
	last unless $index;
	my $time = $times{$fftlen} or 0.0;
	printf "%8s    %.3f\n", $fftlen, $time;
	$index--;
}
