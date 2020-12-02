#!/usr/bin/perl

my $decimalregex = qr/[+-]?(\d+\.\d+|\d+\.|\.\d+)/;

my $results_filename = shift // 'lresults.txt';

my @fftlens = ();
my $maxlen_filename = 'maxlen.txt';
open(my $maxlen_file, '<', $maxlen_filename);
while (<$maxlen_file>) {
	chomp;
	my @pair = split;
	push @fftlens, 0 + $pair[0];
}
close($maxlen_file);

my %times = ();

my $index = 0;
open(my $results_file, '<', $results_filename);
while (<$results_file>) {
	chomp;
	my $result = $_;
	my ($n) = ($result =~ /\d+\*2\^(\d+)-1 is not prime\./);
	my ($ms) = ($result =~ /Time : $decimalregex ms\./);
	unless ($ms) {
		my ($seconds) = ($result =~ /Time : $decimalregex sec\./);
		$ms = $seconds * 1000;
	}
	$times{$fftlens[$index]} = $ms / $n;
	$index++;
}
close($results_file);

printf "%8s%8s\n", 'FFT len', 'Time';
for my $fftlen (@fftlens) {
	last unless $index;
	my $time = $times{$fftlen} or 0.0;
	printf "%8s    %.3f\n", $fftlen, $time;
	$index--;
}
