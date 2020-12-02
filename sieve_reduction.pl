#!/usr/bin/env perl

my $p_min = shift;
my $p_max = shift;
my $terms = shift // 1;

my $reduction = 1.0 - log($p_min)/log($p_max);
my $reduction_terms = $terms * $reduction;

print "${reduction_terms}\n";
