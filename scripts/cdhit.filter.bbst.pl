#!/usr/bin/env perl
use warnings;
use strict;
use constant USAGE =><<EOH;

Usage: perl  cdhit.clstr st.output

input if repeat sequences from SS, St, BB

output repeats shared by BB and St

Version: 20251224

EOH
die USAGE if (scalar(@ARGV)!=2 or $ARGV[0] eq '-h' or $ARGV[0] eq '--help');

###input
my $cdhit_input=$ARGV[0];
my $output=$ARGV[1];
my %hash=();
my $clstrnum=-1;

open(INPUT, "<", $cdhit_input) || die "Error: can not open cdhit input\n";
while (my $line=<INPUT>) {
	chomp $line;
	if ($line=~/^>/) {
		if ($line=~/^>Cluster (\d+)$/) {
			$clstrnum=$1;
		}
		else {
			die "Error: unknown cluster number\n";
		}
	}
	else {
		$line=~s/^\d+\s+\d+aa,\s+>//;
		$line=~s/:\d+:\d+.*\.\.\.\s+.*$//;
		if ($line=~/^Chr\d+B$/) {### Chr1B
			$hash{$clstrnum}{'B'}++;
		}
		elsif ($line=~/^Chr\d+$/) {
			$hash{$clstrnum}{'S'}++;
		}
		elsif ($line=~/^chr\d+St$/) {
			$hash{$clstrnum}{'St'}++;
		}
		else {
			die "Error: unknown chromosome: $line\n";
		}
	}
}
close INPUT;

$clstrnum=-1;
my $clstrprt=0;

open(INPUT, "<", $cdhit_input) || die "Error: can not open cdhit input2\n";
open(OUTPUT, ">", $output) || die "Error: can not write output\n";
while (my $line=<INPUT>) {
	chomp $line;
	if ($line=~/^>/) {
		if ($line=~/^>Cluster (\d+)$/) {
			$clstrnum=$1;
		}
		if (exists $hash{$clstrnum}{'St'} and exists $hash{$clstrnum}{'B'} and ! exists $hash{$clstrnum}{'S'}) {
			$clstrprt=1;
		}
		else {
			$clstrprt=0;
		}
	}
	if ($clstrprt) {
		print OUTPUT $line."\n";
	}
}
close INPUT;
close OUTPUT;
print "Info: Done\n";
