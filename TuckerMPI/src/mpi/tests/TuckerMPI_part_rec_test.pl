#!/usr/bin/perl
use strict;
use warnings;

my $exitval = 0;
my $temp = $ARGV[0];

print "$temp ../drivers/bin/reconstruct --parameter-file ../tests/input_files/reconstruction/paramfile.txt\n";
system("$temp ../drivers/bin/reconstruct --parameter-file ../tests/input_files/reconstruction/paramfile.txt");
$exitval = $exitval | $?;
print "running comparison\n";
system('../../serial/compare/bin/compare 45 reconstructed.mpi input_files/reconstruction/sthosvd_reconstructed.mpi 1e-10');
$exitval = $exitval | $?;

exit $exitval;
