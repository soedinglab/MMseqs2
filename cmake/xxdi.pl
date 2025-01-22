#!/usr/bin/env perl

use strict;
use warnings;

my $file = shift;
open my $input, '<', $file or die "Can't open file for read: $file $!";
my $text = do { local $/; <$input> };
close $input;

my @hex_values = map { "0x$_" } unpack("(H2)*", $text);
my $hex_data = join(",", map { ($_ % 16 == 0 ? "\n\t" : "") . $hex_values[$_] } 0 .. $#hex_values);
my $len_data = length($text);

my $varname = $file;
$varname =~ s/[\/.]/_/g;
print "unsigned char $varname\[\] = { $hex_data \n};\n";
print "unsigned int ${varname}_len = $len_data;\n";
