#!/usr/bin/env perl

sub print_usage{
	print STDERR "extends region by xbp in both directions\n";
	print STDERR "file bp\n";
	exit;
}

my $file = "";
my $bp = 0;


if($ARGV[0]){ $file = $ARGV[0];
}else{print_usage();}
if($ARGV[1]){ $bp = $ARGV[1];
}else{print_usage();}

open FHD, $file;
while (<FHD>){
	chomp;
	my @l = split('\t');
	print $l[0]."\t";
	my $start = $l[1] - $bp;
	if ($start < 1){ $start = 1; }
	print $start."\t";
	my $stop = $l[2] + $bp;
	print $stop."\n";
}
close FHD;


