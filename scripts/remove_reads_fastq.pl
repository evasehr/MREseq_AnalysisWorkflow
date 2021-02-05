#!/opt/bin/env perl

sub print_usage{
	print STDERR "parameter missing\n";
	print STDERR "bed fastq [col-index for bed; default:3]\n";
	exit;
}

my $file = "",
my $reads = "";
my $col = 3;

if ($ARGV[0]){ $file = $ARGV[0]; }
else{ print_usage(); }
if ($ARGV[1]){ $reads = $ARGV[1]; }
else{ print_usage(); }
if ($ARGV[2]){ $col = $ARGV[2]; }

my %toremove = ();

open FHD, $file;
while(<FHD>){
	chomp;
	my @line = split('\t');
	if ($line[$col]){
		$toremove{$line[$col]} = 1;
	}
}
close FHD;
print localtime().": reads to remove: ";
print scalar(keys(%toremove));
print "\n";

print localtime().": removing reads from file: $reads\n";

open FHD, $reads;
open FW, ">$reads.rem.fastq";
my $c = 0;
my $i = 4;
my $print = 1;
while(<FHD>){
	unless ($_ eq "\n"){
		if (/^@(\S+)/){  #@id
			my $id = $1;
			if ($i >= 4){
				$i = 0;
				$print = 1;
				if ($toremove{$id}){ $print = 0;}
			}
			$c++;
			if ($c > 500000){ 
				print ".";
				$c = 0;
			}
		}

		if ($print){ print FW $_; }
		$i++;
	}
}
print "\n";
close FW;
close FHD;
