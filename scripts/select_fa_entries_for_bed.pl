#!/opt/bin/env perl

sub print_usage{
	print STDERR "parameter missing\n";
	print STDERR "fastafile bedfile\n";
	exit;
}

my $fasta = "";
my $bed = "";

if ($ARGV[0]){ $fasta = $ARGV[0];
}else{ print_usage(); }
if ($ARGV[1]){ $bed = $ARGV[1];
}else{ print_usage(); }


my %contigs = ();
open FHD, $bed;
while(<FHD>){
	chomp;
	my @line = split('\s');
	$contigs{$line[0]} = 1;
}
close FHD;
print STDERR scalar(keys(%contigs));
print STDERR "\n";
open FHD, $fasta;
while(<FHD>){
	if (/^>(\S+)/){
		if ($contigs{$1}){ 
			print $_;
			$p = 1;
		}else {$p = 0; }
	}else{
		if ($p){ print $_; }
	}
}
close FHD;

