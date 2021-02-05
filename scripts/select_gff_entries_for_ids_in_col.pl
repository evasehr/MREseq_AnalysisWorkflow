#!/opt/bin/env perl

sub print_usage{
	print STDERR "parameter missing\n";
	print STDERR "gfffile file col\n";
	exit;
}

my $file = "";
my $bed = "";
my $col = 0;

if ($ARGV[0]){ $file = $ARGV[0];
}else{ print_usage(); }
if ($ARGV[1]){ $bed = $ARGV[1];
}else{ print_usage(); }
if ($ARGV[2]){ $col = $ARGV[2];
}


my %contigs = ();
open FHD, $bed;
while(<FHD>){
	chomp;
	my @line = split('\t');
	$contigs{$line[$col]} = 1;
}
close FHD;
print STDERR scalar(keys(%contigs));
print STDERR "\n";
open FHD, $file;
while(<FHD>){
	chomp;
	my @line = split('\t');
	if ($line[8] =~ /(\w+)$/){
		#print $1;
		if ($contigs{$1}){
			print $_."\n";
		}
	}
	if ($line[8] =~ /(\w+)\./){
		#print $1;
		if ($contigs{$1}){
			print $_."\n";
		}
	}
}
close FHD;

