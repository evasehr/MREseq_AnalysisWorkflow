#!/opt/bin/env perl

sub print_usage{
	print STDERR "parameter missing\n";
	print STDERR "samfile fastafile\n";
	exit;
}

my $sam = "";
my $bed = "";

if ($ARGV[0]){ $sam = $ARGV[0];
}else{ print_usage(); }
if ($ARGV[1]){ $file = $ARGV[1];
}else{ print_usage(); }


my %contigs = ();
open FHD, $file;
while(<FHD>){
	chomp;
	if (/^>(\S*)/){
		$contigs{$1} = 1;
	}
}
close FHD;
print STDERR scalar(keys(%contigs));
print STDERR "\n";
open FHD, $sam;
while(<FHD>){
	if (/^@/){
		if (/\@HD/){
			print $_;
		}else{
			my @line = split('\s');
			$line[1] =~ s/^SN://;
			if ($contigs{$line[1]}){ print $_; }
		}
	}else{
		my @line = split('\s');
		if ($contigs{$line[2]}){ print $_; }
	}
}
close FHD;

