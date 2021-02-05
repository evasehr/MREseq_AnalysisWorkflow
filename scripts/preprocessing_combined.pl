#!/usr/bin/env perl

use Getopt::Std;


sub print_usage{
	print STDERR "parameter missing\n";
	print STDERR "\t-f\tfile (required)\n";
	print STDERR "\t-g\tfile2 (optional)\n";
	print STDERR "\t-N\t1: removes Ns at both ends, sequences with Ns within the sequence are removes completely\n";
	print STDERR "\t-a\tadapter suffixes\n";
	print STDERR "\t-s\tadapter sequence (default is methylation adapter, metadCC will remove Psp adapater)\n";
	print STDERR "\t-q\tlow quality ends min quality (recommended 30)\n";
	print STDERR "\t-l\tshort sequences min length (recommended 50)\n";
	exit;
}

my %opts = ();

getopt('fgNasql', \%opts);

foreach my $key(%opts){
		print "$key\t";
		print $opts{$key}."\n";
}

my $file = "";
if ($opts{"f"}){ $file = $opts{"f"};
}else{ print_usage(); }

my $adapter = "CACGACTGTTTAAA";
my @fragments = ();
my $start = "";
if ($opts->{"a"}){
	if ($opts->{"s"}){
		if ($opts->{"s"} eq "metadCC"){
			$start = "GC";
			@fragments = ("ACGACTGTTTAAA","CGACTGTTTAAA","GACTGTTTAAA","ACTGTTTAAA","CTGTTTAAA","TGTTTAAA","GTTTAAA","TTTAAA","TTAAA","TAAA","AAA","AA","A");
		}else{
			$adapter = $opts->{"s"};
			@fragments = @{create_suffixes($opts{"s"})}; 
		}
	}else{
		$start = "GC";
		
		@fragments = ("ACGACTGTTTAAA","CGACTGTTTAAA","GACTGTTTAAA","ACTGTTTAAA","CTGTTTAAA","TGTTTAAA","GTTTAAA","TTTAAA","TTAAA","TAAA","AAA","AA","A");
	}
}


run($file, \%opts, \@fragments, $start, $adapter);

if ($opts{"g"}){
	run($opts{"g"}, \%opts);
	compare_files($file, $opts{"g"});
}

sub run{
	my $file = shift;
	$opts = shift;
	$fragments = shift;
	$start = shift;
	$adapter = shift;

	print localtime()." working on file $file\n";
	
	
		
	open FHD, $file;
	open FW, ">$file.Q".$opts->{"q"}.".L".$opts->{"l"};

	my $i = 4;
	my $ind = 0;
	my $p = 0;
	my %entry = ();
	
	# stats
	
	my $processed = 0;
	my $kept = 0;
	my $removed_allN = 0;
	my $removed_Ninseq = 0;
	my $removed_low_qual = 0;
	my $removed_too_short = 0;
	my $removed_too_shortN = 0;
	my $removed_too_shortA = 0;
	
	my $short_N_start = 0;
	my $short_N_end = 0;
	my $short_complete_adapter = 0;
	my $short_partial_adapter = 0;
	my $short_Q_start = 0;
	my $short_Q_end = 0;
	
	
	
	while (<FHD>){
		chomp;	
		unless ($_ eq ""){ # ignore empty lines
			if (/^@/){	#@id
				
				if ($i >= 4){
					# reinitialize
					$i = 0;
					%entry = ();
					$entry{"id"} = $_;		
					($ind, $p) = print_progress ($ind, $p);		
				}
			}
			if ($i == 1){ $entry{"seq"} = $_; }
			if ($i == 2){ $entry{"info"} = $_; }	
			if ($i == 3){ 
				$entry{"quality"} = $_; 
			
				$processed ++;
				my $rem = 0;
				if ($opts->{"N"}){
					
					my ($rt, $r, $ss, $se) = remove_Ns(\%entry);
					
					$rem += $rt + $r;
					$removed_allN += $rt;
					$removed_Ninseq += $r;
					$short_N_start += $ss;
					$short_N_end += $se;
				}
				
				if ($rem == 0){
					if ($opts->{"l"}){
						if (length($entry{"seq"}) < $opts->{"l"}){ # too short
							$rem = 1;
							$removed_too_shortN ++;
						}
					}
				}
				
				if ($rem == 0){
					if ($opts->{"a"}){
						my ($adcomp, $adpart) = remove_adapter($adapter, \%entry, \@fragments, $start);
						$short_complete_adapter += $adcomp;
						$short_partial_adapter += $adpart;
					}
				}
				
				if ($rem == 0){
					if ($opts->{"l"}){
						if (length($entry{"seq"}) < $opts->{"l"}){ # too short
							$rem = 1;
							$removed_too_shortA ++;
						}
					}
				}
				if ($rem == 0){
					if ($opts->{"q"}){
						my ($r, $ss, $se) = remove_low_quality(\%entry, $opts->{"q"});
						$rem += $r;
						$removed_low_qual += $r;
						$short_Q_start += $ss;
						$short_Q_end += $se;
					}
				}
				if ($rem == 0){
					if ($opts->{"l"}){
						if (length($entry{"seq"}) < $opts->{"l"}){ # too short
							$rem = 1;
							$removed_too_short ++;
						}
					}
				}
				if ($rem == 0){
					$kept ++;
					print FW $entry{"id"}."\n";
					print FW $entry{"seq"}."\n";
					print FW $entry{"info"}."\n";
					print FW $entry{"quality"}."\n";
				}

			}
			$i++;
		}
	}
	close FW;
	close FHD;
	print localtime()." results for file $file:\n";
	
	print "processed reads:\t$processed\n";
	print "remaining reads:\t$kept\n";
	print "\n";
	print "undefined bases N:\n";
	print "\tseq complete N (removed):\t$removed_allN\n";
	print "\tN in seq (removed):\t$removed_Ninseq\n";
	print "\tN at beginning:\t$short_N_start\n";
	print "\tn at end:\t$short_N_end\n";
	print "\tremoved due to length min".$opts->{"l"}." after N removal:\t$removed_too_shortN\n";
	print "\n";
	print "adapter:\n";
	print "\tcomplete adapter:\t$short_complete_adapter\n";
	print "\tpartial adapter:\t$short_partial_adapter\n";
	print "\tremoved due to length min".$opts->{"l"}." after adapter removal:\t$removed_too_shortA\n";
	print "\n";
	print "low quality filter Q".$opts->{"q"}.":\n";
	print "\tremoved due to general low quality:\t$removed_low_qual\n";
	print "\tlow quality at beginning:\t$short_Q_start\n";
	print "\tlow quality at end:\t$short_Q_end\n";
	print "\n";
	print "too short sequence:\n";
	print "\tremoved due to length min".$opts->{"l"}." after low quality filter:\t$removed_too_short\n";
}

sub print_progress {
	my $ind = shift;
	my $p   = shift;

	$ind++;
	if ( $ind > 50000 ) {
		print ".";
		$p++;
		$ind = 0;
	}
	if ( $p >= 20 ) {
		print "\t".localtime()."\n";
		$p = 0;
	}
	return ($ind, $p);
}


sub remove_Ns{
	$entry = shift;
		
	my $seq = $entry->{"seq"};
	my $qual = $entry->{"quality"};
	
	my @res = (0,0,0,0);
	
	if ($seq =~ /N/){
		$c1 ++;
		if ($seq =~ /^N+$/){ # all N -> remove entry
			@res = (1,0,0,0);
		}elsif ($seq =~ /^.*[^N]+N+[^N].*$/){ # internal N -> remove entry
			@res = (0,1,0,0);
		}else{
			if ($seq =~ /^(\w*[^N])(N{1,})$/){ # ends with N -> shorten seq
				$rem_end = length($2);
				$seq = "$1";
				$replaced_end ++;
				$res[3] = 1;
				
				if ($qual =~ /^(.*)(.{$rem_end})$/){ $qual = $1; }
				$entry->{"info"} .= " $rem_end Ns removed from end";
					
			}
			if($seq =~ /^(N+)([^N].*)$/){ # starts with N -> shorten seq
				$rem_start = length($1);
				$seq = "$2";
				$replaced_start ++;
				$res[2] = 1;
				if ($qual =~ /^(.{$rem_start})(.*)$/){ $qual = $2; }
				$entry->{"info"} .= " $rem_start Ns removed from beginning";
				
			}
			$entry->{"seq"} = $seq;
			$entry->{"quality"} = $qual;

		}	
	}
	return @res;
			
}

sub create_suffixes{
	my $seq = shift;
	 
	my @suffixes = ();
 
	my $i = length($seq);
	while($i >= 1){
		if ($seq =~ /(.*)(.{$i})$/){
			push (@suffixes, $2);
			$i --;
		}else{
			last;
		}
	}
	return @suffixes;
}

sub remove_adapter{
	my $complete = shift;
	$entry = shift;
	$fragments = shift;
	my $start = shift;
	
	my $seq = $entry->{"seq"};
	my $c = 0;
	my $f = 0;
		
	if (/^(.*$complete)($start.*)$/){ # whole adapter
		$entry->{"seq"} = $2;
		$to_cut = length($1);
		$c ++;

	}else{ # look for fragments
		foreach my $frag (@{$fragments}){
			if (/^(.{0,5}$frag)($start.*)$/){
				$entry->{"seq"} = $2;
				$to_cut = length($1);
				$f ++;
				last;
			}
		}
	}
				
	if ($to_cut > 0){
		$entry->{"info"} .= "due to adapter $to_cut nucleotides removed from beginning";
		
		my $tmp = $entry->{"quality"};
		if ($tmp =~ /^(.{$to_cut})(.+)$/){
			$entry->{"quality"} = $2;
		}
	}
	return ($c, $f);
}

sub remove_low_quality{
	
	$entry = shift;
	my $qual_limit = shift;
	
	my $se = 0;
	my $ss = 0;
	my $r = 0;

	my $seq = $entry->{"seq"};
	my $qual_line = $entry->{"quality"};
	my @qualities = split('', $qual_line);
	
	# look at front
	my $rem_start = 0;
	foreach my $q (@qualities){
		my $Q = ord($q) - 33;
		if ($Q < $qual_limit){
			$rem_start ++;
			$ss = 1;
		}else{
			last;
		}	
	}
	
	# look at end
	@qualities = reverse @qualities;
	my $rem_end = 0;
	foreach my $q (@qualities){
		my $Q = ord($q) - 33;
		if ($Q < $qual_limit){
			$rem_end ++;
			$se = 1;
		}else{
			last;
		}	
	}
	
	# remove from beginning (seq and quality)
	if ($rem_start > 0){
		if ($seq =~ /^(.{$rem_start})(.*)$/){
			$seq = $2;
		}
		if ($qual_line =~ /^(.{$rem_start})(.*)$/){
			$qual_line = $2
		}
		$a ++;
	}
	# remove from end (seq and quality)
	if($rem_end > 0){
		if ($seq =~ /^(.*)(.{$rem_end})$/){
			$seq = $1;
		}
		if ($qual_line =~ /^(.*)(.{$rem_end})$/){
			$qual_line = $1;
		}
		$a++;
	}
	
	if (length($seq) > 0){
		$entry->{"seq"} = $seq;
		$entry->{"quality"} = $qual_line;
		if ($se + $ss > 0){
			$entry->{"info"} .= " low quality front: $rem_start end: $rem_end";
		}
	}else{
		$r = 1;
	}
	
	return ($r, $ss, $se);
}

sub compare_files{
	my $file1 = shift;
	my $file2 = shift;

	my %content = ();
	my $i = 4;
	
	my $p = 0;
	my $ind = 0;
	
	my $entry = "";
	my $id = "";
	open FHD1, $file1;
	while (<FHD1>){
		chomp;
		unless ($_ eq ""){ # ignore empty lines
			if (/^@([\w:-]+)\s?/){	#@id
			# @HWI-ST225:505:C11VVACXX:2:1101:1416:2084 1:N:0:AGTCAA
			# @HWI-ST225:505:C11VVACXX:2:1101:1416:2084 2:N:0:AGTCAA
				if ($i >= 4){
					$i = 0;
					$id = $1;
					$entry = "$_\n";
					$content{$id} = 1; 
					
					# print progress
					($ind, $p) = print_progress($ind, $p);
				}
			}
			$i++;
		}
	}
	close FHD1;
	
	my %found = ();
	
	$p = 0;
	$ind = 0;
	$id = "";
	
	my $prefix = $file1;
	if ($file1 =~ /^(.+).fastq/){
		$prefix = $1;
	}
	open FHD2, $file2;
	open FW2, ">".$prefix."_pairs_2.fastq";
	open FW3, ">".$prefix."_single_reads.fastq";
	my $c = 1;
	while (<FHD2>){
		chomp;
		unless ($_ eq ""){ # ignore empty lines
			if (/^@([\w:-]+)\s?/){	#@id
			# @HWI-ST225:505:C11VVACXX:2:1101:1416:2084 1:N:0:AGTCAA
			# @HWI-ST225:505:C11VVACXX:2:1101:1416:2084 2:N:0:AGTCAA
			
				if ($i >= 4){
					$i = 0;
					$id = $1;
					if ($content{$id}){
						print FW2 "$_\n";
						$found{$id} = $c;
						$c ++;
					}else{ print FW3 "$_\n"; }
					
					# print progress
					($ind, $p) = print_progress($ind, $p);
				}
			}
			if ($i >= 1 && $i <= 3){
				if ($content{$id}){ print FW2 "$_\n"; 
				}else{ print FW3 "$_\n"; }
			}
			$i++;
		}
	}
	
	close FW2;
	close FHD2;
	
	open FHD, $file1;
	my $c1 = 1;
	my $pr = 0;
	open FW1, ">".$prefix."_pairs_1.fastq";
	while (<FHD>){
		chomp;
		unless ($_ eq ""){ # ignore empty lines
			if (/^@([\w:-]+)\s?/){	#@id
				if ($i >= 4){
					$i = 0;
					$id = $1;
					if ($found{$id}){
						if ($c1 == $found{$id}){
							print FW1 "$_\n";
							$pr = 1;
						}else{ print "bla"; }
						$c1 ++;
					}else{
						print FW3 "$_\n";
						$pr = 0;
					}
					
					# print progress
					($ind, $p) = print_progress($ind, $p);
				}
			}
			if ($i >= 1 && $i <= 3){
				if ($pr){ print FW1 "$_\n"; 
				}else{ print FW3 "$_\n"; }
			}
			$i++;
		}
	}
	close FW1;
	close FW3;
}

	
	
