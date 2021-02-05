#!/usr/bin/env perl

#!/usr/bin/env perl

sub print_usage{
	print STDERR "parameters needed:\n";
	print STDERR "genome.fasta out_name reads [reads2 reads3]\n";
	print STDERR "\tif more than one reads file is given the following order is assumed:\n";
	print STDERR "\tpaired_end1 paired_end2 single_reads\n";
	exit;
}

# collect parameter
my $genome_file = "";
my $out_name = "";
my $reads_1 = "";
my $reads_2 = "";
my $reads_3 = "";

if ($ARGV[0]){$genome_file = $ARGV[0];
}else{print_usage();}
if ($ARGV[1]){$out_name = $ARGV[1];
}else{print_usage();}
if ($ARGV[2]){$reads_1 = $ARGV[2];
}else{print_usage();}
if ($ARGV[3]){$reads_2 = $ARGV[3];}
if ($ARGV[4]){$reads_3 = $ARGV[4];}


my $index = 0;
if ($genome_file =~ /index/){ $index = 1; }

# print progress
print localtime().":\tstarting to map reads to genome\n";
if ($index){ print "index:\t$genome_file\n";
}else{ print "genome:\t$genome_file\n";}
print "reads:\t$reads_1\n";
if ($reads_2 ne ""){ print "\t$reads_2\n";}
if ($reads_3 ne ""){ print "\t$reads_3\n";}


my $bt2_dir = "/opt/bin/bowtie2-2.0.5";
my $samtool_dir = "/opt/bin/samtools-0.1.18";
my $bedtool_dir = "/opt/bin/bedtools-2.17.0/bin";
my $script_dir = "/data/vnx/analysis/analysis_scripts";

# index preparations
my $gf_name = $genome_file;
if($genome_file =~ /\/([\w\.-]+)$/){ $gf_name = $1; }

unless($index){	unless (-e "bt2_index"){ `mkdir bt2_index`; }}

my $res = 0;
# create index
my $index_name = "bt2_index/".$gf_name."_index";
if ($index){ $index_name = $genome_file;
}else{
	print localtime().":\tcreate index..\n";
	$res = `$bt2_dir/bowtie2-build $genome_file $index_name`;
	if ($? == -1){ print_err("create index",$res);	}
}

# create result directory
my $out_dir = $out_name."_dir";

print localtime().":\tcreate result directory $out_dir\n";
unless (-e $out_dir){ `mkdir $out_dir`; }

# map reads

print localtime().":\tmapping reads..\n";

my $params = "";
if ($reads_2 ne ""){ 
	$params = "-1 $reads_1 -2 $reads_2";
	if ($reads_3 ne ""){ $params .= " -U $reads_3"; }
}else{ $params = "-U $reads_1";}

#print "$bt2_dir/bowtie2 -x $index_name -q $params -S $out_dir/$out_name.sam\n";

$res = `$bt2_dir/bowtie2 -x $index_name -q $params -S $out_dir/$out_name.sam`;
if ($? == -1){ print_err("map reads", $res);
}else{ print localtime().":\tdone\n"; }


# reduce sam to only hits...
open FHD, "$out_dir/$out_name.sam";
open FW, ">".$out_dir."/".$out_name.".red.sam";

while(<FHD>){
	if(/^@/){
		print FW $_;
	}else{
		my @l = split("\t");
		if($l[2] ne "*"){
			print FW $_;
		}
	}
}
close FW;
close FHD;

# postprocess

print localtime().":\tstart postprocessing..\n";

my $os = $out_dir."/".$out_name."_sorted";
my $om = $out_dir."/".$out_name."_merged";

print localtime().":\tstart sam to bam..\n";
$res = `$samtool_dir/samtools view -b -S $out_dir/$out_name.red.sam > $out_dir/$out_name.bam`;
if ($? == -1){ print_err("sam to bam",$res); }

print localtime().":\tstart sort bam..\n";
$res = `$samtool_dir/samtools sort $out_dir/$out_name.bam $os`;
if ($? == -1){ print_err("sort bam",$res); }

print localtime().":\tremove old bam..\n";
$res = `rm $out_dir/$out_name.bam`;
if ($? == -1){ print_err("remove old bam",$res); }

print localtime().":\tstart index bam..\n";
$res = `$samtool_dir/samtools index $os.bam $os.bai`;
if ($? == -1){ print_err("index bam",$res); }

print localtime().":\tstart bam to bed..\n";
$res =`$bedtool_dir/bamToBed -i $os.bam > $os.bed`;
if ($? == -1){ print_err("bam to bed",$res); }

print localtime().":\tstart merging bed to info..\n";
$res = `$bedtool_dir/mergeBed -n -nms -scores mean -i $os.bed > $om.info`;
if ($? == -1){ print_err("merge bed to info",$res); }

print localtime().":\tstart add length to info -> info_len..\n";
$res = `perl $script_dir/add_length_to_bed.pl $om.info > $om.info_len`;
if ($? == -1){ print_err("add length to info",$res); }

print localtime().":\tfinished\n";
print localtime().":\tall results are located in $out_dir\n";





#filter

# parameter
my $file = "$om.info_len";
my $col = 5;
my $min = 5;



my $sum_r = 0;
my $max_r = "";
my $min_r = "";
my $total_r = 0;
my $sum_l = 0;
my $max_l = "";
my $min_l = "";
my $total_l = 0;


my %readsupport = ();
my %min5reads = ();

open FHD, $file;
open FW, ">".$file.".col".$col."_min".$min;
while(<FHD>){
	chomp;
	my @line = split('\t');
	
	$readsupport{$line[0]} = 1;
	
	if($line[$col] >= $min){
		print FW $_."\n";
		$min5reads{$line[0]} = 1;
			
		my $val = $line[5];

		if($max_r eq ""){$max_r = $val}
		elsif($val > $max_r){$max_r = $val; }
	
		if($min_r eq ""){$min_r = $val}
		elsif($val < $min_r){$min_r = $val; }

		$sum_r += $val;
		$total_r ++;
		
		$val = $line[6];

		if($max_l eq ""){$max_l = $val}
		elsif($val > $max_l){$max_l = $val; }
	
		if($min_l eq ""){$min_l = $val}
		elsif($val < $min_l){$min_l = $val; }

		$sum_l += $val;
		$total_l ++;
				
	}
}
close FW;
close FHD;

print "regions: $total_r\t$total_l\n";
print "reads\n";
print "Sum:\t$sum_r\n";
print "Min:\t$min_r\n";
print "Max:\t$max_r\n";
print "Avg:\t";
print sprintf("%.2f", $sum_r/$total_r);
print "\n";

print "length\n";
print "Sum:\t$sum_l\n";
print "Min:\t$min_l\n";
print "Max:\t$max_l\n";
print "Avg:\t";
print sprintf("%.2f", $sum_l/$total_l);
print "\n";

open FHD, $genome_file;
open FW_RS, ">".$genome_file.".readsupport";
open FW_5R, ">".$genome_file.".min5reads";

my $p_rs = 0;
my $p_5r = 0;

my $countRS = 0;
my $count5R = 0;
while(<FHD>){
	if (/^>(\S+)/){
		if ($readsupport{$1}){ 
			print FW_RS $_;
			$p_rs = 1;
			$countRS ++;
			if ($min5reads{$1}){
				print FW_5R $_;
				$p_5r = 1;
				$count5R ++;
			}else{ $p_5r = 0; }
		}else { 
			$p_rs = 0;
			$p_5r = 0;
		}
	}else{
		if ($p_rs){ print FW_RS $_; }
		if ($p_5r){ print FW_5R $_; }
	}
}

close FW_RS;
close FW_5R;
close FHD;

print "filter stat:\n";
print "\treadsupport:\t$countRS\n";
print "\tmin5reads:\t$count5R\n";


sub print_err{
	my $mes = shift;
	my $res = shift;

	print STDERR localtime();
	print STDERR ":\tcould not $mes, exiting with output:\n";
	print STDERR $res."\n";
	exit;
}
