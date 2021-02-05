#!/usr/bin/env perl

sub print_usage{
	print STDERR "parameters needed:\n";
	print STDERR "genome.fasta/index out_name reads [reads2 reads3]\n";
	print STDERR "if more than one reads file is given the following order is assumed:\n";
	print STDERR "paired_end1 paired_end2 single_reads\n";
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
my $script_dir = "/data/AnalysisTemp2/NGS_analysis/Analysis_scripts";

# index preparations
my $gf_name = $genome_file;
if($genome_file =~ /\/([\w\.-]+)$/){ $gf_name = $1; }

unless($index){	unless (-e "bt2_index"){ `mkdir bt2_index`; }}

my $res = 0;
# create index
my $index_name = "bt2_index/".$gf_name."_index";
if ($index){ $index_name = $genome_file;
}else{
	$res = `$bt2_dir/bowtie2-build $genome_file $index_name`;
	if ($? == -1){ print_err("create index",$res); }
}

# create result directory
my $out_dir = $out_name."_dir";

print localtime().":\tcreate result directory $out_dir\n";
unless (-e $out_dir){ `mkdir $out_dir`; }

# map reads

print localtime().":\tmapping reads..\n";

my $params = "";
if ($reads_2 ne ""){ $params = "-1 $reads_1 -2 $reads_2 -U $reads_3";
}else{ $params = "-U $reads_1";}

$res = `$bt2_dir/bowtie2 -x $index_name -q $params -S $out_dir/$out_name.sam`;
if ($? == -1){ print_err("map reads", $res);
}else{ print localtime().":\tdone\n"; }



# postprocess

print localtime().":\tstart postprocessing..\n";

my $os = $out_dir."/".$out_name."_sorted";
my $om = $out_dir."/".$out_name."_merged";

print localtime().":\tstart sam to bam..\n";
$res = `$samtool_dir/samtools view -b -S $out_dir/$out_name.sam > $out_dir/$out_name.bam`;
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




sub print_err{
	my $mes = shift;
	my $res = shift;

	print STDERR localtime();
	print STDERR ":\tcould not $mes, exiting with output:\n";
	print STDERR $res."\n";
	exit;
}

