#!/usr/bin/env perl
use strict;
use warnings;
use Getopt::Long;
use constant USAGE=><<EOH;

SYNOPSIS:

perl $0 --input my.fa [Options]
Version: v20250320

Requirements:
    Programs: bgzip in HTSlib
    Modules: Getopt::Long

Descriptions:
    Determine the insert size given pairs of seqing data by
    mapping them to a reference.

Options:
    --help | -h
        Print this help/usage;
    --input | -i
        Input fasta, support bgzipped
    --transcript | -s <first/longest>
        Use first/longest transcript, default: first
    --gzip | -z
        Gziped output
    --verbose
        Detailed output for trouble-shooting;
    --version | -v!
        Print current SCRIPT version;
    --output | -o
        Output fasta file path

Example:
    perl $0

Author:
    Fu-Hao Lu
    Professor, PhD
    State Key Labortory of Crop Stress Adaptation and Improvement
    College of Life Science
    Jinming Campus, Henan University
    Kaifeng 475004, P.R.China
    E-mail: lufuhao\@henu.edu.cn
EOH
###HELP ends#########################################################
die USAGE unless @ARGV;



###Receving parameter################################################
my ($help, $version);
my (@input, $output);
my $verbose=0;
my $debug=0;
my $trans="first";
my $out_gzip=0;

GetOptions(
    "help|h!" => \$help,
    "input|i:s" => \@input,
    "transcript|s:s" => \$trans,
    "gzip|z!" => \$out_gzip,
    "output|o:s" => \$output,
    "debug!" => \$debug,
    "verbose!" => \$verbose,
    "version|v!" => \$version) or die USAGE;
($help or $version) and die USAGE;

### 检查是否提供了输出文件路径
if (!defined $output && @input > 1) {
    die "Error: When multiple input files are provided, --output option is required.\n";
}

### Defaults ########################################################
if ($trans=~/^[fF][iI][rR][sS][tT]$/) {
    $trans="first";
}
elsif ($trans=~/^[lL][oO][nG][gG][eE][sS][tT]$/) {
    $trans="longest";
}
else {
    die "Error: unknown options --transcript/-s: $trans\n";
}

### input and output ################################################
chomp $trans;

### Main ############################################################
foreach my $indf (@input) {
    my $outtmp;
    if (defined $output) {
        $outtmp = $output;
        if (@input > 1) {
            # 如果有多个输入文件，在输出文件名中添加输入文件名的标识
            (my $input_base = $indf) =~ s/^.*\///;
            $outtmp =~ s/(\.\w+)$/.$input_base$1/;
        }
    } else {
        ($outtmp=$indf)=~s/^.*\///;
    }
    print "\n\n\n######\n";
    print "Info: input: $indf\n";
    if ($indf=~/\.[gG][zZ]$/) {
        open(FASTAIN, "bgzip -dc $indf | ") || die "Error: can not open bgzipped FASTA: $indf\n";
    }
    elsif ($indf=~/\.[fF][aA]$/ or $indf=~/\.[fF][aA][sS]$/ or $indf=~/\.[fF][aA][sS][tT][aA]$/) {
        open(FASTAIN, "<", $indf) || die "Error: can not open FASTA: $indf\n";
    }
    else {
        open(FASTAIN, "<", $indf) || die "Error: can not guess input format: $indf\n";
    }
    my $num_t=0;
    my $num_c=0;
    my %seqid=();
    my $seqlen=0;
    my $seqname="";
    while (my $line=<FASTAIN>) {
        chomp $line;
        if ($line=~/^>/) {
            if ($num_t != 0) {
                if ($seqname=~/^(.*)\.(\d+)$/) {
                    my $geneid=$1;
                    my $isoform=$2;
                    if (exists $seqid{$geneid} and exists $seqid{$geneid}{$isoform}) {
                        die "Error: duplicated seqID: $seqname\n";
                    }
                    $seqid{$geneid}{$isoform}=$seqlen;
                }
                else {
                    die "Error: unknown transcript name: $seqname\n";
                }
                $seqlen=0;
            }
            ($seqname=$line)=~s/^>//;
            $seqname=~s/\s+.*$//;
            $num_t++;
        }
        else {
            $seqlen+=length($line);
        }
    }
    if ($seqname=~/^(.*)\.(\d+)$/) {
        my $geneid=$1;
        my $isoform=$2;
        if (exists $seqid{$geneid} and exists $seqid{$geneid}{$isoform}) {
            die "Error: duplicated seqID: $seqname\n";
        }
        $seqid{$geneid}{$isoform}=$seqlen;
    }
    else {
        die "Error: unknown transcript name: $seqname\n";
    }
    close FASTAIN;
    ### QC
    print "Total seqs: $num_t\n";
    print "     Genes: ".scalar(keys %seqid)."\n";
    ### select seqid and output to a list
    my %select=();
    foreach my $geneid (sort keys %seqid) {
        my @tsarr=sort {$a<=>$b} (keys %{$seqid{$geneid}});
        if ($debug or $verbose) {
            print "Debug: gene $geneid isoform: ".join("\t",@tsarr)."\n";
        }
        my $tsid=$tsarr[0];
        if (($trans eq "longest") and scalar(@tsarr)>1) {
            my $maxlen=$seqid{$geneid}{$tsarr[0]};
            if ($debug or $verbose) {
                print "Debug: seq: $geneid.$tsid length: $maxlen\n";
            }
            for (my $x=1; $x<scalar(@tsarr); $x++) {
                if ($debug or $verbose) {
                    print "Debug: seq: $geneid.$tsarr[$x] length: $seqid{$geneid}{$tsarr[$x]}\n";
                }
                if ($seqid{$geneid}{$tsarr[$x]} > $seqid{$geneid}{$tsarr[0]}) {
                    $tsid=$tsarr[$x];
                    $maxlen=$seqid{$geneid}{$tsarr[$x]};
                }
            }
        }
        if ($debug or $verbose) {
            print "Debug: selected seq: $geneid.$tsid length: $seqid{$geneid}{$tsid}\n";
        }
        $select{$geneid.".".$tsid}++;
    }
    print "  selected: ".scalar(keys %select)."\n";
    %seqid=();
    ### output selected sequences
    if ($indf=~/\.[gG][zZ]$/) {
        open(FASTAIN2, "bgzip -dc $indf | ") || die "Error: can not open bgzipped FASTA: $indf\n";
        if (!defined $output) {
            $outtmp=~s/\.(\w+)\.[gG][zZ]$/.primary.$1/;
        }
    }
    elsif ($indf=~/\.[fF][aA]$/ or $indf=~/\.[fF][aA][sS]$/ or $indf=~/\.[fF][aA][sS][tT][aA]$/) {
        open(FASTAIN2, "<", $indf) || die "Error: can not open FASTA: $indf\n";
        if (!defined $output) {
            $outtmp=~s/\.(\w+)$/.primary.$1/;
        }
    }
    else {
        open(FASTAIN2, "<", $indf) || die "Error: can not guess input format: $indf\n";
    }
    if ($out_gzip==0) {
        open(FASTAOUT, " > ", $outtmp) || die "Error: can not write FASTA: $outtmp\n";
    }
    else {
        $outtmp.=".gz";
        open(FASTAOUT, " | bgzip > $outtmp") || die "Error: can not open bgzipped FASTA: $outtmp\n";
    }
    print "\nInfo: output: $outtmp\n";
    my $out_switch=0;
    my $outcount=0;
    while (my $line=<FASTAIN2>) {
        chomp $line;
        if ($line=~/^>/) {
            (my $seqname=$line)=~s/\s+.*$//;
            $seqname=~s/^>//;
            if (exists $select{$seqname}) {
                $out_switch=1;
                $outcount++;
            }
            else {
                $out_switch=0;
            }
        }
        if ($out_switch==1) {
            print FASTAOUT $line."\n";
        }
    }
    close FASTAIN2;
    close FASTAOUT;
    print "   num out: $outcount\n";
}

print "\nInfo: Gracefully done\n";    
