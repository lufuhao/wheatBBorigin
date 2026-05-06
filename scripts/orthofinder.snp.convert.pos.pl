#!/usr/bin/env perl
use warnings;
use strict;
use Data::Dumper qw/Dumper/;
use constant USAGE =><<EOH;

Usage: perl gp2gene bb.gff3 out.vcf

Version: 20251225

EOH
die USAGE if (scalar(@ARGV)!=3 or $ARGV[0] eq '-h' or $ARGV[0] eq '--help');



my $gp2gene=$ARGV[0];
my $bbgff=$ARGV[1];
my $outvcf=$ARGV[2];



sub readGroup () {
	my $RGfile=shift;
	my $RGhash={};
	open(RGGROUP, "<", $RGfile) || die "Error: can not open $RGfile\n";
	while (my $RGline=<RGGROUP>) {
		chomp $RGline;
		my @arr=split(/\t/,$RGline);
		die "Error: duplicated OG key: $RGline\n" if (exists ${$RGhash}{$arr[0]});
		${$RGhash}{$arr[0]}=$arr[1];
	}
	close RGGROUP;
	print "Info: Read OG ".scalar(keys(%{$RGhash}))."\n";
	return $RGhash;
}



sub readGff3 () {
	my $RGfile=shift;
	my $RGgene={};
	print "Info: read GFF3\n";
	if ($RGfile=~/\.gz$/i) {
		open (RGGFF3, "zcat $RGfile | ") || die "Error: can not open gzipped gff3: $RGfile\n";
	}
	else {
		open (RGGFF3, "<", $RGfile) || die "Error: can not open gff3: $RGfile\n";
	}
	while (my $RGline=<RGGFF3>) {
		chomp $RGline;
		next if ($RGline=~/^#/);
		my @RGarr=split(/\t/,$RGline);
		next unless ($RGarr[2] eq "exon");
		(my $RGtsid=$RGarr[8])=~s/^.*;Parent=//;
		$RGtsid=~s/;.*$//;
		unless (exists ${$RGgene}{$RGtsid}{'chr'}) {
			${$RGgene}{$RGtsid}{'chr'}=$RGarr[0];
		}
		unless (exists ${$RGgene}{$RGtsid}{'std'}) {
			${$RGgene}{$RGtsid}{'std'}=$RGarr[6];
		}
		die "Error: duplicated start: $RGline\n" if (exists ${$RGgene}{$RGtsid} and exists ${$RGgene}{$RGtsid}{'exon'} and exists ${$RGgene}{$RGtsid}{'exon'}{$RGarr[3]});
		${$RGgene}{$RGtsid}{'exon'}{$RGarr[3]}{"end"}=$RGarr[4];
		${$RGgene}{$RGtsid}{'exon'}{$RGarr[3]}{"len"}=$RGarr[4]-$RGarr[3]+1;
	}
	close RGGFF3;
	print "Info: read genes ".scalar(keys(%{$RGgene}))."\n";
#	print Dumper $RGgene; ### for test ###
	return $RGgene;
}



sub readFasta () {
	my ($RFfile, $RFid)=@_;
#	print "Info: read fasta\n";### for test ###
	open(FASTA, "<", $RFfile) || die "Error: can not open FASTA: $RFfile\n";
	my $RGseq="";
	my $RGtest=0;
	while (my $RFline=<FASTA>) {
		chomp $RFline;
		if ($RFline=~/^>/) {
			(my $RFts=$RFline)=~s/^>//;
			$RFts=~s/\s\+.*$//;
			if ($RFts eq $RFid) {
				$RGseq="";
				$RGtest=1;
			}
			else {
				$RGtest=0;
			}
		}
		else {
			if ($RGtest==1) {
				$RGseq.=$RFline;
			}
		}
	}
	close FASTA;
	die "Error: empty seq for $RFid in fasta $RFfile\n" if (length($RGseq)==0);
	return $RGseq;
}



sub detectCol () {
	my $str=shift;
	my ($col_hh,$col_st,$col_rr,$col_ee,$col_tt,$col_ss,$col_bb)=(-1,-1,-1,-1,-1,-1,-1);
	my @arr=split(/\t/,$str);
	for (my $i=9;$i<scalar(@arr); $i++) {
		if ($arr[$i]=~/^HORVU.MOREX/) {$col_hh=$i;next;}
		elsif ($arr[$i]=~/^Es\d+StG/) {$col_st=$i;next;}
		elsif ($arr[$i]=~/^SECCE/) {$col_rr=$i;next;}
		elsif ($arr[$i]=~/^GWHPABKY/) {$col_ee=$i;next;}
		elsif ($arr[$i]=~/^Ammut_EIv1\.0/) {$col_tt=$i;next;}
		elsif ($arr[$i]=~/^GWHPBFXR/) {$col_ss=$i;next;}
		elsif ($arr[$i]=~/^TraesCS/) {$col_bb=$i;next;}
		else {
			die "Error: detectCol $arr[$i] $str\n";
		}
	}
	if ($col_hh ==-1 or $col_st ==-1 or $col_rr ==-1 or $col_ee ==-1 or $col_tt ==-1 or $col_ss ==-1 or $col_bb ==-1) {
		die "Error: detectCol $col_hh,$col_st,$col_rr,$col_ee,$col_tt,$col_ss,$col_bb\n";
	}
#	print "Test: $col_hh,$col_st,$col_rr,$col_ee,$col_tt,$col_ss,$col_bb\n";### for TEST ###
	return ($col_hh,$col_st,$col_rr,$col_ee,$col_tt,$col_ss,$col_bb);
}
sub getLen () {
	my ($str,$len)=@_;
#	print "    DNA:\n".length($str)."\n$len\n";### for TEST ###
#	print $str."\n";### for TEST ###
	my $s=substr($str,0,$len);
	$s=~s/-//g;
#	print length($s)."\n";### for TEST ###
	return length($s);
}
sub getPos () {
	my ($hash, $ts, $pos)=@_;
	my $newpos=-1;
#	print "    Old pos: $pos\n"; ### for TEST ###
	if (${$hash}{$ts}{'std'} eq "+") {
#		print "Strand: +\n";### for TEST ###
		foreach my $start (sort {$a<=>$b} keys(%{${$hash}{$ts}{'exon'}})) {
			if (${$hash}{$ts}{'exon'}{$start}{'len'}<$pos) {
				$pos=$pos-${$hash}{$ts}{'exon'}{$start}{'len'};
			}
			else {
#				print "$newpos=$start+$pos-1\n";### for TEST ###
				$newpos=$start+$pos-1;
				last;
			}
		}
	}
	elsif (${$hash}{$ts}{'std'} eq "-") {
#		print "Strand: -\n";### for TEST ###
		foreach my $start (sort {$b<=>$a} keys(%{${$hash}{$ts}{'exon'}})) {
			if (${$hash}{$ts}{'exon'}{$start}{'len'}<$pos) {
				$pos=$pos-${$hash}{$ts}{'exon'}{$start}{'len'};
			}
			else {
				$newpos=${$hash}{$ts}{'exon'}{$start}{'end'}-$pos+1;
				last;
			}
		}
	}
	if ($newpos==-1) {
		print Dumper ${$hash}{$ts};
		die "Error: getPos newpos : $newpos\n";
	}
	return $newpos;
}



sub main () {
	my $og=&readGroup($gp2gene);
	my $gn=&readGff3($bbgff);
	open(OUTVCF, ">", $outvcf) || die "Error: can not write output\n";
	my $num=0;
	foreach my $gp (sort keys(%{$og})) {
		$num++;
		my $tsid=${$og}{$gp};
		print "### OG$num: $gp mRNA: $tsid\n";
		my $fafile="Single_Copy_Orthologue_Sequences_mRNA_MAFFT/".$gp."_mRNA.aligned.fa";
		my $vcffile="Single_Copy_Orthologue_Sequences_mRNA_MAFFT_VCF/".$gp."_mRNA.aligned.vcf";
		my $seq=&readFasta($fafile, $tsid);
		open(VCF, "<", $vcffile) || die "Error: can not open vcf file: $vcffile\n";
		my ($col_hh,$col_st,$col_rr,$col_ee,$col_tt,$col_ss,$col_bb)=(-1,-1,-1,-1,-1,-1,-1);
		while (my $line=<VCF>) {
			chomp $line;
			if ($line=~/^#CHROM/) {
				($col_hh,$col_st,$col_rr,$col_ee,$col_tt,$col_ss,$col_bb)=&detectCol($line);
				my @arr=split(/\t/, $line);
				print OUTVCF "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t".
				             $arr[$col_hh]."\t".
				             $arr[$col_st]."\t".
				             $arr[$col_rr]."\t".
				             $arr[$col_ee]."\t".
				             $arr[$col_tt]."\t".
				             $arr[$col_ss]."\t".
				             $arr[$col_bb]."\n";
			}
			next if ($line=~/^#/);
			my @arr2=split(/\t/, $line);
			next if ($arr2[1]==0); ### starting input ###
			next if ($arr2[4]=~/,/);### ignore complex ###
			next if ($arr2[3] eq 'N' or $arr2[3] eq 'n'); ### ignore gap region ###
			next if (length($arr2[3])>1 or length($arr2[4])>1);
			(my $ss_allele=$arr2[$col_ss])=~s/\/.*$//;
			(my $st_allele=$arr2[$col_st])=~s/\/.*$//;
#			next unless ($ss_allele eq $st_allele);
			die "Error: OG $gp: REF ALT: $line\n" unless ($arr2[3]=~/^[ACGT]{1}$/ and $arr2[4]=~/^[ACGT]{1}$/);
			my $pos1=&getLen($seq, $arr2[1]);
			die "Error: gff3 no mRNA: $tsid\n" unless (exists ${$gn}{$tsid});
			my $pos2=&getPos($gn, $tsid, $pos1);
			print OUTVCF ${$gn}{$tsid}{'chr'}."\t".$pos2."\t.\t".$arr2[3]."\t".$arr2[4]."\t.\t.\t".$arr2[7]."\t".$arr2[8]."\t".
				          $arr2[$col_hh]."\t".
				          $arr2[$col_st]."\t".
				          $arr2[$col_rr]."\t".
				          $arr2[$col_ee]."\t".
				          $arr2[$col_tt]."\t".
				          $arr2[$col_ss]."\t".
				          $arr2[$col_bb]."\n";
		}
		close VCF;
	}
	close OUTVCF;
}


&main()
