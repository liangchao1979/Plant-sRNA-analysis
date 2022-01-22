#!/usr/bin/perl
# Copyright (c)  2010-
# Program:			sRNA_dist_lib
# Author:			Gaolei <highlei@gmail.com or leigao@ucr.edu>
# Program Date:		2010.03.11
# Modifier:			Gaolei <highlei@gmail.com or leigao@ucr.edu>
# Last Modified:	2014.03.24
# Description:	the distribution of miRNA in every lib
#**************************
# Version: 1.1	use the soap results and normalize the data.
# Version: 1.2	fix bugs; use $c
# Version: 1.3	use eq replace the =~ when compare mature miRNAs with small RNAs.
# Version: 1.4	use the input library_size
# Version: 1.5	$inclusion: perfect match or inclusion
# Version: 2.0	use the mature miRNAs
# Version: 2.1	can read fasta file without format >ath1_89_90x 
# Version: 3.0	the distribution of miRNA and its truncation and tailing in every lib
#**************************
# e-mail:highlei@gmail.com

my $version="3.0";
print STDERR ("\n==================| $0 start |==========================================\n");

my $start = time();
my $Time_Start = sub_format_datetime(localtime(time())); #运行开始时间
print STDERR "Now = $Time_Start\n\n";


use Getopt::Std;
getopts("hi:d:s:e:m:b:g:a:n:c:S:N:L:O:I:A:l:0:");
my $flag0		= (defined $opt_0) ? $opt_0 : 1;
my $infile		= $opt_i;
my $libFile		= (defined $opt_d) ? $opt_d : "";
#my $startPos	= (defined $opt_s) ? $opt_s : 2;
#my $endPos		= (defined $opt_e) ? $opt_e : 13;
#my $mismatch	= (defined $opt_m) ? $opt_m : 1;	 
#my $bulge		= (defined $opt_b) ? $opt_b : 1;	# = gap
#my $GUpair		= (defined $opt_g) ? $opt_g : 0.5;
#my $double		= (defined $opt_a) ? $opt_a : 2;	# 
my $start_len	= (defined $opt_l) ? $opt_l : 8;	# max truncation length
my $soap_result	= (defined $opt_S) ? $opt_S : "";
my $normal_base	= (defined $opt_N) ? $opt_N : 1000000;
my $lib_size	= (defined $opt_s) ? $opt_s : "";	## the libraries size
my $library		= (defined $opt_L) ? $opt_L : "";	## the libraries size
my $opformat	= (defined $opt_O) ? $opt_O : 1;	# 0: +, -; 1: + plus -; 2: +; 3: -
my $inclusion	= (defined $opt_I) ? $opt_I : 0;	# 0: match; 1: include; 2: 4kind T&T; 3: all T&T
my $ask_input	= (defined $opt_A) ? $opt_A : 1;
#my $trunction	= (defined $opt_T) ? $opt_T : 0; # 0: no truncation; 1: do truncation

if ($opt_h || $infile eq ""){
	usage();
}


sub numerically{$a<=>$b};

use FileHandle;
use strict;

my ($i,$j,$k,$m,$n,$k1,$k2,$k3,$k4,$file,$line,$in,$match,$omatch,$a,$b,$end);
my (@buf,@tmp,@genome,@gnmName,@gnmLen);
my (%gnm,%seg,%num,%matureSeq);
my $key="";
my ($endLen,$minLen,$addLen,$maxLen,$foldFileNum,$bfile);
my ($foldFile,$foldFileOut);

#===========================================================================================================
#====================                  main
#===========================================================================================================
#my $flag0	= 1;
my $yesorno	= "y";
while ($flag0) {
	print STDERR ("\n------------------------------------------------------------\n");
	print STDERR ("\n $0 version $version\n\n");
	print STDERR ("Settings for this run:");
	printf STDERR ("\n i  %55s : %-25s","input miRNA file",$infile);#%45s
	printf STDERR ("\n d  %55s : %-25s","input library file(s)",$libFile);
	printf STDERR ("\n S  %55s : %-25s","input soap result file(s)",$soap_result);
	printf STDERR ("\n N  %55s : %-25s","input the base for normalization",$normal_base);
	printf STDERR ("\n s  %55s : %-25s","input library size",$lib_size);
	printf STDERR ("\n L  %55s : %-25s","input library",$library);
	printf STDERR ("\n O  %55s : %-25s","output format",$opformat);
	printf STDERR ("\n I  %55s : %-25s","perfect match or inclusion or Truncation & Tailing",$inclusion);
	printf STDERR ("\n A  %55s : %-25s","ask input or not",$ask_input);
	printf STDERR ("\n l  %55s : %-25s","input the max truncation length",$start_len);
#	if($zero==1) {printf STDERR ("\n z  %45s : %-25s","output coverage region?","1");}
#	elsif ($zero == 0) {printf STDERR ("\n z  %45s : %-25s","concise output","0");}
#	else {printf STDERR ("\n z  %45s : %-25s","output all genome","2");}
	printf STDERR ("\n x  %55s","exit the program!");
	print STDERR ("\n\n");
	print STDERR "y to accept these or type the letter for one to change!\n";
	$yesorno = <STDIN>;	$yesorno =~s/[\s|\t|\r|\n]+$//g;	$yesorno = lc($yesorno);
	if ($yesorno eq "y") {print STDERR ("\n------------------------------------------------------------\n\n\n"); $flag0 = 0;}
	elsif($yesorno eq "i") {print STDERR "please input miRNA file:\n"; $infile	= <STDIN>;	$infile	=~s/[\s|\t|\r|\n]+$//g;}
	elsif($yesorno eq "d") {print STDERR "please input library file(s):\n"; $libFile	= <STDIN>;$libFile	=~s/[\s|\t|\r|\n]+$//g;}
	elsif($yesorno eq "S") {print STDERR "please input soap result file(s):\n"; $soap_result	= <STDIN>;$soap_result	=~s/[\s|\t|\r|\n]+$//g;}
	elsif($yesorno eq "N") {print STDERR "please input the base for normalization:\n";$normal_base	= <STDIN>;$normal_base	=~s/[\s|\t|\r|\n]+$//g;}
	elsif($yesorno eq "s") {print STDERR "please input library size:\n";$lib_size	= <STDIN>;$lib_size	=~s/[\s|\t|\r|\n]+$//g;}
	elsif($yesorno eq "L") {print STDERR "please input library:\n";$library	= <STDIN>;$library	=~s/[\s|\t|\r|\n]+$//g;}
	elsif($yesorno eq "O") {print STDERR "please output format (0: +,-; 1:+ plus -; 2:+; 3:-):\n";$opformat	= <STDIN>;$opformat	=~s/[\s|\t|\r|\n]+$//g;}
	elsif($yesorno eq "I") {print STDERR "please input 0: match; 1: include; 2: 4kind T&T; 3: all T&T:\n";$inclusion	= <STDIN>;$inclusion	=~s/[\s|\t|\r|\n]+$//g;}
	elsif($yesorno eq "A") {print STDERR "please input ask iput (1) or not (0):\n";$ask_input	= <STDIN>;$ask_input	=~s/[\s|\t|\r|\n]+$//g;}
	elsif($yesorno eq "l") {print STDERR "please input the max trunction length:\n";$start_len	= <STDIN>;$start_len	=~s/[\s|\t|\r|\n]+$//g;}

	elsif($yesorno eq "x") {print STDERR ("============================================================\n");exit(0);}
}

if ($lib_size ne "") {
	@tmp = split(/\,/,$lib_size);
	print STDERR "\tlib_size: @tmp\t","\n";
	for ($k1 = 0; $k1	< @tmp ;$k1++) {
		$i	= $k1;
		$num{$k1}->[0]	= $k1;
		$num{$k1}->[1]	= $tmp[$k1];
		print "Soap $k1,\tnum=",$num{$i}->[1];	print STDERR "\tSoap $i, $k1,\tnum=",$num{$i}->[1];
		$num{$i}->[1]	= 1.0*$num{$i}->[1]/$normal_base;
		print "\t normalize=",$num{$i}->[1],"\n";	print STDERR "\t normalize=",$num{$i}->[1],"\n";
	}
}

############################################ read file ######################################################


$k1	= 0;	$k3	= 0;	$k4	= 0;
$file = new FileHandle ("$infile") || die("Cannot open miRNA file: $infile\n");
while(<$file>)
{
	$_=~s/^[\s|\t]+//g;
	$_=~s/[\s|\t|\r|\n]+$//g;
	if ($_ =~/^>(\S+)/) {
		$i	= $1;
		$gnm{$i}{"Seq"}	= "";
		$k1++;
	} else {
		$gnm{$i}{"Seq"}	.= uc($_);
	}
}
close $file || die;

print STDERR "\nNow = [",sub_format_datetime(localtime(time())),"]\tLoad file: $infile OK\t$k1\n\n";

foreach $i (keys %gnm) {
	$j = $gnm{$i}{"Seq"};
	$matureSeq{$j}{$i}	= 0;
	if ($inclusion >= 2) {
		for ($k1 = 0; $k1 <= $start_len ;$k1++) {
			$matureSeq{substr($j,0,length($j)-$k1)}{$i}	= $k1;
		}
	}
}
############################################ read soap result file ##########################################
$i	= 0;	$k3	= 0;	$m	= 0;	my $j2 = "";	my $j3	= "";
if ($soap_result ne "") {
	if ($soap_result =~/^(.+\/)([^\/]*)\*([^\*]+)$/) {								# data/a*b data/*b
		$k	= $3;	$j	= $1;	$j2	= $2;print STDERR "1.j2=$j2,k=$k,j=$j\n";	# $k = b;	$j	= data/
	} elsif ($soap_result =~/^([^\/]*)\*([^\*]+)$/) {								# a*b, *b
		$k	= $2;	$j	= "./";	$j2	= $1;print STDERR "2.j2=$j2,k=$k,j=$j\n";	# $k = b;	$j	= ./
	} elsif ($soap_result =~/^(.+\/)([^\/]+)$/) {									# data/a
		$k	= $2;	$j	= $1;	$j2 = "";print STDERR "3.j2=$j2,k=$k,j=$j\n";	# $k = a;	$j	= data/
	} elsif ($soap_result =~/^([^\*]+)$/) {											# a
		$k	= $1;	$j	= "./";	$j2 = "";print STDERR "4.j2=$j2,k=$k,j=$j\n";	# $k = a;	$j	= ./
	} else {
		print STDERR	$soap_result;
		die("reinput $soap_result!\n");
	}
	$k	= quotemeta($k);	if ($j2 ne "") {$j2	= quotemeta($j2)};
	opendir(FDIR, $j) || die("Can not open dir: $j\n");
	foreach $file (sort readdir(FDIR)) {
#	}
#	while ($file=readdir(FDIR)) {
		if ($file=~/$j2(.*)$k$/) {
			if ($1 ne "") {
				$num{$i}->[0]	= $1;
	#			print STDERR "$i=$1, num=";
			} else {
				$num{$i}->[0]	= $i;
	#			print STDERR "$i=$1";
			}
			$num{$i}->[1]	= 0;
			$bfile	= $j ."/$file";
			print STDERR "\n\tthe $i file lib $num{$i}->[0]:\t$bfile\n";
			
#--------------------------------------------------------------------------------------------------
$k1	= 0;	$b	= "";	$a	= "";	my %nam;
$file = new FileHandle ("$bfile") || die("Cannot open miRNA file: $bfile\n");
while(<$file>)
{
#	$_=~s/^[\s|\t]+//g;
#	$_=~s/[\s|\t|\r|\n]+$//g;
	if ($_ =~/^(\S+)\s+/) {
		$a	= $1;
		if (!exists($nam{$a})) {
#			$a	=~/^(\S+)\_(\S+)\_(\d+)x?/;
#			$num{$i}->[1]	+= $3;
#			$k1++;
			if ($a	=~/^(\S+)\_(\S+)\_(\d+)x/) {
				$num{$i}->[1]	+= $3;
				$k1++;
			} else {
				$num{$i}->[1]	+= 1;
				$k1++;
			}
		}
	} else {
		print STDERR "wrong format $_\n";
	}
}
undef(%nam);
close $file || die;
print "Soap $1,\tnum=",$num{$i}->[1];	print STDERR "\tSoap $i, $1,\tnum=",$num{$i}->[1];
$num{$i}->[1]	= 1.0*$num{$i}->[1]/$normal_base;
print "\t normalize=",$num{$i}->[1],"\n";	print STDERR "\t normalize=",$num{$i}->[1];

print	STDERR "\n\tLoad file: $bfile OK\t$k1\n";
	$k3+=$k1;
#--------------------------------------------------------------------------------------------------
			$i++;	#print STDERR "i=$i\n"
		}

	}
	closedir(FDIR);
print STDERR "\nNow = [",sub_format_datetime(localtime(time())),"]\tLoad all soap result files OK. $i,$k1\n\n";
}

############################################ read library ##########################################
$i	= 0;	$k3	= 0;	$m	= 0;	my $j2 = "";	my $j3	= "";
if ($library ne "") {
	if ($library =~/^(.+\/)([^\/]*)\*([^\*]+)$/) {								# data/a*b data/*b
		$k	= $3;	$j	= $1;	$j2	= $2;print STDERR "1.j2=$j2,k=$k,j=$j\n";	# $k = b;	$j	= data/
	} elsif ($soap_result =~/^([^\/]*)\*([^\*]+)$/) {								# a*b, *b
		$k	= $2;	$j	= "./";	$j2	= $1;print STDERR "2.j2=$j2,k=$k,j=$j\n";	# $k = b;	$j	= ./
	} elsif ($soap_result =~/^(.+\/)([^\/]+)$/) {									# data/a
		$k	= $2;	$j	= $1;	$j2 = "";print STDERR "3.j2=$j2,k=$k,j=$j\n";	# $k = a;	$j	= data/
	} elsif ($soap_result =~/^([^\*]+)$/) {											# a
		$k	= $1;	$j	= "./";	$j2 = "";print STDERR "4.j2=$j2,k=$k,j=$j\n";	# $k = a;	$j	= ./
	} else {
		print STDERR	$library;
		die("reinput $library!\n");
	}
	$k	= quotemeta($k);	if ($j2 ne "") {$j2	= quotemeta($j2)};
	opendir(FDIR, $j) || die("Can not open dir: $j\n");
	foreach $file (sort readdir(FDIR)) {
#	}
#	while ($file=readdir(FDIR)) {
		if ($file=~/$j2(.*)$k$/) {
			if ($1 ne "") {
				$num{$i}->[0]	= $1;
	#			print STDERR "$i=$1, num=";
			} else {
				$num{$i}->[0]	= $i;
	#			print STDERR "$i=$1";
			}
			$num{$i}->[1]	= 0;
			$bfile	= $j ."/$file";
			print STDERR "\n\tthe $i file lib $num{$i}->[0]:\t$bfile\n";
			
#--------------------------------------------------------------------------------------------------
$k1	= 0;	$b	= "";	$a	= "";	my %nam;
$file = new FileHandle ("$bfile") || die("Cannot open miRNA file: $bfile\n");
while(<$file>)
{
#	$_=~s/^[\s|\t]+//g;
#	$_=~s/[\s|\t|\r|\n]+$//g;
	if ($_ =~/^>(\S+)\s+/) {
		$a	= $1;
		if (!exists($nam{$a})) {
			if ($a	=~/^(\S+)\_(\S+)\_(\d+)x/) {
				$num{$i}->[1]	+= $3;
				$k1++;
			} else {
				$num{$i}->[1]	+= 1;
				$k1++;
			}
		}
	} else {
#		print STDERR "wrong format $_\n";
	}
}
undef(%nam);
close $file || die;
print "Soap $1,\tnum=",$num{$i}->[1];	print STDERR "\tSoap $i, $1,\tnum=",$num{$i}->[1];
$num{$i}->[1]	= 1.0*$num{$i}->[1]/$normal_base;
print "\t normalize=",$num{$i}->[1],"\n";	print STDERR "\t normalize=",$num{$i}->[1];

print	STDERR "\n\tLoad file: $bfile OK\t$k1\n";
	$k3+=$k1;
#--------------------------------------------------------------------------------------------------
			$i++;	#print STDERR "i=$i\n"
		}

	}
	closedir(FDIR);
print STDERR "\nNow = [",sub_format_datetime(localtime(time())),"]\tLoad all library OK. $i,$k1\n\n";
}


#sub_end_program();
############################################ read file ######################################################
$i	= 0;	$k3	= 0;	$m	= 0;	my	$ss1;	
if ($libFile ne "") {
	if ($libFile =~/^(.+\/)([^\/]*)\*([^\*]+)$/) {								# data/a*b data/*b
		$k	= $3;	$j	= $1;	$j2	= $2;print STDERR "1.j2=$j2,k=$k,j=$j\n";	# $k = b;	$j	= data/
	} elsif ($libFile =~/^([^\/]*)\*([^\*]+)$/) {								# a*b, *b
		$k	= $2;	$j	= "./";	$j2	= $1;print STDERR "2.j2=$j2,k=$k,j=$j\n";	# $k = b;	$j	= ./
	} elsif ($libFile =~/^(.+\/)([^\/]+)$/) {									# data/a
		$k	= $2;	$j	= $1;	$j2 = "";print STDERR "3.j2=$j2,k=$k,j=$j\n";	# $k = a;	$j	= data/
	} elsif ($libFile =~/^([^\*]+)$/) {											# a
		$k	= $1;	$j	= "./";	$j2 = "";print STDERR "4.j2=$j2,k=$k,j=$j\n";	# $k = a;	$j	= ./
	} else {
		print STDERR	$libFile;
		die("reinput $libFile!\n");
	}
	$k	= quotemeta($k);	if ($j2 ne "") {$j2	= quotemeta($j2)};
	opendir(FDIR, $j) || die("Can not open dir: $j\n");
	foreach $file (sort readdir(FDIR)) {
#	}
#	while ($file=readdir(FDIR)) {
		if ($file=~/$j2(.*)$k$/) {
			if ($1 ne "") {
				$seg{$i}	= $1; # number ==> name
			} else {
				$seg{$i}	= $i; # number ==> number
			}
			$bfile	= $j ."/$file";
			print STDERR "\n\tthe $i file lib $seg{$i}:\t$bfile\n";
			

			foreach $key (keys %gnm) {
				if (!exists($gnm{$key}{$i})) {
					if ($inclusion == 2) {
						$gnm{$key}{$i}{0}	= 0; # match
						$gnm{$key}{$i}{1}	= 0; # truncation
						$gnm{$key}{$i}{2}	= 0; # tailing
						$gnm{$key}{$i}{3}	= 0; # truncation and tailing
					} elsif ($inclusion > 2) {
						for ($k4 = 0; $k4 <= $start_len ;$k4++) {
							for ($n = 0; $n <= $start_len ;$n++) {
								$gnm{$key}{$i}{$k4."_".$n}	= 0;
							}
						}
					} else {
						$gnm{$key}{$i}	= 0;
						$gnm{$key}{",$i"}	= 0;
					}

				} else {
					print STDERR "lib name is wrong:$key,$i!\n";
				}
			}
#--------------------------------------------------------------------------------------------------
$k1	= 0;	$b	= "";	$a	= "";	my $c	= 0;
$file = new FileHandle ("$bfile") || die("Cannot open miRNA file: $bfile\n");
while(<$file>)
{
	$_=~s/^[\s|\t]+//g;
	$_=~s/[\s|\t|\r|\n]+$//g;
	if ($_ =~/^>(\S+)\_(\S+)\_(\d+)x/) {
		$a	= $c;	$c	= $3;
		if ($b	ne "") {
			my $b1	= reverseDNAString($b);
			if ($inclusion	== 0) {
				if (exists($matureSeq{$b})) {
					foreach $key (keys %{$matureSeq{$b}}) {
						$gnm{$key}{$i}	+= $a;
					}
				}
				if (exists($matureSeq{$b1})) {
					foreach $key (keys %{$matureSeq{$b1}}) {
						$gnm{$key}{",$i"}	+= $a;
					}
				}
#				foreach $key (keys %gnm) {
#		#			if ($b=~/$gnm{$key}{"Seq"}/i || $gnm{$key}{"Seq"} =~ /$b/i) {
#					if ($b eq $gnm{$key}{"Seq"}) {
#						$gnm{$key}{$i}	+= $a;
#		#			} elsif ($b1=~/$gnm{$key}{"Seq"}/i || $gnm{$key}{"Seq"} =~ /$b1/i) {
#					} elsif ($b1 eq $gnm{$key}{"Seq"}) {
#						$gnm{$key}{",$i"}	+= $a;
#					} 
#				}
			} elsif ($inclusion	== 1) {
				foreach $key (keys %gnm) {
					if ($b=~/$gnm{$key}{"Seq"}/i || $gnm{$key}{"Seq"} =~ /$b/i) {
						$gnm{$key}{$i}	+= $a;
					} elsif ($b1=~/$gnm{$key}{"Seq"}/i || $gnm{$key}{"Seq"} =~ /$b1/i) {
						$gnm{$key}{",$i"}	+= $a;
					} 
				}
			} elsif ($inclusion	== 2) {
				my %tmpkey=();
				for ($k4 = 0; $k4 <=$start_len ;$k4++) {
					$n = substr($b,0,length($b)-$k4);
					if (exists($matureSeq{$n})) {
						foreach $key (keys %{$matureSeq{$n}}) {
							if (exists($tmpkey{$key})) {
								next;
							} else {
								$tmpkey{$key} = 1;
							}
				#			print STDERR "$b\t$a\n";
							if ($matureSeq{$n}{$key}==0 && $k4 ==0 ) {
								$gnm{$key}{$i}{0}	+= $a;
							} elsif ($matureSeq{$n}{$key}>0 && $k4 ==0 ) {
								$gnm{$key}{$i}{1}	+= $a;
							} elsif ($matureSeq{$n}{$key}==0 && $k4 > 0 ) {
								$gnm{$key}{$i}{2}	+= $a;
							} elsif ($matureSeq{$n}{$key}>0 && $k4 >0 ) {
								$gnm{$key}{$i}{3}	+= $a;
							}
						}
					}
				}
			} else {
				my %tmpkey=();
				for ($k4 = 0; $k4 <=$start_len ;$k4++) {
					$n = substr($b,0,length($b)-$k4);
					if (exists($matureSeq{$n})) {
						foreach $key (keys %{$matureSeq{$n}}) {
							if (exists($tmpkey{$key})) {
								next;
							} else {
								$tmpkey{$key} = 1;
							}
							$gnm{$key}{$i}{$matureSeq{$n}{$key}."_".$k4}	+= $a; # truncation _ tail
						}
					}
				}
			}
		}
		$b	= "";
		$k1++;	#	$k1%10000 != 0 || print STDERR "$k1,";
	} elsif ($_=~/^>/) {
		$a	= $c;	$c	= 1;
		if ($b	ne "") {
			my $b1	= reverseDNAString($b);
			if ($inclusion	== 0) {
				if (exists($matureSeq{$b})) {
					foreach $key (keys %{$matureSeq{$b}}) {
						$gnm{$key}{$i}	+= $a;
					}
				}
				if (exists($matureSeq{$b1})) {
					foreach $key (keys %{$matureSeq{$b1}}) {
						$gnm{$key}{",$i"}	+= $a;
					}
				}
			} elsif ($inclusion	== 1) {
				foreach $key (keys %gnm) {
					if ($b=~/$gnm{$key}{"Seq"}/i || $gnm{$key}{"Seq"} =~ /$b/i) {
						$gnm{$key}{$i}	+= $a;
					} elsif ($b1=~/$gnm{$key}{"Seq"}/i || $gnm{$key}{"Seq"} =~ /$b1/i) {
						$gnm{$key}{",$i"}	+= $a;
					} 
				}
			} elsif ($inclusion	== 2) {
				my %tmpkey=();
				for ($k4 = 0; $k4 <=$start_len ;$k4++) {
					$n = substr($b,0,length($b)-$k4);
					if (exists($matureSeq{$n})) {
						foreach $key (keys %{$matureSeq{$n}}) {
							if (exists($tmpkey{$key})) {
								next;
							} else {
								$tmpkey{$key} = 1;
							}
							if ($matureSeq{$n}{$key}==0 && $k4 ==0 ) {
								$gnm{$key}{$i}{0}	+= $a;
							} elsif ($matureSeq{$n}{$key}>0 && $k4 ==0 ) {
								$gnm{$key}{$i}{1}	+= $a;
							} elsif ($matureSeq{$n}{$key}==0 && $k4 > 0 ) {
								$gnm{$key}{$i}{2}	+= $a;
							} elsif ($matureSeq{$n}{$key}>0 && $k4 >0 ) {
								$gnm{$key}{$i}{3}	+= $a;
							}
						}
					}
				}
			} else {
				my %tmpkey=();
				for ($k4 = 0; $k4 <=$start_len ;$k4++) {
					$n = substr($b,0,length($b)-$k4);
					if (exists($matureSeq{$n})) {
						foreach $key (keys %{$matureSeq{$n}}) {
							if (exists($tmpkey{$key})) {
								next;
							} else {
								$tmpkey{$key} = 1;
							}
							$gnm{$key}{$i}{$matureSeq{$n}{$key}."_".$k4}	+= $a; # truncation _ tail
						}
					}
				}
			}
		}
		$b	= "";
		$k1++;	#	$k1%10000 != 0 || print STDERR "$k1,";
	}else {
		$b	.= uc($_);
	}
}
$a	= $c;
if ($b	ne "") {
	my $b1	= reverseDNAString($b);
#	foreach $key (keys %gnm) {
#	#	if ($b=~/$gnm{$key}{"Seq"}/i || $gnm{$key}{"Seq"} =~ /$b/i) {
#		if ($b eq $gnm{$key}{"Seq"}) {
#			$gnm{$key}{$i}	+= $a;
#	#	} elsif ($b1=~/$gnm{$key}{"Seq"}/i || $gnm{$key}{"Seq"} =~ /$b1/i) {
#		} elsif ($b1 eq $gnm{$key}{"Seq"}) {
#			$gnm{$key}{",$i"}	+= $a;
#		} 
#	}
	if ($inclusion == 0) {
		if (exists($matureSeq{$b})) {
			foreach $key (keys %{$matureSeq{$b}}) {
				$gnm{$key}{$i}	+= $a;
			}
		}
		if (exists($matureSeq{$b1})) {
			foreach $key (keys %{$matureSeq{$b1}}) {
				$gnm{$key}{",$i"}	+= $a;
			}
		}
#		foreach $key (keys %gnm) {
#		#	if ($b=~/$gnm{$key}{"Seq"}/i || $gnm{$key}{"Seq"} =~ /$b/i) {
#			if ($b eq $gnm{$key}{"Seq"}) {
#				$gnm{$key}{$i}	+= $a;
#		#	} elsif ($b1=~/$gnm{$key}{"Seq"}/i || $gnm{$key}{"Seq"} =~ /$b1/i) {
#			} elsif ($b1 eq $gnm{$key}{"Seq"}) {
#				$gnm{$key}{",$i"}	+= $a;
#			} 
#		}
	} elsif ($inclusion	== 1) {
		foreach $key (keys %gnm) {
			if ($b=~/$gnm{$key}{"Seq"}/i || $gnm{$key}{"Seq"} =~ /$b/i) {
				$gnm{$key}{$i}	+= $a;
			} elsif ($b1=~/$gnm{$key}{"Seq"}/i || $gnm{$key}{"Seq"} =~ /$b1/i) {
				$gnm{$key}{",$i"}	+= $a;
			} 
		}
	} elsif ($inclusion	== 2) {
		my %tmpkey=();
		for ($k4 = 0; $k4 <=$start_len ;$k4++) {
			$n = substr($b,0,length($b)-$k4);
			if (exists($matureSeq{$n})) {
				foreach $key (keys %{$matureSeq{$n}}) {
					if (exists($tmpkey{$key})) {
						next;
					} else {
						$tmpkey{$key} = 1;
					}
					if ($matureSeq{$n}{$key}==0 && $k4 ==0 ) {
						$gnm{$key}{$i}{0}	+= $a;
					} elsif ($matureSeq{$n}{$key}>0 && $k4 ==0 ) {
						$gnm{$key}{$i}{1}	+= $a;
					} elsif ($matureSeq{$n}{$key}==0 && $k4 > 0 ) {
						$gnm{$key}{$i}{2}	+= $a;
					} elsif ($matureSeq{$n}{$key}>0 && $k4 >0 ) {
						$gnm{$key}{$i}{3}	+= $a;
					}
				}
			}
		}
	} else {
		my %tmpkey=();
		for ($k4 = 0; $k4 <=$start_len ;$k4++) {
			$n = substr($b,0,length($b)-$k4);
			if (exists($matureSeq{$n})) {
				foreach $key (keys %{$matureSeq{$n}}) {
					if (exists($tmpkey{$key})) {
						next;
					} else {
						$tmpkey{$key} = 1;
					}
					$gnm{$key}{$i}{$matureSeq{$n}{$key}."_".$k4}	+= $a; # truncation _ tail
				}
			}
		}
	}
}
close $file || die;

print	STDERR "\tLoad file: $bfile OK\t$k1\n";
	$k3+=$k1;
#--------------------------------------------------------------------------------------------------
			$i++;
		}

	}
	closedir(FDIR);
print STDERR "\nNow = [",sub_format_datetime(localtime(time())),"]\tLoad all lib files: $i,$k1\n";
}

############################################ output ######################################################
print "Name\tLength";	my %libst=();	$i	= 0;
foreach $key (sort numerically keys %seg) {
	if (exists($libst{$seg{$key}})) {
		$i	= 1;
	} else {
		$libst{$seg{$key}}	= $key;# name ==> number
	}
}
$i = 1;
if ($i	== 0) {
	print STDERR "\nLib:";
	foreach $key (sort keys %libst) {
		if ($opformat == 0) {
			print "\tRaw_$key\tRaw_$key*";
#			print STDERR "\tRaw_$key\tRaw*";
		} elsif ($opformat == 1) {
			print "\tRaw_",$key;
#			print STDERR "\tRaw_",$key;
		} elsif ($opformat == 2) {
			if ($inclusion == 2) {
				print "\t$key\t$key\_Tr\t$key\_Ta\t$key\_T&T";
			} else {
				print "\t$key+";
			}
#			print STDERR "\tRaw_$key+";
		} elsif ($opformat == 3) {
			print "\tRaw_$key-";
#			print STDERR "\tRaw_$key-";
		}
	}
	foreach $key (sort keys %libst) {
		if ($opformat == 0) {
			print "\t",$key,"\t$key*";
			print STDERR "\t",$key,"\t*";
		} elsif ($opformat == 1) {
			print "\t",$key;
			print STDERR "\t",$key;
		} elsif ($opformat == 2) {
			if ($inclusion == 2) {
				print "\t$key\t$key\_Tr\t$key\_Ta\t$key\_T&T";
			} else {
				print "\t$key+";
			}
			print STDERR "\t$key+";
		} elsif ($opformat == 3) {
			print "\t$key-";
			print STDERR "\t$key-";
		}
	}
	print "\n";	print STDERR "\n";
	if ($soap_result ne "" || $library ne "" || $lib_size ne "") {
		$j2	= 0;
		print STDERR "\nSoap:";
		foreach $key (sort numerically keys %num) {
			print STDERR "\t$num{$key}->[0]\t$key";
			if ($normal_base > 1) {
				$num{$key}->[1]	= 1.0*$num{$key}->[1];#/$normal_base;
			}
			$buf[$j2]	= $j2;#	print STDERR "buf[$j2]=$buf[$j2]\t";
			$j2++;
		}
		if ($ask_input == 1) {
			print STDERR "\n\n1.please input the num in soap result: e.g. 3,2,0,1; \[0,1,2,3\]\n";
			$ss1	= <STDIN>; $ss1	=~s/[\s|\t|\r|\n]+$//g;
			if ($ss1 ne "") {
				@buf=split(/\,/,$ss1);
			}
		}
		foreach $key (sort keys %gnm) {
			print $key,"\t",length($gnm{$key}{"Seq"});	$j2	= 0;
			foreach $k1 (sort keys %libst) {
				if ($opformat == 0) {
					print "\t",$gnm{$key}{$libst{$k1}},"\t",$gnm{$key}{",$libst{$k1}"};
				} elsif ($opformat == 1) {
					print "\t",$gnm{$key}{$libst{$k1}}+$gnm{$key}{",$libst{$k1}"};
				} elsif ($opformat == 2) {
					if ($inclusion < 2) {
						print "\t",$gnm{$key}{$libst{$k1}};
					} elsif ($inclusion == 2) {
						print "\t",$gnm{$key}{$libst{$k1}}{0},"\t",$gnm{$key}{$libst{$k1}}{1},"\t",
							$gnm{$key}{$libst{$k1}}{2},"\t",$gnm{$key}{$libst{$k1}}{3};
					} else {
						print "\t";
						for ($k4 = 0; $k4 <= $start_len ;$k4++) {
							for ($n = 0; $n <= $start_len ;$n++) {
								print ";",$gnm{$key}{$libst{$k1}}{$k4."_".$n};
							}
						}
					}
				} elsif ($opformat == 3) {
					print "\t",$gnm{$key}{",$libst{$k1}"};
				}
			}
			foreach $k1 (sort keys %libst) {
				if ($opformat == 0) {
					printf("\t%.1f\t%.1f",$gnm{$key}{$libst{$k1}}/$num{$buf[$j2]}->[1],$gnm{$key}{",$libst{$k1}"}/$num{$buf[$j2]}->[1]);
				} elsif ($opformat == 1) {
					printf("\t%.1f",$gnm{$key}{$libst{$k1}}/$num{$buf[$j2]}->[1]+$gnm{$key}{",$libst{$k1}"}/$num{$buf[$j2]}->[1]);
				} elsif ($opformat == 2) {
					#	printf("\t%.1f",$gnm{$key}{$libst{$k1}}/$num{$buf[$j2]}->[1]);
					if ($inclusion < 2) {
					#	print "\t",$gnm{$key}{$libst{$k1}};
						printf("\t%.1f",$gnm{$key}{$libst{$k1}}/$num{$buf[$j2]}->[1]);
					} elsif ($inclusion == 2) {
						printf("\t%.1f\t%.1f\t%.1f\t%.1f",$gnm{$key}{$libst{$k1}}{0}/$num{$buf[$j2]}->[1],$gnm{$key}{$libst{$k1}}{1}/$num{$buf[$j2]}->[1],
							$gnm{$key}{$libst{$k1}}{2}/$num{$buf[$j2]}->[1],$gnm{$key}{$libst{$k1}}{3}/$num{$buf[$j2]}->[1]);
			#			print "\t",$gnm{$key}{$libst{$k1}}{0}"\t",$gnm{$key}{$libst{$k1}}{1}"\t",$gnm{$key}{$libst{$k1}}{2}"\t",$gnm{$key}{$libst{$k1}}{3};
					} else {
						$m = "\t";
						for ($k4 = 0; $k4 <= $start_len ;$k4++) {
							for ($n = 0; $n <= $start_len ;$n++) {
								$m.= ";",$gnm{$key}{$libst{$k1}}{$k4."_".$n}/$num{$buf[$j2]}->[1];
							}
						}
						print $m;
					}
				} elsif ($opformat == 3) {
					printf("\t%.1f",$gnm{$key}{",$libst{$k1}"}/$num{$buf[$j2]}->[1]);
				}
				$j2++;
			}
			print "\n";
		}
	} else {
		foreach $key (sort keys %gnm) {
			print $key,"\t",length($gnm{$key}{"Seq"});
			foreach $k1 (sort keys %libst) {
				if ($opformat == 0) {
					print "\t",$gnm{$key}{$libst{$k1}},"\t",$gnm{$key}{",$libst{$k1}"};
				} elsif ($opformat == 1) {
					print "\t",$gnm{$key}{$libst{$k1}}+$gnm{$key}{",$libst{$k1}"};
				} elsif ($opformat == 2) {
				#	print "\t",$gnm{$key}{$libst{$k1}};
					if ($inclusion < 2) {
						print "\t",$gnm{$key}{$libst{$k1}};
					} elsif ($inclusion == 2) {
						print "\t",$gnm{$key}{$libst{$k1}}{0},"\t",$gnm{$key}{$libst{$k1}}{1},"\t",
							$gnm{$key}{$libst{$k1}}{2},"\t",$gnm{$key}{$libst{$k1}}{3};
					} else {
						print "\t";
						for ($k4 = 0; $k4 <= $start_len ;$k4++) {
							for ($n = 0; $n <= $start_len ;$n++) {
								print ";",$gnm{$key}{$libst{$k1}}{$k4."_".$n};
							}
						}
					}
				} elsif ($opformat == 3) {
					print "\t",$gnm{$key}{",$libst{$k1}"};
				}
			}
			print "\n";
		}
	}
} else {
	foreach $key (sort numerically keys %seg) {
		if ($opformat == 0) {
			print "\tRaw_$seg{$key}\tRaw*";
#			print STDERR "\tRaw_$seg{$key}\tRaw*";
		} elsif ($opformat == 1) {
			print "\tRaw_",$seg{$key};
#			print STDERR "\tRaw_",$seg{$key};
		} elsif ($opformat == 2) {
		#	print "\tRaw_$seg{$key}+";
			if ($inclusion == 2) {
				print "\t$seg{$key}\t$seg{$key}\_Tr\t$seg{$key}\_Ta\t$seg{$key}\_T&T";
			} else {
				print "\t$seg{$key}";
			}
#			print STDERR "\tRaw_$seg{$key}+";
		} elsif ($opformat == 3) {
			print "\tRaw_$seg{$key}-";
#			print STDERR "\tRaw_$seg{$key}-";
		}
	}
	foreach $key (sort numerically keys %seg) {
		if ($opformat == 0) {
			print "\t",$seg{$key},"\t*";
			print STDERR "\t",$seg{$key},"\t*";
		} elsif ($opformat == 1) {
			print "\t",$seg{$key};
			print STDERR "\t",$seg{$key};
		} elsif ($opformat == 2) {
			if ($inclusion == 2) {
				print "\t$seg{$key}\t$seg{$key}\_Tr\t$seg{$key}\_Ta\t$seg{$key}\_T&T";
			} else {
				print "\t$seg{$key}";
			}
		#	print "\t$seg{$key}+";
		#	print STDERR "\t$seg{$key}+";
		} elsif ($opformat == 3) {
			print "\t$seg{$key}-";
			print STDERR "\t$seg{$key}-";
		}
	}
	print "\n";	print STDERR "\n";
	if ($soap_result ne "" || $library ne "" || $lib_size ne "") {
		$j2	= 0;
#		for (;$j2 < 4 ;$j2++) {
#			$buf[$j2]	= $j2;
#		}
		foreach $key (sort numerically keys %num) {
			print STDERR "\t$key\t",$num{$key}->[0];
			if ($normal_base > 1) {
				$num{$key}->[1]	= 1.0*$num{$key}->[1];#/$normal_base;
			}
			$buf[$j2]	= $j2;	#print STDERR "buf[$j2]=$buf[$j2]\t";
			$j2++;
		}
		if ($ask_input == 1) {
			print STDERR "\n2.please input the num in soap result: e.g. 3,2,0,1; [0,1,2,3]\n";
			$ss1	= <STDIN>; $ss1	=~s/[\s|\t|\r|\n]+$//g;
			if ($ss1 ne "") {
	#			@buf=split(/\,/,$ss1);
			}
		}
		
		foreach $key (sort keys %gnm) {
			print $key,"\t",length($gnm{$key}{"Seq"});	$j2	= 0;
			foreach $k1 (sort numerically keys %seg) {
				if ($opformat == 0) {
					print "\t",$gnm{$key}{$k1},"\t",$gnm{$key}{",$k1"};
				} elsif ($opformat == 1) {
					print "\t",$gnm{$key}{$k1}+$gnm{$key}{",$k1"};
				} elsif ($opformat == 2) {
					#print "\t",$gnm{$key}{$k1};
					if ($inclusion < 2) {
						print "\t",$gnm{$key}{$k1};
					} elsif ($inclusion == 2) {
						print "\t",$gnm{$key}{$k1}{0},"\t",$gnm{$key}{$k1}{1},"\t",$gnm{$key}{$k1}{2},"\t",$gnm{$key}{$k1}{3};
					} else {
						print "\t";
						for ($k4 = 0; $k4 <= $start_len ;$k4++) {
							for ($n = 0; $n <= $start_len ;$n++) {
								print ";",$gnm{$key}{$k1}{$k4."_".$n};
							}
						}
					}
				} elsif ($opformat == 3) {
					print "\t",$gnm{$key}{",$k1"};
				}
			}
			foreach $k1 (sort numerically keys %seg) {
				if ($opformat == 0) {
					printf("\t%.1f\t%.1f",$gnm{$key}{$k1}/$num{$buf[$j2]}->[1],$gnm{$key}{",$k1"}/$num{$buf[$j2]}->[1]);
				} elsif ($opformat == 1) {
					printf("\t%.1f",$gnm{$key}{$k1}/$num{$buf[$j2]}->[1]+$gnm{$key}{",$k1"}/$num{$buf[$j2]}->[1]);
				} elsif ($opformat == 2) {
					#printf("\t%.1f",$gnm{$key}{$k1}/$num{$buf[$j2]}->[1]);
					if ($inclusion < 2) {
					#	print "\t",$gnm{$key}{$k1};
						printf("\t%.1f",$gnm{$key}{$k1}/$num{$buf[$j2]}->[1]);
					} elsif ($inclusion == 2) {
						printf("\t%.1f\t%.1f\t%.1f\t%.1f",$gnm{$key}{$k1}{0}/$num{$buf[$j2]}->[1],$gnm{$key}{$k1}{1}/$num{$buf[$j2]}->[1],
							$gnm{$key}{$k1}{2}/$num{$buf[$j2]}->[1],$gnm{$key}{$k1}{3}/$num{$buf[$j2]}->[1]);
			#			print "\t",$gnm{$key}{$k1}{0}"\t",$gnm{$key}{$k1}{1}"\t",$gnm{$key}{$k1}{2}"\t",$gnm{$key}{$k1}{3};
					} else {
						$m = "\t";
						for ($k4 = 0; $k4 <= $start_len ;$k4++) {
							for ($n = 0; $n <= $start_len ;$n++) {
								$m.= ";",$gnm{$key}{$k1}{$k4."_".$n}/$num{$buf[$j2]}->[1];
							}
						}
						print $m;
					}
				} elsif ($opformat == 3) {
					printf("\t%.1f",$gnm{$key}{",$k1"}/$num{$buf[$j2]}->[1]);
				}
				$j2++;
			}
			print "\n";
		}
	} else {
		foreach $key (sort keys %gnm) {
			print $key,"\t",length($gnm{$key}{"Seq"});
			foreach $k1 (sort numerically keys %seg) {
				if ($opformat == 0) {
					print "\t",$gnm{$key}{$k1},"\t",$gnm{$key}{",$k1"};
				} elsif ($opformat == 1) {
					print "\t",$gnm{$key}{$k1}+$gnm{$key}{",$k1"};
				} elsif ($opformat == 2) {
				#	print "\t",$gnm{$key}{$k1};
					if ($inclusion < 2) {
						print "\t",$gnm{$key}{$k1};
					} elsif ($inclusion == 2) {
						print "\t",$gnm{$key}{$k1}{0},"\t",$gnm{$key}{$k1}{1},"\t",$gnm{$key}{$k1}{2},"\t",$gnm{$key}{$k1}{3};
					} else {
						print "\t";
						for ($k4 = 0; $k4 <= $start_len ;$k4++) {
							for ($n = 0; $n <= $start_len ;$n++) {
								print ";",$gnm{$key}{$k1}{$k4."_".$n};
							}
						}
					}
				} elsif ($opformat == 3) {
					print "\t",$gnm{$key}{",$k1"};
				}
			}
			print "\n";
		}
	}
}
sub_end_program();


#############################################################################################################
####################################                                         ################################
####################################              "main end"                 ################################
####################################                                         ################################
#############################################################################################################


sub usage
{
	print "Program :\t$0\n";
	print "Version :\t$version\n";
	print "Author  :\tLei Gao, UC,Riverside\n";
	print "Contact :\tLei Gao <highlei\@gmail.com>\n";
	print "\nUsage:	$0 [options]\n";
	print "\t-i	<str>	input the miRNA file.";
	print " eg: maize10.2cn20.-CDS.1.all.mature.fa\n";
	print "\t-d	<str>	input the library file(s).";
	print " eg: \"data/zea.*fa\"\n";
	print "\t-S	<str>	input soap result file(s).";
	print " eg: \"data/zea.*soap\"\n";
	print "\t-N	<int>	input the base for normalization.";
	print " [$normal_base]\n";
	print "\t-s	<int>	the library size.";
	print " [$lib_size]\n";
	print "\t-L	<str>	input library.";
	print " [$library]\n";
	print "\t-O	<int>	output format (0: +,-; 1:+ plus -; 2:+; 3:-).";
	print " [$opformat]\n";
	print "\t-I	<int>	0: match; 1: include; 2: 4kind Truncation&Tailing; 3: all T&T.";
	print " [$inclusion]\n";
	print "\t-A	<int>	ask iput (1) or not (0).";
	print " [$ask_input]\n";
	print "\t-l	<int>	the max truncation length.";
	print " [$start_len]\n";
#	print "\t-c	<int>	the score is allowed.";
#	print " [$score]\n";
	print "\n\t-h	display this help\n";
#	print "		Note: please add quotation mark, if you input parameter in command line!\n";
	print "\nExample:\n";
	print "$0 -i maize10.2cn20.-CDS.1.all.mature.fa -d \"data/zea.*fa\"\n";
	print ("==========================| $0  end  |==================================\n\n");

    exit(0);
}
############################################################################################################
######################                  sub_format_datetime
############################################################################################################

sub sub_format_datetime #时间子程序
{
    my ($sec, $min, $hour, $day, $mon, $year, $wday, $yday, $isdst) = @_;
    sprintf("%4d-%02d-%02d %02d:%02d:%02d", $year+1900, $mon+1, $day, $hour, $min, $sec);
}

############################################################################################################
######################                  sub_end_program
############################################################################################################
sub sub_end_program
{
	print STDERR ("\n............................................................\n");
	my $Time_End = sub_format_datetime(localtime(time()));
	print STDERR "Running from [$Time_Start] to [$Time_End]\n";
	$end = time();
	printf STDERR ("Total execute time : %.2f s\n",$end-$start);
	print STDERR ("==========================| $0  end  |==================================\n\n");
	exit(0);

}
############################################################################################################
######################                  reverseDNAString
############################################################################################################

sub reverseDNAString 
{
	my($rdstr)	= @_;
	my ($sr1,$sr2);
	$rdstr	= reverse($rdstr);
	$rdstr	=~tr/ACGTRYMKacgtrymk/TGCAYRKMtgcayrkm/;
	return $rdstr;
}
