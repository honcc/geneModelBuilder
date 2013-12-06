#!/usr/bin/perl/ -w

use strict;
#use Math::Random; #---For generating normal distribution random numbers
use File::Path; #---for removing tmp GNUPLOT dat files

######################################################################################################################################################
#
#	Description
#		This is a perl script to build the gene models by combining the transfrags and the splicing information; its would also look for alternative 
#	splicing;
#
#	Input
#		--NGSGffPath=			the transfrag GFF; normally generated from "transfragDiscoverer.pl";
#		--junctionBedPath=		the bed of the junction info; normally generated from HMMSplicerBEDToSAMParser.pl;
#		--pileupPath=			the pileup file generated from pileup counter, used to determine the ration of splice and unspliced read flanking the splicing jucntions;
#		--minSplicingRatio=		the minimum ratio of spliced to unspliced read for a junction to be regarded as a "major" junction;
#		--minMajorSupportRead=	the minimum number of supporting read for a junction to be regarded as a "major" junction;
#		--intronBoundCovPath=	"no" or the path of a file that contains the intron bound coverage of the junctions. If "no", the coverage will be calculated from the pileu file;
#		--junctionBedType=		"HMMSplicer" or "tophat"; as the bed is a bit different between the two, mainly on supporting read numnber; "HMMSplicer" fits the output from HMMSplicerBEDToSAMParser, "tophat" fits for junctions.bed from tophat;
#		--outDir= 				output directory
#
#	Output
#
#
#	Usage
#		perl geneModelBuilder_v0.2.pl --NGSGffPath=/Volumes/A_MPro2TB/NGS/full/s2/pileupCounter/transfrag/G5_C2_L100_W40_F5.both.gff --junctionBedPath=/Volumes/A_MPro2TB/NGS/full/s2/HMMSplicer/junctionInfo/HMMSplicer.filter.can.unique.dup.bed ---outDir=/Volumes/A_MPro2TB/NGS/full/s2/geneModel/builder/ --pileupPath=/Volumes/A_MPro2TB/NGS/full/s2/pileupCounter/HMMSplicer.combined.sorted.sam_c.1/featureAndReadPileup.txt --intronBoundCovPath=/Volumes/A_MPro2TB/NGS/full/s2/geneModel/builder/intronBoundCov.txt --minSplicingRatio=1 --minMajorSupportRead=5
#
#
#	Version
#
#		v0.2
#		- The recursive search of fragment junction overlapping now changed to fragment junction cluster overlapping, since the former case will not pick up junction-junction overlapping as the search is based on fragment junction overlapping;
# 
#		v0.3
#		- added junctionBedType option;
######################################################################################################################################################

#==========================================================Main body starts==========================================================================#

#1----------Read the parameters----------#
use vars qw ($NGSGffPath $junctionBedPath $pileupPath $minSplicingRatio $intronBoundCovPath $minMajorSupportRead $junctionBedType $outDir);
my ($NGSGffPath, $junctionBedPath, $pileupPath, $minSplicingRatio, $intronBoundCovPath, $minMajorSupportRead, $junctionBedType, $outDir) = readParameters();
printCMDLogOrFinishMessage("CMDLog");

#2----------Read the Gff----------#
my ($trnsfgStrdHsh_ref, $trnsfgRngHsh_ref, $trnsfgCntgHsh_ref) = readGffNoAltNoIntronNoCDS($NGSGffPath);

#3----------Read the BED----------#
my ($junctStrdHsh_ref, $junctRngHsh_ref, $junctCntgHsh_ref, $junctScoreHsh_ref, $junctReadNumHsh_ref, $junctIntronRngHsh_ref) = readJunctBED($junctionBedPath);

#4---------Compute the overlapping junctions to define overlapping whole BED ranges (implemented in v0.2)
my ($SSOvrlpJunctBEDRngHsh_ref, $dummy1, $XSOvrlpJunctBEDRngHsh_ref, $dummy2) = findFturOverlap($junctRngHsh_ref, $junctRngHsh_ref, $junctStrdHsh_ref, $junctStrdHsh_ref, "no", "no", "all", "all");

#5-------define the overlapping junct BED ranges that is used for clustering with the transfrags
my ($junctBEDRngClusterJunctHsh_ref, $junctBEDRngClusterNameByJunctHsh_ref, $junctBEDClusterRngHsh_ref, $junctBEDClusterStrdHsh_ref, $junctBEDClusterCntgHsh_ref) = defineOverlappingJunctBEDRng($SSOvrlpJunctBEDRngHsh_ref, $junctRngHsh_ref, $junctCntgHsh_ref, $junctStrdHsh_ref);

#6---------Compute the overlapping transfrag and junction BED clusters
my ($plusSSHitByTrnsfgHsh_ref, $plusSSHitByJunctHsh_ref, $plusXSHitByTrnsfgHsh_ref, $plusXSHitByJunctHsh_ref) = findFturOverlap($trnsfgRngHsh_ref, $junctBEDClusterRngHsh_ref, $trnsfgStrdHsh_ref, $junctBEDClusterStrdHsh_ref, "yes", "no", "all", "+");
my ($minusSSHitByTrnsfgHsh_ref, $minusSSHitByJunctHsh_ref, $minusXSHitByTrnsfgHsh_ref, $minusXSHitByJunctHsh_ref) = findFturOverlap($trnsfgRngHsh_ref, $junctBEDClusterRngHsh_ref, $trnsfgStrdHsh_ref, $junctBEDClusterStrdHsh_ref, "yes", "no", "all", "-");

#7---------Collapse the overlapping transfrag and jucntion clusters
my $clusterNum = 0;
our %clusterJunctTrnsfgHsh = ();
our $clusterJunctTrnsfgHsh_ref = \%clusterJunctTrnsfgHsh;
($clusterJunctTrnsfgHsh_ref, $clusterNum)= collpaseOverlappingTransfragJunctionBEDRngCluster($plusXSHitByTrnsfgHsh_ref, $plusXSHitByJunctHsh_ref, $trnsfgStrdHsh_ref, $trnsfgRngHsh_ref, $trnsfgCntgHsh_ref, $junctStrdHsh_ref, $junctRngHsh_ref, $junctCntgHsh_ref, $clusterJunctTrnsfgHsh_ref, $clusterNum, "+");
%clusterJunctTrnsfgHsh = %{$clusterJunctTrnsfgHsh_ref};
$clusterJunctTrnsfgHsh_ref = \%clusterJunctTrnsfgHsh;
($clusterJunctTrnsfgHsh_ref, $clusterNum)= collpaseOverlappingTransfragJunctionBEDRngCluster($minusXSHitByTrnsfgHsh_ref, $minusXSHitByJunctHsh_ref, $trnsfgStrdHsh_ref, $trnsfgRngHsh_ref, $trnsfgCntgHsh_ref, $junctStrdHsh_ref, $junctRngHsh_ref, $junctCntgHsh_ref, $clusterJunctTrnsfgHsh_ref, $clusterNum, "-");
%clusterJunctTrnsfgHsh = %{$clusterJunctTrnsfgHsh_ref};

#8---print the junct BED cluster name in the JT cluster back to junctStr and print temporary cluster info
convertAndPrintJunctTransfragCluster(\%clusterJunctTrnsfgHsh, $junctBEDRngClusterJunctHsh_ref);

#9---------Compute the overlapping junctions (only intron but not the whole range)
my ($SSOvrlpJunctIntronHsh_ref, $dummy3, $XSOvrlpJunctHsh_ref, $dummy4) = findFturOverlap($junctIntronRngHsh_ref, $junctIntronRngHsh_ref, $junctStrdHsh_ref, $junctStrdHsh_ref, "no", "no", "all", "all");

#10---scanJunctionIntronBoundCoverage
my $junctIntronBoundCovHsh_ref = scanJunctionIntronBoundCoverage($junctIntronRngHsh_ref);

#11----------define the intron clusters
my ($junctSplicingRatioHsh_ref, $majorJunctAvgSplcingRatioHsh_ref, $superJunctHsh_ref, $clusterAllJunctHsh_ref, $junctClusterNameByJunctHsh_ref, $junctClusterInfoHsh_ref, $junctClusterSSOverlapHsh_ref, $totalExonSkippingTypeHsh_ref, $majorExonSkippingTypeHsh_ref, $exactExonSkippingClusterHsh_ref, $splicingSiteDiffAllHsh_ref, $splicingSiteDiffHshByClusterHsh_ref, $junctOnProminentIsofmHsh_ref) = defineIntronClusters($SSOvrlpJunctIntronHsh_ref, $junctStrdHsh_ref, $junctScoreHsh_ref, $junctReadNumHsh_ref, $junctIntronBoundCovHsh_ref, $junctIntronRngHsh_ref, $XSOvrlpJunctHsh_ref);

#12---------build the gene models
buildGeneModel ($clusterJunctTrnsfgHsh_ref,  $junctOnProminentIsofmHsh_ref,  $junctClusterInfoHsh_ref,  $junctClusterSSOverlapHsh_ref,  $trnsfgRngHsh_ref,  $junctRngHsh_ref, $trnsfgCntgHsh_ref, $trnsfgStrdHsh_ref, $junctClusterNameByJunctHsh_ref);

#------print hit log----------#
printCMDLogOrFinishMessage("finishMessage");

exit;
#========================================================= Main body ends ===========================================================================#

########################################################################## readParameters
sub readParameters {
	
	$junctionBedType = "tophat";

	foreach my $param (@ARGV) {

		if ($param =~ m/--NGSGffPath=/) {$NGSGffPath = substr ($param, index ($param, "=")+1);}
		elsif ($param =~ m/--junctionBedPath=/) {$junctionBedPath = substr ($param, index ($param, "=")+1);}
		elsif ($param =~ m/--pileupPath=/) {$pileupPath = substr ($param, index ($param, "=")+1);}
		elsif ($param =~ m/--minSplicingRatio=/) {$minSplicingRatio = substr ($param, index ($param, "=")+1);}
		elsif ($param =~ m/--intronBoundCovPath=/) {$intronBoundCovPath = substr ($param, index ($param, "=")+1);}
		elsif ($param =~ m/--minMajorSupportRead=/) {$minMajorSupportRead = substr ($param, index ($param, "=")+1);}
		elsif ($param =~ m/--junctionBedType=/) {$junctionBedType = substr ($param, index ($param, "=")+1);}
		elsif ($param =~ m/--outDir=/) {$outDir = substr ($param, index ($param, "=")+1);}
	}
	
	open (TEST, "$NGSGffPath") || die "Cannot open $NGSGffPath: $!\n"; close TEST;
	open (TEST, "$junctionBedPath") || die "Cannot open $junctionBedPath: $!\n"; close TEST;
	open (TEST, "$pileupPath") || die "Cannot open $pileupPath: $!\n"; close TEST;
	if ($intronBoundCovPath ne "no") {
		open (TEST, "$intronBoundCovPath") || die "Cannot open $intronBoundCovPath: $!\n"; close TEST;
	}
	
	chop $outDir if ($outDir =~ m/\/$/); #---remove the last slash
	system "mkdir -p -m 777 $outDir/";
	system "mkdir -p -m 777 $outDir/log/";
	system "mkdir -p -m 777 $outDir/plotData/";
	system "mkdir -p -m 777 $outDir/plotPdf/";
	system "mkdir -p -m 777 $outDir/cov/";
	system "mkdir -p -m 777 $outDir/finalGTF/";
	
	return ($NGSGffPath, $junctionBedPath, $pileupPath, $minSplicingRatio, $intronBoundCovPath, $minMajorSupportRead, $junctionBedType, $outDir);
}
########################################################################## findFturOverlap
sub findFturOverlap {

#					The 7 scenes of overlapping and proximity 
#
#
#     case 0: complete overlapp (($refStart == $qryStart) && ($refEnd == $qryEnd))
#			
#     case 1: overlapHead    case 2: overlapTail      case 3: cover		      case 4: within		case 5: prxmtyHead	     case 6: prxmtyTail
#
#r	     |--------|		        |---------|	        |-------------|	              |-----|			    |-----|				       	                |-------|
#q	<=========>	           <==========>		          <=========>	            <==========>		     	    <==========>	      <==========>
#
#   ($refStart<$qryStart)&&	($refStart>=$qryStart)&&  ($refStart<$qryStart)&&  ($refStart>$qryStart)&&   ($refStart<$qryStart)&&	   ($refEnd>$qryStart)&&
#   ($refEnd>=$qryStart)&&	($refStart<=$qryEnd)&&	  ($refEnd>$qryEnd)	       ($refEnd<$qryEnd)	     ($refStart<$qryEnd)		   ($refStart>$qryEnd)
#   ($refEnd<=$qryEnd)	    ($refEnd>$qryEnd)												 
#

	my %refRangeXStrdHsh = %{$_[0]};
	my %qryRangeXStrdHsh = %{$_[1]}; 
	my %refStrdByFturHsh = %{$_[2]};
	my %qryStrdByFturHsh = %{$_[3]};
	my $reportExactOverlap = $_[4]; #----yes or no, if no, exactly overlap i.e. case 0 will not be reported. This is designed to prevent self-hit if reference and query are the same hash;
	my $verbose = $_[5]; #----yes or no
	my $refStrdFilter = $_[6];#--- only the indicated strd will be process, use all to indicate all
	my $qryStrdFilter = $_[7];#--- only the indicated strd will be process, use all to indicate all
	
	my (%SSHitByRefHsh, %SSHitByQryHsh, %XStrdHitByRefHsh, %XStrdHitByQryHsh);

	foreach my $cntg (sort {$a cmp $b} keys %refRangeXStrdHsh) {
		print "Finding overlapping features on $cntg\n" if ($verbose eq "yes");
		if (exists $qryRangeXStrdHsh{$cntg}) {#---if there are ftur on the $strd of $cntg of qryGffPath
			foreach my $refFtur (sort {$a cmp $b} keys %{$refRangeXStrdHsh{$cntg}}) {#--- all ftur on the $strd of $cntg of refGff
				
				#---skip if the strd of the refFtur doesnt match the filter
				next if (($refStrdFilter ne "all") and ($refStrdByFturHsh{$refFtur} ne $refStrdFilter));
				
				my $refStart = ${${$refRangeXStrdHsh{$cntg}}{$refFtur}}{"start"};
				my $refEnd = ${${$refRangeXStrdHsh{$cntg}}{$refFtur}}{"end"};
				foreach  my $qryFtur (sort {$a cmp $b} keys %{$qryRangeXStrdHsh{$cntg}}) {#--- all ftur on the $strd of $cntg of qryGffPath
					
					#---skip if the strd of the refFtur doesnt match the filter
					next if (($qryStrdFilter ne "all") and ($qryStrdByFturHsh{$qryFtur} ne $qryStrdFilter));

					my $sameStrd = "yes";
					
					#---two futrs are not on the same strd and strd has to be + or -;  (i.e. not . or * for no strand with be hit with both + or - in the SS manner)
					if (($refStrdByFturHsh{$refFtur} ne $qryStrdByFturHsh{$qryFtur}) and (($qryStrdByFturHsh{$qryFtur} eq "+") or ($qryStrdByFturHsh{$qryFtur} eq "-")) and (($refStrdByFturHsh{$refFtur} eq "+") or ($refStrdByFturHsh{$refFtur} eq "-"))) {
						$sameStrd = "no";
					}
					
					my $qryStart = ${${$qryRangeXStrdHsh{$cntg}}{$qryFtur}}{"start"};
					my $qryEnd = ${${$qryRangeXStrdHsh{$cntg}}{$qryFtur}}{"end"};
					
					if  (($refStart == $qryStart) && ($refEnd == $qryEnd)) {#---scene 0

						if ($reportExactOverlap eq "yes") {#----to prevent selfhit when ref and qry is the same hash.
							${$XStrdHitByRefHsh{$refFtur}}{$qryFtur} = 0;
							${$XStrdHitByQryHsh{$qryFtur}}{$refFtur} = 0;

							if ($sameStrd eq "yes") {	
								${$SSHitByRefHsh{$refFtur}}{$qryFtur} = 0;
								${$SSHitByQryHsh{$qryFtur}}{$refFtur} = 0;
							}

						} else {

							if ($sameStrd eq "no") {#---no on the same strand	
								${$XStrdHitByRefHsh{$refFtur}}{$qryFtur} = 0;
								${$XStrdHitByQryHsh{$qryFtur}}{$refFtur} = 0;
							}
						}
						
					} elsif (($refStart<=$qryStart)&&($refEnd>=$qryStart)&&($refEnd<=$qryEnd)) {#---scene 1		

						${$XStrdHitByRefHsh{$refFtur}}{$qryFtur} = 1;
						${$XStrdHitByQryHsh{$qryFtur}}{$refFtur} = 1;

						if ($sameStrd eq "yes") {
							${$SSHitByRefHsh{$refFtur}}{$qryFtur} = 1;
							${$SSHitByQryHsh{$qryFtur}}{$refFtur} = 1;
						}
											
					} elsif (($refStart>=$qryStart)&&($refStart<=$qryEnd)&&($refEnd>=$qryEnd)) {#---scene 2					

						${$XStrdHitByRefHsh{$refFtur}}{$qryFtur} = 2;
						${$XStrdHitByQryHsh{$qryFtur}}{$refFtur} = 2;

						if ($sameStrd eq "yes") {
							${$SSHitByRefHsh{$refFtur}}{$qryFtur} = 2;
							${$SSHitByQryHsh{$qryFtur}}{$refFtur} = 2;
						}
					
					} elsif (($refStart<=$qryStart)&&($refEnd>=$qryEnd)) {#---scene 3		

						${$XStrdHitByRefHsh{$refFtur}}{$qryFtur} = 3;
						${$XStrdHitByQryHsh{$qryFtur}}{$refFtur} = 3;

						if ($sameStrd eq "yes") {
							${$SSHitByRefHsh{$refFtur}}{$qryFtur} = 3;
							${$SSHitByQryHsh{$qryFtur}}{$refFtur} = 3;
						}
					
					} elsif (($refStart>=$qryStart)&&($refEnd<=$qryEnd)) {#---scene 4						

						${$XStrdHitByRefHsh{$refFtur}}{$qryFtur} = 4;
						${$XStrdHitByQryHsh{$qryFtur}}{$refFtur} = 4;

						if ($sameStrd eq "yes") {
							${$SSHitByRefHsh{$refFtur}}{$qryFtur} = 4;
							${$SSHitByQryHsh{$qryFtur}}{$refFtur} = 4;
						}
					
					} elsif (($refStart<$qryStart)&&($refStart<$qryEnd)) {#---scene 5
						#---non-overlapping, nothing to be done;						
					} elsif (($refEnd>$qryStart)&&($refStart>$qryEnd)) {#---scene 6
						#---non-overlapping, nothing to be done;						
					} else {#---BUG! possibly other scene?
						die "Unexpected overlapping scene between $refFtur and $qryFtur. It's a Bug. Program qutting.\n";
					}
				}
			}
		}
	}
	
	return (\%SSHitByRefHsh, \%SSHitByQryHsh, \%XStrdHitByRefHsh, \%XStrdHitByQryHsh);
}
########################################################################## readJunctBED
sub readJunctBED {

	my $junctBEDPath = $_[0];
	
	my (%junctRngHsh, %junctStrdHsh, %junctCtngHsh, %junctReadNumHsh, %junctScoreHsh, %junctIntronRngHsh);
	
	my $dupJunctNum = 0;
	
	open (INFILE, "$junctBEDPath");
	open (TMPLOG, ">$outDir/log/junctionInBedMoreThanOnce.txt");
	print "Reading $junctBEDPath\n";
	while (my $theLine = <INFILE>) {
		chomp $theLine;
		next if ($theLine =~ m/^track name/);
		my @theLineSplt = split (/\t/, $theLine);
		my $cntg = $theLineSplt[0];
		my $bedStart = $theLineSplt[1];
		my $bedEnd = $theLineSplt[2];
		my $strd = $theLineSplt[5];
		my @blkSizesSplt = split /,/, $theLineSplt[10];
		my $blk1Size = $blkSizesSplt[0];
		my $blk2Size = $blkSizesSplt[1];
		my $readNum = my $score = 0;
		
		if ($junctionBedType eq "HMMSplicer") {
			my @readNumAndScoreSplt = split /\|/, $theLineSplt[3];
			$readNum = substr ($readNumAndScoreSplt[0], index ($readNumAndScoreSplt[0], "=")+1);
			$score = substr ($readNumAndScoreSplt[1], index ($readNumAndScoreSplt[1], "=")+1);
		} elsif ($junctionBedType eq "tophat") {
			$readNum = $score = $theLineSplt[4];
		} else {
			die "unspecified junction bed type\n";
		}
		
		my $intronStart = $bedStart + $blk1Size + 1;
		my $intronEnd = $bedEnd - $blk2Size;

		my $junctStr = $cntg.":".$intronStart.":".$intronEnd; #---assumed to be unique

		${${$junctIntronRngHsh{$cntg}}{$junctStr}}{"start"} = $intronStart;
		${${$junctIntronRngHsh{$cntg}}{$junctStr}}{"end"} = $intronEnd;
		
		#---multiple $junctStr may exist in HMMSplicer as unique and duplicated individual junctions are collapsed seperately
		if (not exists $junctCtngHsh{$junctStr}) { #---most of the case
			$junctCtngHsh{$junctStr} = $cntg;
			${${$junctRngHsh{$cntg}}{$junctStr}}{"start"} = $bedStart;
			${${$junctRngHsh{$cntg}}{$junctStr}}{"end"} = $bedEnd;
			
			$junctStrdHsh{$junctStr} = $strd;
			$junctReadNumHsh{$junctStr} = $readNum;
			$junctScoreHsh{$junctStr} = $score;

		} else { #---appeared twice, extend the range 
			
			$dupJunctNum++;
			
			my $storedStart = ${${$junctRngHsh{$cntg}}{$junctStr}}{"start"}; 
			my $storedEnd = ${${$junctRngHsh{$cntg}}{$junctStr}}{"end"};
			my $storedScore = $junctScoreHsh{$junctStr};
			my $storedReadNum = $junctReadNumHsh{$junctStr};
			
			$junctReadNumHsh{$junctStr} = $readNum + $storedReadNum;
			$junctScoreHsh{$junctStr} = $score if ($score > $storedScore);;
			
			${${$junctRngHsh{$cntg}}{$junctStr}}{"start"} = $bedStart if ($bedStart < $storedStart);
			${${$junctRngHsh{$cntg}}{$junctStr}}{"end"} = $bedEnd if ($bedEnd > $storedEnd);	

			print TMPLOG "$junctStr appeared in the bed file twice. Collpasing the two BED lines. Score = $storedScore|$score, read=$storedReadNum|$readNum\n";

		}

	}	
	close INFILE;
	close TMPLOG;
	
	my $junctNum = keys %junctStrdHsh;
	
	print "Totally $junctNum junctions have been stored. With $dupJunctNum of them appeared twice and collapsed.\n";

	return (\%junctStrdHsh, \%junctRngHsh, \%junctCtngHsh, \%junctScoreHsh, \%junctReadNumHsh, \%junctIntronRngHsh);
	
}
########################################################################## printCMDLogOrFinishMessage
sub printCMDLogOrFinishMessage {

	my $CMDLogOrFinishMessage = $_[0];

	if ($CMDLogOrFinishMessage eq "CMDLog") {
		#---open a log file if it doesnt exists
		my $scriptNameXext = $0;
		$scriptNameXext =~ s/\.\w+$//;
		open (CMDLOG, ">>$scriptNameXext.cmd.log.txt"); #---append the CMD log file
		open (OUTDIRCMDLOG, ">$outDir/run.cmd.log.txt"); #---make the CMD log file
		my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst)=localtime(time);
		my $runTime = sprintf "%04d-%02d-%02d %02d:%02d", $year+1900, $mon+1,$mday,$hour,$min;
		print CMDLOG "[".$runTime."]\t"."perl $0 ".(join " ", @ARGV)."\n";
		close CMDLOG;
		print OUTDIRCMDLOG "[".$runTime."]\t"."perl $0 ".(join " ", @ARGV)."\n";
		close OUTDIRCMDLOG;
		print "\n=========================================================================\n";
		print "$0 starts running at [$runTime]\n";
		print "=========================================================================\n\n";

	} elsif ($CMDLogOrFinishMessage eq "finishMessage") {
		my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst)=localtime(time);
		my $runTime = sprintf "%04d-%02d-%02d %02d:%02d", $year+1900, $mon+1,$mday,$hour,$min;
		print "\n=========================================================================\n";
		print "$0 finished running at [$runTime]\n";
		print "=========================================================================\n\n";
	}

}
########################################################################## readGffNoAltNoIntronNoCDS
sub readGffNoAltNoIntronNoCDS {#--- read the Gff assuming there no alternative splicing, no intron, one gene one exon, which is the case in NGSGff

	my $gffPath = $_[0];

	#---variables to retun
	my (%strndByGeneHsh, %cntgByGeneHsh, %oneExonRngByCntgByGeneHsh, %multiExonRngByCntgByGeneHsh, %geneCtgryHsh);

	#---tmp variables
	my (%geneByRNAHsh, %IDHsh);

	#---read the gff
	open (INFILE, $gffPath);
	print "Reading $gffPath\n";
	while (my $theLine = <INFILE>) {
		chomp $theLine;
		
		if (($theLine !~ m/^\#|^\@/) and ($theLine !~ m/\tsupercontig\t/)) {
			my @theLineSplt = split (/\t/, $theLine);
			my $cntg = $theLineSplt[0];
			my $geneCategory = $theLineSplt[2];
			my $featureStart = $theLineSplt[3];
			my $featureEnd = $theLineSplt[4];
			my $geneStrd = $theLineSplt[6];
			my $allAttributes = $theLineSplt[8];
			my @allAttributesSplt = split /;/, $allAttributes;
			

			#---get the ID and parent
			my ($unqID, $parent);
			foreach my $theAttribute (@allAttributesSplt) {
				if ($theAttribute =~ m/^ID=/) {$unqID = substr ($theAttribute, index ($theAttribute, "=")+1);}
				if ($theAttribute =~ m/^Parent=/) {$parent = substr ($theAttribute, index ($theAttribute, "=")+1);}
			}

			$IDHsh{$unqID}++;
			die "$unqID is not a unique ID. Make sure all ID attributes are unique\n" if ($IDHsh{$unqID} > 1);
			
			if ($geneCategory eq "gene") {#---gene
				my $geneID = $unqID;
				$strndByGeneHsh{$geneID} = $geneStrd;
				$cntgByGeneHsh{$geneID} = $cntg;

			} elsif ($geneCategory eq "exon") {#---exon, may be exons of alternative transcripts, wiull sort out later
				my $exonID = $unqID;
				my $RNAID = $parent; #---exon has to be coming from a RNA
				my $geneID = $geneByRNAHsh{$RNAID};
				my $cntg = $cntgByGeneHsh{$geneID};

				#---use this variable if single exon is expected;
				${${$oneExonRngByCntgByGeneHsh{$cntg}}{$geneID}}{"start"} = $featureStart;
				${${$oneExonRngByCntgByGeneHsh{$cntg}}{$geneID}}{"end"} = $featureEnd;

				#---use this variable if multiexon is expected;
				${${$multiExonRngByCntgByGeneHsh{$geneID}}{$exonID}}{"start"} = $featureStart;
				${${$multiExonRngByCntgByGeneHsh{$geneID}}{$exonID}}{"end"} = $featureEnd;

			} else {#---can be tRNA, rRNA, mRNA, repRNA, ncRNA
				my $RNAID = $unqID;
				my $geneID = $parent;
				$geneByRNAHsh{$RNAID} = $geneID;
				$geneCtgryHsh{$geneID} = $geneCategory;
			}
		}#---end of if (($theLine !~ m/^\#|^\@/) and ($theLine !~ m/\tsupercontig\t/)) {
	}#---end of while (my $theLine = <INFILE>)
	close INFILE;
	
	#--- check is there any ftur stored;
	my $fturNum =  keys %geneCtgryHsh;
	die "No feature stored. Likely due to conflicts in the contig names or wrong cntgFltr options\n" if ($fturNum == 0);
	
	print "Finished reading.\n";

	my $geneNum = keys %cntgByGeneHsh;
	print "The exon ranges of $geneNum genes have been stored.\n";

	return (\%strndByGeneHsh, \%oneExonRngByCntgByGeneHsh, \%cntgByGeneHsh);
}
########################################################################## defineOverlappingJunctBEDRng
sub defineOverlappingJunctBEDRng {

	my %SSOvrlpJunctBEDRngHsh = %{$_[0]};
	my %junctRngHsh = %{$_[1]};
	my %junctCntgHsh = %{$_[2]};
	my %junctStrdHsh = %{$_[3]};
	
	#--- define the clusters
	my ($junctBEDRngClusterNum, $junctBEDRngClusterNameByJunctHsh_ref, $junctBEDRngClusterJunctHsh_ref) = recursivelyFindJunctionCluster(\%SSOvrlpJunctBEDRngHsh, "junction BED range", "bc_");
	my %junctBEDRngClusterJunctHsh = %{$junctBEDRngClusterJunctHsh_ref};
	my %junctBEDRngClusterNameByJunctHsh = %{$junctBEDRngClusterNameByJunctHsh_ref};

	#--- assign a cluster to the orphan
	my $orphanBEDRngNum = 0;
	
	foreach my $junctStr (sort {$a cmp $b} keys %junctStrdHsh) {
		$orphanBEDRngNum++;
		my $clusterName = "bo_".$orphanBEDRngNum;
		if (not exists $junctBEDRngClusterNameByJunctHsh{$junctStr}) { #--it is an orphan
			$junctBEDRngClusterNameByJunctHsh{$junctStr} = $clusterName;
			push @{$junctBEDRngClusterJunctHsh{$clusterName}}, $junctStr;
		}
	}
	
	print "Defining the junction BED range clusters\n";
	#---define the cluster range
	my (%junctBEDClusterRngHsh, %junctBEDClusterStrdHsh, %junctBEDClusterCntgHsh);
	foreach my $clusterName (sort {$a cmp $b} keys %junctBEDRngClusterJunctHsh) {
		my ($cntg, $strd, @tmpBoundAry);
		foreach my $junctStr (@{$junctBEDRngClusterJunctHsh{$clusterName}})	{
			$cntg = $junctCntgHsh{$junctStr};
			$strd = $junctStrdHsh{$junctStr};
			push @tmpBoundAry, ${${$junctRngHsh{$cntg}}{$junctStr}}{"start"};
			push @tmpBoundAry, ${${$junctRngHsh{$cntg}}{$junctStr}}{"end"};
		}
		my @tmpBoundSortedAry = sort {$a <=> $b} @tmpBoundAry;
		
		${${$junctBEDClusterRngHsh{$cntg}}{$clusterName}}{"start"} = $tmpBoundSortedAry[0];
		${${$junctBEDClusterRngHsh{$cntg}}{$clusterName}}{"end"} = $tmpBoundSortedAry[-1];
		$junctBEDClusterStrdHsh{$clusterName} = $strd;
		$junctBEDClusterCntgHsh{$clusterName} = $cntg;
		#print $clusterName."\t".${${$junctBEDClusterRngHsh{$cntg}}{$clusterName}}{"start"}."\t".${${$junctBEDClusterRngHsh{$cntg}}{$clusterName}}{"end"}."\t".$strd."\t".$cntg."\n";
	}
	
	return (\%junctBEDRngClusterJunctHsh, \%junctBEDRngClusterNameByJunctHsh, \%junctBEDClusterRngHsh, \%junctBEDClusterStrdHsh, \%junctBEDClusterCntgHsh);
}
########################################################################## recursivelyFindJunctionCluster
sub recursivelyFindJunctionCluster {
	
	#---works for both intron range and BED range
	
	our %SSOvrlpRngHsh = %{$_[0]};
	my $reportString = $_[1];
	my $clusterIDTag = $_[2];
	
	#---var to return
	our %clusterJunctRngHsh = ();
	our %junctClusterNameHsh = (); 

	#---check each transfrag
	my $i = my $clusterNum = 0;
	
	print "Start clustering the $reportString.\n";
	
	foreach my $refJunct (keys %SSOvrlpRngHsh) {
	
		if (not exists $junctClusterNameHsh{$refJunct}) {
			$clusterNum++; $i++;
			my $clusterName = $clusterIDTag.$clusterNum;
			my @intitialJunctRecursiveAry = ($refJunct);
			my $junctClusterNameHsh_ref = \%junctClusterNameHsh;
			my $clusterJunctRngHsh_ref = \%clusterJunctRngHsh;
	
			recursiveClusterExtension(\@intitialJunctRecursiveAry, $clusterName);
		
			if ($i eq 1000) {
				$i = 0;	
				print $clusterNum." $reportString clusters identified.\n";
			}
		}
	}

	return ($clusterNum, \%junctClusterNameHsh, \%clusterJunctRngHsh);
	
	####################################################recursiveClusterExtension####################################################
	sub recursiveClusterExtension {
	
		my @recusiveJunctInAry = @{$_[0]};
		my $r_clusterName = $_[1];

		my @recusiveJunctOutAry;
		
		foreach my $recusiveJunctIn (@recusiveJunctInAry) {

			if (not exists $junctClusterNameHsh{$recusiveJunctIn}) { 
				$junctClusterNameHsh{$recusiveJunctIn} = $r_clusterName;
				push @{$clusterJunctRngHsh{$r_clusterName}}, $recusiveJunctIn;

				foreach my $r_refJunct (keys %{$SSOvrlpRngHsh{$recusiveJunctIn}}) {#---check each 
					
					if (not exists $junctClusterNameHsh{$r_refJunct}) { 
						$junctClusterNameHsh{$r_refJunct} = $r_clusterName;
						push @{$clusterJunctRngHsh{$r_clusterName}}, $r_refJunct;
				
						#----look for the junctions hit this junction
						foreach my $r_qryJunct (keys %{$SSOvrlpRngHsh{$r_refJunct}}) {
					
							#---only if the hit is recusiveJunctIn and the hit has not been reached before
							if (($r_qryJunct ne $recusiveJunctIn) and (not exists $junctClusterNameHsh{$r_qryJunct}) and ($r_qryJunct ne $r_refJunct)) {
								push @recusiveJunctOutAry, $r_qryJunct;
							}
						}
					}
				}
			}
		}

		#---if there are junct hit, call itself again
		recursiveClusterExtension(\@recusiveJunctOutAry, $r_clusterName) if (@recusiveJunctOutAry > 0);
	}	
}
########################################################################## collpaseOverlappingTransfragJunctionBEDRngCluster
sub collpaseOverlappingTransfragJunctionBEDRngCluster {
	
	#---common variables for all innner subroutines
	our %XSHitByTrnsfgHsh = %{$_[0]};
	our %XSHitByJunctHsh = %{$_[1]};
	our %trnsfgStrdHsh = %{$_[2]};
	our %trnsfgRngHsh = %{$_[3]};
	our %trnsfgCntgHsh = %{$_[4]};
	our %junctStrdHsh = %{$_[5]};
	our %junctRngHsh = %{$_[6]};
	our %junctCntgHsh = %{$_[7]};
	our $clusterJunctTrnsfgHsh_ref = $_[8];
	our $clusterNum = $_[9];
	our $clusterStrand = $_[10];

	#---use recursiveClusterExtension within findTrnsfgJunctionCluster to collapse all overlapping Jnct and Transfg
	findTrnsfgJunctionCluster();

	our %clusterJunctTrnsfgHsh = %{$clusterJunctTrnsfgHsh_ref};

	return (\%clusterJunctTrnsfgHsh, $clusterNum);
	
	#--- A subroutine contains another recursive routine to find the clusters
	####################################################findTrnsfgJunctionCluster####################################################
	sub findTrnsfgJunctionCluster {

		#---var to return
		our %clusterJunctTrnsfgHsh = %{$clusterJunctTrnsfgHsh_ref};
		our %junctClusterNumHsh = (); 
		our %trnsfgClusterNumHsh = ();
		#---check each transfrag
		my $i = 0;
		
		print "Start clustering the transfrags and $clusterStrand junctions\n";
		
		foreach my $trnsfg (keys %XSHitByTrnsfgHsh) {
		
			if (not exists $trnsfgClusterNumHsh{$trnsfg}) {
				$clusterNum++; $i++;
				${$clusterJunctTrnsfgHsh{"JT_".$clusterNum}}{"strand"} = $clusterStrand;
				${$clusterJunctTrnsfgHsh{"JT_".$clusterNum}}{"cntg"} = $trnsfgCntgHsh{$trnsfg};
				my @intitialTrnsfrgRecursiveAry = ($trnsfg);
				my $trnsfgClusterNumHsh_ref = \%trnsfgClusterNumHsh;
				my $junctClusterNumHsh_ref = \%junctClusterNumHsh;
				my $clusterJunctTrnsfgHsh_ref = \%clusterJunctTrnsfgHsh;
		
				#($trnsfgClusterNumHsh_ref, $junctClusterNumHsh_ref, $clusterJunctTrnsfgHsh_ref) = recursiveClusterExtension($intitialTrnsfrgRecursiveAry_ref, $trnsfgClusterNumHsh_ref, $junctClusterNumHsh_ref, $clusterNum, $clusterJunctTrnsfgHsh_ref);
				recursiveTrnsfgJunctionClusterExtension(\@intitialTrnsfrgRecursiveAry, $clusterNum);
			
				if ($i eq 500) {
					$i = 0;	
					print $clusterNum." clusters identified.\n";
				}
			}
		}

		my $collapsedJunctNum = keys %junctClusterNumHsh;
		my $collapsedTrnsfgNum = keys %trnsfgClusterNumHsh;
		
		print "Totally ".$clusterNum." clusters identified. $collapsedJunctNum junctions and $collapsedTrnsfgNum transfrags were clustered.\n";
	
		####################################################recursiveClusterExtension####################################################
		sub recursiveTrnsfgJunctionClusterExtension {
		
			my @recusiveTrnsfgInAry = @{$_[0]};
			my $r_clusterNum = $_[1];
	
			my @recusiveTrnsfgOutAry;
			
			foreach my $recusiveTrnsfgIn (@recusiveTrnsfgInAry) {
	
				if (not exists $trnsfgClusterNumHsh{$recusiveTrnsfgIn}) { 
					$trnsfgClusterNumHsh{$recusiveTrnsfgIn} = $r_clusterNum;
					push @{${$clusterJunctTrnsfgHsh{"JT_".$r_clusterNum}}{"trnsfg"}}, $recusiveTrnsfgIn;
	
					foreach my $r_junctBEDCluster (keys %{$XSHitByTrnsfgHsh{$recusiveTrnsfgIn}}) {#---check each 
						
						if (not exists $junctClusterNumHsh{$r_junctBEDCluster}) { 
							$junctClusterNumHsh{$r_junctBEDCluster} = $r_clusterNum;
							push @{${$clusterJunctTrnsfgHsh{"JT_".$r_clusterNum}}{"junctBEDCluster"}}, $r_junctBEDCluster;
					
							#----look for the tranfrag that hit this junction
							foreach my $junctHitTrnsfg (keys %{$XSHitByJunctHsh{$r_junctBEDCluster}}) {
						
								#---only if the hit is recusiveTrnsfgIn and the hit has not been reached before
								if (($junctHitTrnsfg ne $recusiveTrnsfgIn) and (not exists $trnsfgClusterNumHsh{$junctHitTrnsfg})) {
									push @recusiveTrnsfgOutAry, $junctHitTrnsfg;
								}
							}
						}
					}
				}
			}
	
			#---if there are trnsfg hit, call itself again
			recursiveTrnsfgJunctionClusterExtension(\@recusiveTrnsfgOutAry, $r_clusterNum) if (@recusiveTrnsfgOutAry > 0);
		}	
	}
}
########################################################################## convertAndPrintJunctTransfragCluster
sub convertAndPrintJunctTransfragCluster {

	my %clusterJunctTrnsfgHsh = %{$_[0]};
	my %junctBEDRngClusterJunctHsh = %{$_[1]};

	foreach my $cluster (sort {$a cmp $b} keys %clusterJunctTrnsfgHsh) {
		foreach my $BEDCluster (@{${$clusterJunctTrnsfgHsh{$cluster}}{"junctBEDCluster"}}) {
			foreach my $junctStr (@{$junctBEDRngClusterJunctHsh{$BEDCluster}}) {
				push @{${$clusterJunctTrnsfgHsh{$cluster}}{"junct"}}, $junctStr;
			}
		}
		delete ${$clusterJunctTrnsfgHsh{$cluster}}{"junctBEDCluster"};
	}
	
	open (TMPLOG, ">$outDir/log/transfragJunctionCluster.txt");
	foreach my $cluster (sort {$a cmp $b} keys %clusterJunctTrnsfgHsh) {
		print TMPLOG $cluster;
		print TMPLOG "\t".${$clusterJunctTrnsfgHsh{$cluster}}{"strand"};
		foreach my $trnsfg (@{${$clusterJunctTrnsfgHsh{$cluster}}{"trnsfg"}}) {print TMPLOG "\t".$trnsfg;}
		foreach my $junctStr (@{${$clusterJunctTrnsfgHsh{$cluster}}{"junct"}}) {print TMPLOG "\t".$junctStr;}
		print TMPLOG "\n";
	}
	close TMPLOG;
}
########################################################################## scanJunctionIntronBoundCoverage
sub scanJunctionIntronBoundCoverage {#----scan for the coverage flanking the intron bounds and record the splicing ratio 

	my %junctIntronRngHsh = %{$_[0]};
	
	#---var to return
	my %junctIntronBoundCovHsh = ();

	#---read or calculate the intron bound coverage
	if ($intronBoundCovPath ne "no") {#---read
		
		print "Reading $intronBoundCovPath for intron bound coverage\n";
		
		open (INTRNBOUNDCOV, "$intronBoundCovPath")|| die "Can't read $intronBoundCovPath :$!\n";
		
		while (my $theLine = <INTRNBOUNDCOV>) {
			chomp $theLine;
			my @theLineSplt = split /\t/, $theLine;
			my $junctStr = $theLineSplt[0];
			${${$junctIntronBoundCovHsh{$junctStr}}{"start"}}{"+"} = $theLineSplt[1];
			${${$junctIntronBoundCovHsh{$junctStr}}{"start"}}{"-"} = $theLineSplt[2];
			${${$junctIntronBoundCovHsh{$junctStr}}{"end"}}{"+"} = $theLineSplt[3];
			${${$junctIntronBoundCovHsh{$junctStr}}{"end"}}{"-"} = $theLineSplt[4];
		}
		close INTRNBOUNDCOV;
		
	} else {#calculate
		
		open (PILEUPFILE, "$pileupPath");
		open (INTRNBOUNDCOV, ">$outDir/cov/intronBoundCov.txt");
	
		print "Reading $pileupPath.\n";
	
		my %tmpCovByPosHsh = ();
		my $flankSize = 30;
		
		#--get the first line
		my $theCurntLine;
		while ($theCurntLine = <PILEUPFILE>) {
			if ($theCurntLine !~ m/^[\@|\#]/) {#---ignore comments
				chomp $theCurntLine;
				last;
			}
		}
		
		#--go through the rest of the file
		while (my $theNextLine = <PILEUPFILE>) {
			next if ($theNextLine =~ m/^[\@|\#]/);
			chomp $theNextLine;
		
			my @theNextLineSplt = split /\t/, $theNextLine;
			my $nextCtng = $theNextLineSplt[0];
	
			my @theCurntLineSplt = split /\t/, $theCurntLine;
			my $curntCtng = $theCurntLineSplt[0];
	
			my $pos = $theCurntLineSplt[1];
			my $plusCov = $theCurntLineSplt[3];
			my $minusCov = $theCurntLineSplt[4];
		
			#---store the cov if plus or minus strand > 0;
			if (($plusCov > 0) or ($minusCov > 0)) {
				${$tmpCovByPosHsh{$pos}}{"+"} = $plusCov;
				${$tmpCovByPosHsh{$pos}}{"-"} = $minusCov;
			}
			
			if (eof(PILEUPFILE)) {#---last line of the file
				my $pos = $theNextLineSplt[1];
				my $plusCov = $theNextLineSplt[3];
				my $minusCov = $theNextLineSplt[4];
	
				#---store the cov if plus or minus strand > 0;
				if (($plusCov > 0) or ($minusCov > 0)) {
					${$tmpCovByPosHsh{$pos}}{"+"} = $plusCov;
					${$tmpCovByPosHsh{$pos}}{"-"} = $minusCov;
				}
			}
			
			#---change contig or end of file
			if (($curntCtng ne $nextCtng) or (eof(PILEUPFILE))) {
				
				print "Scanning intron bound coverage in contig $curntCtng.\n";
				
				foreach my $junctStr (keys %{$junctIntronRngHsh{$curntCtng}}) {
					
					my $intronStart = ${${$junctIntronRngHsh{$curntCtng}}{$junctStr}}{"start"};
					my $intronEnd = ${${$junctIntronRngHsh{$curntCtng}}{$junctStr}}{"end"};
					
					my $startPlusCovSum = my $startMinusCovSum = my $endPlusCovSum = my $endMinusCovSum = 0;
					
					#---intron start position flank region cov
					for my $pos (($intronStart)..($intronStart + $flankSize)) {
						if (exists $tmpCovByPosHsh{$pos}) {
							$startPlusCovSum += ${$tmpCovByPosHsh{$pos}}{"+"};
							$startMinusCovSum += ${$tmpCovByPosHsh{$pos}}{"-"};
						}
					}
					
					${${$junctIntronBoundCovHsh{$junctStr}}{"start"}}{"+"} = sprintf "%.02f", $startPlusCovSum/$flankSize;
					${${$junctIntronBoundCovHsh{$junctStr}}{"start"}}{"-"} = sprintf "%.02f", $startMinusCovSum/$flankSize;
	
					#---intron end position flank region cov
					for my $pos (($intronEnd - $flankSize)..$intronEnd) {
						if (exists $tmpCovByPosHsh{$pos}) {
							$endPlusCovSum += ${$tmpCovByPosHsh{$pos}}{"+"};
							$endMinusCovSum += ${$tmpCovByPosHsh{$pos}}{"-"};
						}
					}
					
					${${$junctIntronBoundCovHsh{$junctStr}}{"end"}}{"+"} = sprintf "%.02f", $endPlusCovSum/$flankSize;
					${${$junctIntronBoundCovHsh{$junctStr}}{"end"}}{"-"} = sprintf "%.02f", $endMinusCovSum/$flankSize;
	
					print INTRNBOUNDCOV $junctStr."\t";
					print INTRNBOUNDCOV ${${$junctIntronBoundCovHsh{$junctStr}}{"start"}}{"+"}."\t";
					print INTRNBOUNDCOV ${${$junctIntronBoundCovHsh{$junctStr}}{"start"}}{"-"}."\t";
					print INTRNBOUNDCOV ${${$junctIntronBoundCovHsh{$junctStr}}{"end"}}{"+"}."\t";
					print INTRNBOUNDCOV ${${$junctIntronBoundCovHsh{$junctStr}}{"end"}}{"-"}."\n";
	
				} #---end of foreach my $junctStr (keys %{$junctIntronRngHsh{$curntCtng}})
				
				%tmpCovByPosHsh = (); #---empty the hash
	
			} #---end of if (($curntCtng ne $nextCtng) or (eof(PILEUPFILE))) {
			
			$theCurntLine = $theNextLine;
			
		}#---end of while (my $theNextLine = <PILEUPFILE>) {
	
		close PILEUPFILE;
		close INTRNBOUNDCOV;
	}
	
	return (\%junctIntronBoundCovHsh);
	
}
########################################################################## defineIntronClusters
sub defineIntronClusters {
	
	our %SSOvrlpJunctIntronHsh = %{$_[0]};
	our %junctStrdHsh = %{$_[1]};
	our %junctScoreHsh = %{$_[2]};
	our %junctReadNumHsh = %{$_[3]};
	our %junctIntronBoundCovHsh = %{$_[4]};
	our %junctIntronRngHsh = %{$_[5]};
	our %XSOvrlpJunctHsh = %{$_[6]};
	
	#---define splcing raio and the major junctions
	my ($junctSplicingRatioHsh_ref, $majorJunctAvgSplcingRatioHsh_ref) = defineMajorJunctions();

	#---define superJunctions
	my $superJunctHsh_ref = defineSuperJunctions();
	
	#---define overlapping intron clusters excluding superJunctions
	my ($clusterAllJunctHsh_ref, $junctClusterNameByJunctHsh_ref, $junctClusterInfoHsh_ref, $junctClusterSSOverlapHsh_ref) = defineIntronCluster($superJunctHsh_ref, $majorJunctAvgSplcingRatioHsh_ref);
	
	#---define the exon skipping events
	my ($totalExonSkippingTypeHsh_ref, $majorExonSkippingTypeHsh_ref, $exactExonSkippingClusterHsh_ref) = defineExonSkipping($superJunctHsh_ref, $majorJunctAvgSplcingRatioHsh_ref, $junctClusterNameByJunctHsh_ref);
	
	#---define the alternative 5' and 3' splicing sites
	my ($splicingSiteDiffAllHsh_ref, $splicingSiteDiffHshByClusterHsh_ref) = defineAltSplicingSiteWithinCluster($clusterAllJunctHsh_ref, $majorJunctAvgSplcingRatioHsh_ref);

	#---define the junctions on the prominent isoform
	my $junctOnProminentIsofmHsh_ref; 
	($junctClusterInfoHsh_ref, $junctOnProminentIsofmHsh_ref) = defineJunctionsOfProminentIsoform($junctClusterNameByJunctHsh_ref, $junctClusterInfoHsh_ref, $junctClusterSSOverlapHsh_ref);
	my %junctClusterInfoHsh = %{$junctClusterInfoHsh_ref}; $junctClusterInfoHsh_ref = \%junctClusterInfoHsh;
	
	return ($junctSplicingRatioHsh_ref, $majorJunctAvgSplcingRatioHsh_ref, $superJunctHsh_ref, $clusterAllJunctHsh_ref, $junctClusterNameByJunctHsh_ref, $junctClusterInfoHsh_ref, $junctClusterSSOverlapHsh_ref, $totalExonSkippingTypeHsh_ref, $majorExonSkippingTypeHsh_ref, $exactExonSkippingClusterHsh_ref, $splicingSiteDiffAllHsh_ref, $splicingSiteDiffHshByClusterHsh_ref, $junctOnProminentIsofmHsh_ref);
	
	################################################## defineMajorJunctions  #########################################################################
	sub defineMajorJunctions {
	
		print "Defining major junctions.\n";

		#---define the major junctions;
		my $histogramCutOff = 20;
		my (@endSplicingRatioToPlotAry, @startSplicingRatioToPlotAry, %startSplicingRatioToPlotHsh, %endSplicingRatioToPlotHsh, %junctSplicingRatioHsh, %majorJunctAvgSplcingRatioHsh);
		
		open (TMPLOG, ">$outDir/log/junctIntronBoundCov.M$minSplicingRatio.log.txt");
		foreach my $junctStr (keys %junctIntronBoundCovHsh) {
			my $supportReadNum = $junctReadNumHsh{$junctStr};		
			my $startBothSum =  ${${$junctIntronBoundCovHsh{$junctStr}}{"start"}}{"+"} + ${${$junctIntronBoundCovHsh{$junctStr}}{"start"}}{"-"};
			my $endBothSum = ${${$junctIntronBoundCovHsh{$junctStr}}{"end"}}{"+"} + ${${$junctIntronBoundCovHsh{$junctStr}}{"end"}}{"-"};
			my $startSplicingRatio = my $endSplicingRatio = 99999;
			$startSplicingRatio = sprintf "%.01f", $supportReadNum/$startBothSum if ($startBothSum > 0); 
			$endSplicingRatio = sprintf "%.01f", $supportReadNum/$endBothSum if ($endBothSum > 0);
			${$junctSplicingRatioHsh{$junctStr}}{"start"} = $startSplicingRatio;
			${$junctSplicingRatioHsh{$junctStr}}{"end"} = $endSplicingRatio;
	
			my $startSplicingRatioToPlot = $startSplicingRatio;
			my $endSplicingRatioToPlot = $endSplicingRatio; 
	
			$startSplicingRatioToPlot = $histogramCutOff if ($startSplicingRatioToPlot > $histogramCutOff);
			$endSplicingRatioToPlot = $histogramCutOff if ($endSplicingRatioToPlot > $histogramCutOff);
	
			$startSplicingRatioToPlotHsh{$startSplicingRatioToPlot}++;
			$endSplicingRatioToPlotHsh{$endSplicingRatioToPlot}++;
	
			push @endSplicingRatioToPlotAry, $endSplicingRatioToPlot;
			push @startSplicingRatioToPlotAry, $startSplicingRatioToPlot;
			
			my $majorJunct = "no";
	
			if (($startSplicingRatio >= $minSplicingRatio) and ($endSplicingRatio >= $minSplicingRatio) and ($supportReadNum >= $minMajorSupportRead)) {
				$majorJunct = "yes";
				$majorJunctAvgSplcingRatioHsh{$junctStr} = sprintf "%0.2f", ($endSplicingRatio + $startSplicingRatio)/2;
			}
			print TMPLOG $junctStr."\t".$startBothSum."\t".$endBothSum."\t".$startSplicingRatio."\t".$endSplicingRatio."\t".$majorJunct."\n";
		}
		close TMPLOG;
		
		my $totalJunctNum = keys %junctIntronBoundCovHsh;
		my $majorJunctNum = keys %majorJunctAvgSplcingRatioHsh;
		
		print $majorJunctNum." out of $totalJunctNum junctions having boundary splicing ratio >= $minSplicingRatio and supporting read >= $minMajorSupportRead.\n";
		
		#---plot the splicing ratio info
		my $cumulStartSplicingRatioToPlotHsh_ref = generateCumulativePctHash(\%startSplicingRatioToPlotHsh, "a");
		my $cumulEndSplicingRatioToPlotHsh_ref = generateCumulativePctHash(\%endSplicingRatioToPlotHsh, "a");
		GNUPLOTAryHistogram(\@startSplicingRatioToPlotAry, "$outDir/plotPdf/startSplicingRatioCount.pdf", "$outDir/plotData/startSplicingRatioCount.dat", "splicing ratio", "frequency", "splicing ratio of intron start bound", "no");
		GNUPLOTAryHistogram(\@endSplicingRatioToPlotAry, "$outDir/plotPdf/endSplicingRatioCount.pdf", "$outDir/plotData/endSplicingRatioCount.dat", "splicing ratio", "frequency", "splicing ratio of intron end bound", "no");
		GNUPlotXYScatterWithLines($cumulEndSplicingRatioToPlotHsh_ref, "$outDir/plotPdf/cumulEndSplicingRatioCount.pdf", "$outDir/plotData/cumulEndSplicingRatioCount.dat", "smaller than splicing ratio", "percentage", "linear", "linear", "splicing ratio of intron end bound", "no");
		GNUPlotXYScatterWithLines($cumulStartSplicingRatioToPlotHsh_ref, "$outDir/plotPdf/cumulStartSplicingRatioCount.pdf", "$outDir/plotData/cumulStartSplicingRatioCount.dat", "smaller than splicing ratio", "percentage", "linear", "linear", "splicing ratio of intron start bound", "no");
		
		return (\%junctSplicingRatioHsh, \%majorJunctAvgSplcingRatioHsh);
	}
	################################################## defineSuperJunctions  #########################################################################
	sub defineSuperJunctions {#---A superjunction is defined as a junctiuon that is overlapping with two or more non-overlapping junctions
		
		print "Defining superJunctions.\n";

		my %superJunctHsh;
		
		open (TMPLOG, ">$outDir/log/allOverlappingJunctions.txt");
		foreach my $refJunctStr (sort {$a cmp $b} keys %SSOvrlpJunctIntronHsh) {
			print TMPLOG $refJunctStr;
			foreach my $qryJunctStr (sort {$a cmp $b} keys %{$SSOvrlpJunctIntronHsh{$refJunctStr}}) {
				print TMPLOG "\t".$qryJunctStr;
				foreach my $refOverlapJunct (keys %{$SSOvrlpJunctIntronHsh{$refJunctStr}}) {
					$superJunctHsh{$refJunctStr}++ if ((not exists ${$SSOvrlpJunctIntronHsh{$refOverlapJunct}}{$qryJunctStr}) and ($qryJunctStr ne $refOverlapJunct));
				}
			}
			print TMPLOG "\n";
		}
		close TMPLOG;
		
		print "Defining overlapping superJunctions\n";
		#----construct an hash to contain overlapping superJunctions
		my %SSOvrlpingSuperJunctHsh;
		foreach my $superJunct (sort {$a cmp $b} keys %superJunctHsh) {
			foreach my $qryJunctStr (sort {$a cmp $b} keys %{$SSOvrlpJunctIntronHsh{$superJunct}}) {
				if (exists $superJunctHsh{$qryJunctStr}) {
					${$SSOvrlpingSuperJunctHsh{$superJunct}}{$qryJunctStr} = ${$SSOvrlpJunctIntronHsh{$superJunct}}{$qryJunctStr};
				}
			}
		}

		open (TMPLOG, ">$outDir/log/overlappingSuperJunctions.txt");
		foreach my $refJunctStr (sort {$a cmp $b} keys %SSOvrlpingSuperJunctHsh) {
			print TMPLOG $refJunctStr;
			foreach my $qryJunctStr (sort {$a cmp $b} keys %{$SSOvrlpingSuperJunctHsh{$refJunctStr}}) {
				print TMPLOG "\t".$qryJunctStr;
			}
			print TMPLOG "\n";
		}
		close TMPLOG;
		
		open (TMPLOG, ">$outDir/log/allSuperJunctionsOverlappings.txt");
		foreach my $superJunct (sort {$a cmp $b} keys %superJunctHsh) {
			print TMPLOG $superJunct;
			foreach my $qryJunctStr (sort {$a cmp $b} keys %{$SSOvrlpJunctIntronHsh{$superJunct}}) {
				print TMPLOG "\t".$qryJunctStr;
			}
			print TMPLOG "\n";
		}
		close TMPLOG;
		
		return \%superJunctHsh;
	}
	############################################## defineIntronCluster ########################################################
	sub defineIntronCluster {
	
		print "Defining super-/inferior- junction clusters.\n";
	
		my %superJunctHsh = %{$_[0]};
		my %majorJunctAvgSplcingRatioHsh = %{$_[1]};
		
		#----construct an hash to contain overlapping junctions without superJunctions
		my (%SSOvrlpingJunctXSuperHsh, %SSOvrlpingJunctOnlySuperHsh);
		foreach my $refJunctStr (sort {$a cmp $b} keys %SSOvrlpJunctIntronHsh) {
			foreach my $qryJunctStr (sort {$a cmp $b} keys %{$SSOvrlpJunctIntronHsh{$refJunctStr}}) {
				if ((not exists $superJunctHsh{$refJunctStr}) and (not exists $superJunctHsh{$qryJunctStr})) {
					#---non-superJunction overlpping hash
					${$SSOvrlpingJunctXSuperHsh{$refJunctStr}}{$qryJunctStr} = ${$SSOvrlpJunctIntronHsh{$refJunctStr}}{$qryJunctStr};

				} elsif ((exists $superJunctHsh{$refJunctStr}) and (exists $superJunctHsh{$qryJunctStr})) {
					#---Only superJunction overlpping hash
					${$SSOvrlpingJunctOnlySuperHsh{$refJunctStr}}{$qryJunctStr} = ${$SSOvrlpJunctIntronHsh{$refJunctStr}}{$qryJunctStr};
				}
			}
		}
		
		my (%clusterAllJunctHsh, %junctClusterNameByJunctHsh);
		
		#---recursively find the inferior clusters
		my ($inferiorClusterNum, $inferiorClusterNameByJunctHsh, $inferiorClusterJunctHsh_ref) = recursivelyFindJunctionCluster(\%SSOvrlpingJunctXSuperHsh, "inferior junction intron", "ic_");
		my %inferiorClusterJunctHsh = %{$inferiorClusterJunctHsh_ref};

		print "$inferiorClusterNum inferior junction clusters were identified.\n";

		#----transfer the inferior results to the big hash
		foreach my $clusterName (keys %inferiorClusterJunctHsh) {
			foreach my $junctStr (@{$inferiorClusterJunctHsh{$clusterName}}) {
				$junctClusterNameByJunctHsh{$junctStr}  = $clusterName;
				push @{$clusterAllJunctHsh{$clusterName}}, $junctStr;
			}
			delete $inferiorClusterJunctHsh{$clusterName};
		}
		
		#---recursively find the super clusters
		my ($superClusterNum, $superClusterNameByJunctHsh, $superClusterJunctHsh_ref) = recursivelyFindJunctionCluster(\%SSOvrlpingJunctOnlySuperHsh, "super junction intron", "sc_");
		

		my %superClusterJunctHsh = %{$superClusterJunctHsh_ref};

		print "$superClusterNum super junction clusters were identified.\n";

		#----transfer the inferior results to the big hash
		foreach my $clusterName (keys %superClusterJunctHsh) {
			foreach my $junctStr (@{$superClusterJunctHsh{$clusterName}}) {
				$junctClusterNameByJunctHsh{$junctStr}  = $clusterName;
				push @{$clusterAllJunctHsh{$clusterName}}, $junctStr;
			}
			delete $superClusterJunctHsh{$clusterName};
		}

		#---assigne numbers for orphans also
		my $inferiorOphranNum = my $superOphranNum = 0;
		foreach my $junctStr (sort {$a cmp $b} keys %junctStrdHsh) {
			my $strand = $junctStrdHsh{$junctStr};
			if (not exists $junctClusterNameByJunctHsh{$junctStr}) {#---not in the clusters, no overlapping
				if (not exists $superJunctHsh{$junctStr}) {#---non superJunction
					$inferiorOphranNum++;
					my $clusterName = "io_".$inferiorOphranNum.$strand;
					push @{$clusterAllJunctHsh{$clusterName}}, $junctStr;
					$junctClusterNameByJunctHsh{$junctStr}  = $clusterName;
				} else {
					$superOphranNum++;
					my $clusterName = "so_".$superOphranNum.$strand;
					push @{$clusterAllJunctHsh{$clusterName}}, $junctStr;
					$junctClusterNameByJunctHsh{$junctStr}  = $clusterName;
				}
			}
		}

		print "$superOphranNum super junction orphans were identified.\n";
		print "$inferiorOphranNum inferior junction orphans were identified.\n";

		#---defining overlapping junction clusters
		print "Defining strand specific overlapping of clusters\n";
		my %junctClusterSSOverlapHsh;
		foreach my $refJunctStr (sort {$a cmp $b} keys %SSOvrlpJunctIntronHsh) {
			my $refCluster = $junctClusterNameByJunctHsh{$refJunctStr};
			foreach my $qryJunctStr (sort {$a cmp $b} keys %{$SSOvrlpJunctIntronHsh{$refJunctStr}}) {
				my  $qryCluster = $junctClusterNameByJunctHsh{$qryJunctStr};
				${$junctClusterSSOverlapHsh{$refCluster}}{$qryCluster}++ if ($refCluster ne $qryCluster);
			}
		}

		my (%junctClusterAntisenseOverlapHsh, %junctStrAntisenseOverlapHsh);
		
		print "Defining overlapping sense/Antisense jucntions clusters\n";
		
		foreach my $refJunctStr (keys %XSOvrlpJunctHsh) {
			my $refStrand = $junctStrdHsh{$refJunctStr};
			my $refCluster = $junctClusterNameByJunctHsh{$refJunctStr};
			
			#---skip this refJunct if it is not strand specific
			next if (($refStrand ne "+") and  ($refStrand ne "-"));
			foreach my $qryJunctStr (keys %{$XSOvrlpJunctHsh{$refJunctStr}}) {
				my $qryStrand = $junctStrdHsh{$qryJunctStr};
				my $qryCluster = $junctClusterNameByJunctHsh{$qryJunctStr};

				#---skip this refJunct if it is not strand specific
				next if (($qryStrand ne "+") and  ($qryStrand ne "-"));

				if (($refStrand ne $qryStrand) and ($refCluster ne $qryCluster)) {
					${$junctClusterAntisenseOverlapHsh{$refCluster}}{$qryCluster}++;
					${$junctStrAntisenseOverlapHsh{$refJunctStr}}{$qryJunctStr}++;
				}
			}
		}

		#---get all info of the clusters
		my %junctClusterInfoHsh;
		print "Getting all cluster information.\n";
		foreach my $clusterName (sort {$a cmp $b} keys %clusterAllJunctHsh) {
			${$junctClusterInfoHsh{$clusterName}}{"majorJunctNum"} = 0;
			${$junctClusterInfoHsh{$clusterName}}{"totalJunctNum"} = @{$clusterAllJunctHsh{$clusterName}};
			my $prominentJunctStr;
			my $prominentJunctSupportReadNum = 0;
			${$junctClusterInfoHsh{$clusterName}}{"majorCluster"} = "no";
			foreach my $junctStr (@{$clusterAllJunctHsh{$clusterName}}) {
				${$junctClusterInfoHsh{$clusterName}}{"majorJunctNum"}++ if (exists $majorJunctAvgSplcingRatioHsh{$junctStr});			
				if ($junctReadNumHsh{$junctStr} > $prominentJunctSupportReadNum) { 
					$prominentJunctStr = $junctStr;
					$prominentJunctSupportReadNum = $junctReadNumHsh{$junctStr};
				}
			}
			${$junctClusterInfoHsh{$clusterName}}{"majorCluster"} = "yes" if (exists $majorJunctAvgSplcingRatioHsh{$prominentJunctStr});
			${$junctClusterInfoHsh{$clusterName}}{"prominentJunct"} = $prominentJunctStr;
			${$junctClusterInfoHsh{$clusterName}}{"strand"} = $junctStrdHsh{$prominentJunctStr};
			${$junctClusterInfoHsh{$clusterName}}{"mostSupportReadNum"} = $prominentJunctSupportReadNum;
		}

		#---print the overlapping antisense junctions
		open (TMPLOG, ">$outDir/log/overlappingAntisenseJunctions.txt");
		foreach my $refJunctStr (sort {$a cmp $b} keys %junctStrAntisenseOverlapHsh) {
			my $refCluster = $junctClusterNameByJunctHsh{$refJunctStr};
			my $refReadNum = $junctReadNumHsh{$refJunctStr};
			my $refIsMajor = "no";
			$refIsMajor = "yes" if (exists $majorJunctAvgSplcingRatioHsh{$refJunctStr});
			foreach my $qryJunctStr (sort {$a cmp $b} keys %{$junctStrAntisenseOverlapHsh{$refJunctStr}}) {
				my $qryCluster = $junctClusterNameByJunctHsh{$qryJunctStr};
				my $qryReadNum = $junctReadNumHsh{$qryJunctStr};
				my $qryIsMajor = "no";
				$qryIsMajor = "yes" if (exists $majorJunctAvgSplcingRatioHsh{$qryJunctStr});
				my $readRatio = sprintf "%.05f", $refReadNum/$qryReadNum;
				print TMPLOG join "\t", ($refJunctStr, $qryJunctStr, $refCluster, $qryCluster, $refReadNum, $qryReadNum, $refIsMajor, $qryIsMajor, $readRatio."\n");
			}
		}
		close TMPLOG;

		#---print all info of the intronClusters
		open (TMPLOG, ">$outDir/log/intronClusters.txt");
		foreach my $clusterName (sort {$a cmp $b} keys %clusterAllJunctHsh) {
			my @tmpOverlappingSSClusterAry = ();
			if (exists $junctClusterSSOverlapHsh{$clusterName}) {
				foreach my $overlappingSSCluster (keys %{$junctClusterSSOverlapHsh{$clusterName}}) {
					push @tmpOverlappingSSClusterAry, $overlappingSSCluster;
				}
			}

			my @tmpOverlappingAntisenseClusterAry = ();
			if (exists $junctClusterAntisenseOverlapHsh{$clusterName}) {
				foreach my $overlappingAntisenseCluster (keys %{$junctClusterAntisenseOverlapHsh{$clusterName}}) {
					push @tmpOverlappingAntisenseClusterAry, $overlappingAntisenseCluster;
				}
			}

			print TMPLOG join "\t", $clusterName, ${$junctClusterInfoHsh{$clusterName}}{"majorJunctNum"}, ${$junctClusterInfoHsh{$clusterName}}{"totalJunctNum"}, ${$junctClusterInfoHsh{$clusterName}}{"majorCluster"}, ${$junctClusterInfoHsh{$clusterName}}{"prominentJunct"}, ${$junctClusterInfoHsh{$clusterName}}{"mostSupportReadNum"}, (join ";", @{$clusterAllJunctHsh{$clusterName}}),(join ";", @tmpOverlappingSSClusterAry), (join ";", @tmpOverlappingAntisenseClusterAry)."\n";
		}
		close TMPLOG;

		return \%clusterAllJunctHsh, \%junctClusterNameByJunctHsh, \%junctClusterInfoHsh, \%junctClusterSSOverlapHsh;
	}
	############################################## findSuperJunctionOverlapCluster ########################################################
	sub defineExonSkipping {
	
		print "Defining exon skipping events.\n";
	
		my %superJunctHsh = %{$_[0]};
		my %majorJunctAvgSplcingRatioHsh = %{$_[1]};
		my %junctClusterNameByJunctHsh = %{$_[2]};
		
		#---hash to return
		my (%totalExonSkippingTypeHsh, %majorExonSkippingTypeHsh, %exactExonSkippingClusterHsh, %superJunctOvrlpClusterHsh);
		
		foreach my $superJunct (sort {$a cmp $b} keys %superJunctHsh) {
			#--- non-super junctions are defined as inferior junctions here
			my $startTotalHit = my $endTotalHit = my $startMajorHit = my $endMajorHit = 0;
			my @superJunctSplt = split /:/, $superJunct;
			my $superJunctStart = $superJunctSplt[1];
			my $superJunctEnd = $superJunctSplt[2];;

			#---check for exact overlapping of start and end of the inferior junctions
			foreach my $inferiorJunctStr (sort {$a cmp $b} keys %{$SSOvrlpJunctIntronHsh{$superJunct}}) {
				my @inferiorJunctStrSplt = split /:/, $inferiorJunctStr;
				my $inferiorJunctStart = $inferiorJunctStrSplt[1];
				my $inferiorJunctEnd = $inferiorJunctStrSplt[2];
					
				my $inferiorClusterNum = $junctClusterNameByJunctHsh{$inferiorJunctStr};
					 
				${$superJunctOvrlpClusterHsh{$superJunct}}{$inferiorClusterNum}++;
					
				if ($superJunctStart eq $inferiorJunctStart) {
					${$exactExonSkippingClusterHsh{$superJunct}}{$inferiorJunctStr}++;
					$startTotalHit ++;
					$startMajorHit ++ if (exists $majorJunctAvgSplcingRatioHsh{$inferiorJunctStr});
				}
				
				if ($superJunctEnd eq $inferiorJunctEnd) {
					${$exactExonSkippingClusterHsh{$superJunct}}{$inferiorJunctStr}++;
					$endTotalHit ++;
					$endMajorHit ++ if (exists $majorJunctAvgSplcingRatioHsh{$inferiorJunctStr});
				}
			} #---end of foreach my $inferiorJunctStr (sort {$a cmp $b} keys %{$SSOvrlpJunctIntronHsh{$superJunct}}) {
			
			if (($startTotalHit > 0) and ($endTotalHit > 0)) {#---exact exonSkipping
				$totalExonSkippingTypeHsh{$superJunct} = "exact";

			} elsif (($startTotalHit > 0) or ($endTotalHit > 0)) {#---partial exonSkpping
				$totalExonSkippingTypeHsh{$superJunct} = "partial";
			
			} elsif (($startTotalHit == 0) or ($endTotalHit  == 0)) {#---inexact exonSkpping
				$totalExonSkippingTypeHsh{$superJunct} = "inexact";
				
			} else {#---unexpected events
				die "Unexpected sincerio of exon skipping events.\n";
			}
			
			if (($startMajorHit > 0) and ($endMajorHit > 0)) {#---exact exonSkipping
				$majorExonSkippingTypeHsh{$superJunct} = "exact";

			} elsif (($startMajorHit > 0) or ($endMajorHit > 0)) {#---partial exonSkpping
				$majorExonSkippingTypeHsh{$superJunct} = "partial";
			
			} elsif (($startMajorHit == 0) or ($endMajorHit  == 0)) {#---inexact exonSkpping
				$majorExonSkippingTypeHsh{$superJunct} = "inexact";
				
			} else {#---unexpected events
				die "Unexpected sincerio of exon skipping events.\n";
			}

		}
		
		open (TMPLOG, ">$outDir/log/totalExonSkippingEvents.txt");
		foreach my $superJunct (sort {$totalExonSkippingTypeHsh{$a} cmp $totalExonSkippingTypeHsh{$b}} keys %totalExonSkippingTypeHsh) {
			print TMPLOG $superJunct."\t".$totalExonSkippingTypeHsh{$superJunct};
			my $skippedClusterNum = keys %{$superJunctOvrlpClusterHsh{$superJunct}};
			foreach my $inferiorClusterNum (sort {$a cmp $b} keys %{$superJunctOvrlpClusterHsh{$superJunct}}) {
				print TMPLOG "\t".$inferiorClusterNum;
			}
			print TMPLOG "\n";
		}
		close TMPLOG;

		open (TMPLOG, ">$outDir/log/majorExonSkippingEvents.txt");
		foreach my $superJunct (sort {$majorExonSkippingTypeHsh{$a} cmp $majorExonSkippingTypeHsh{$b}} keys %majorExonSkippingTypeHsh) {
			print TMPLOG $superJunct."\t".$majorExonSkippingTypeHsh{$superJunct};
			my $skippedClusterNum = keys %{$superJunctOvrlpClusterHsh{$superJunct}};
			foreach my $inferiorClusterNum (sort {$a cmp $b} keys %{$superJunctOvrlpClusterHsh{$superJunct}}) {
				print TMPLOG "\t".$inferiorClusterNum;
			}
			print TMPLOG "\n";
		}
		close TMPLOG;

		return (\%totalExonSkippingTypeHsh, \%majorExonSkippingTypeHsh, \%exactExonSkippingClusterHsh);
		
	}
	############################################## defineAltSplicingSiteWithinCluster ########################################################
	sub defineAltSplicingSiteWithinCluster {
	
		my %clusterAllJunctHsh = %{$_[0]};
		my %majorJunctAvgSplcingRatioHsh = %{$_[1]};

		print "Defining alternative splicing site within intron clusters\n";
		
		#---hashes to return
		my (%splicingSiteDiffAllHsh, %splicingSiteDiffHshByClusterHsh, %splicingSiteDiffMajorHsh);

		foreach my $clusterName (sort {$a cmp $b} keys %clusterAllJunctHsh) {
			
			if (@{$clusterAllJunctHsh{$clusterName}} > 1) { #---clusters with more than 1 junct, not orphan and superJunctions
				my $lastIndex = $#{$clusterAllJunctHsh{$clusterName}};
				
				#---store the Strd of the first junction, to make sure all junct within cluster is on the same strd
				#---since 0 wont appear as $i,so get the index 0 strd first
				my $clusterStrd = $junctStrdHsh{${$clusterAllJunctHsh{$clusterName}}[0]};
				
				#---pairwise comparison of cluster junctions
				for (my $i=$lastIndex; $i>0; $i--) {#---5, 4, 3, 2, 1 for 6 junctions
					my $refJunctStr = ${$clusterAllJunctHsh{$clusterName}}[$i];
					my @refJunctStrSplt = split /:/, $refJunctStr;
					my $refJunctStart = $refJunctStrSplt[1];
					my $refJunctEnd = $refJunctStrSplt[2];
					
					my $refIsMajor = "no";
					$refIsMajor = "yes" if (exists $majorJunctAvgSplcingRatioHsh{$refJunctStr});
					
					die "Intron cluster contain jucntions of diifferent strands. Unexpected scinerio.\n" if ($clusterStrd ne $junctStrdHsh{$refJunctStr});

					for (my $j=$i-1; $j>=0; $j--) {#----4, 3, 2, 1, 0|| 3, 2, 1, 0 || 2, 1, 0 || 1, 0|| 0 for 6 junctions
						my $qryJunctStr = ${$clusterAllJunctHsh{$clusterName}}[$j];
						my @qryJunctStrSplt = split /:/, $qryJunctStr;
						my $qryJunctStart = $qryJunctStrSplt[1];
						my $qryJunctEnd = $qryJunctStrSplt[2];
						
						my $startDiff = abs ($refJunctStart - $qryJunctStart);
						my $endDiff = abs ($refJunctEnd - $qryJunctEnd);
						my $fullDiff = abs (($qryJunctEnd - $qryJunctStart + 1) - ($refJunctEnd - $refJunctStart + 1));
						${$splicingSiteDiffAllHsh{"full"}}{$fullDiff}++;
						${${$splicingSiteDiffHshByClusterHsh{$clusterName}}{"full"}}{$fullDiff}++;
						if ($refIsMajor eq "yes") {
							${$splicingSiteDiffMajorHsh{"full"}}{$fullDiff}++;
						}

						if ($clusterStrd eq "-") {#---for minus strand
		
							${$splicingSiteDiffAllHsh{"3"}}{$startDiff}++;
							${$splicingSiteDiffAllHsh{"5"}}{$endDiff}++;
							${${$splicingSiteDiffHshByClusterHsh{$clusterName}}{"3"}}{$startDiff}++;
							${${$splicingSiteDiffHshByClusterHsh{$clusterName}}{"5"}}{$endDiff}++;
							
							if ($refIsMajor eq "yes") {
								${$splicingSiteDiffMajorHsh{"3"}}{$startDiff}++;
								${$splicingSiteDiffMajorHsh{"5"}}{$endDiff}++;							
							}
							
						} else {#---for both no strand and plus strand

							${$splicingSiteDiffAllHsh{"5"}}{$startDiff}++;
							${$splicingSiteDiffAllHsh{"3"}}{$endDiff}++;
							${${$splicingSiteDiffHshByClusterHsh{$clusterName}}{"5"}}{$startDiff}++;
							${${$splicingSiteDiffHshByClusterHsh{$clusterName}}{"3"}}{$endDiff}++;
							
							if ($refIsMajor eq "yes") {
								${$splicingSiteDiffMajorHsh{"5"}}{$startDiff}++;
								${$splicingSiteDiffMajorHsh{"3"}}{$endDiff}++;					
							}
						}
					}
				}
			}#---end of if (@{$clusterAllJunctHsh{$clusterName}} > 1)
		}#---end of foreach my $clusterName (sort {$a cmp $b} keys %clusterAllJunctHsh) {

		#---plot the bar chart for all
		my $plotCutoff = 30;
		foreach my $bound53 (keys %splicingSiteDiffAllHsh) {
			
			my (%tmpForPlotFullHsh, %tmpForPlotCutoffHsh);
			
			foreach my $diff (keys  %{$splicingSiteDiffAllHsh{$bound53}}) {
				
				$tmpForPlotFullHsh{$diff} = ${$splicingSiteDiffAllHsh{$bound53}}{$diff};
				
				if (($diff <= $plotCutoff) and ($diff > 0)) {
					$tmpForPlotCutoffHsh{$diff} = ${$splicingSiteDiffAllHsh{$bound53}}{$diff};
				} elsif ($diff > $plotCutoff) {
					$tmpForPlotCutoffHsh{$plotCutoff} = 0 if (not exists $tmpForPlotCutoffHsh{$plotCutoff});
					$tmpForPlotCutoffHsh{$plotCutoff} += ${$splicingSiteDiffAllHsh{$bound53}}{$diff};
				}
			}
			
			GNUPlotBarChartNumberItem(\%tmpForPlotFullHsh, "intron $bound53 end altSplicingSiteDistance", "$outDir/plotData/intron$bound53.altSplicingSiteDistance.full.dat", "$outDir/plotPdf/intron$bound53.altSplicingSiteDistance.full.pdf", "Distance (nt)", "Frequency", "no", "no" );
			GNUPlotBarChartNumberItem(\%tmpForPlotCutoffHsh, "intron $bound53 end altSplicingSiteDistance", "$outDir/plotData/intron$bound53.altSplicingSiteDistance.cutoff.dat", "$outDir/plotPdf/intron$bound53.altSplicingSiteDistance.cutoff.pdf", "Distance (nt)", "Frequency", "no", "no");
		}

		#---plot the bar chart for major
		foreach my $bound53 (keys %splicingSiteDiffMajorHsh) {
			
			my (%tmpForPlotFullHsh, %tmpForPlotCutoffHsh);
			
			foreach my $diff (keys  %{$splicingSiteDiffMajorHsh{$bound53}}) {
				
				$tmpForPlotFullHsh{$diff} = ${$splicingSiteDiffMajorHsh{$bound53}}{$diff};
				
				if (($diff <= $plotCutoff) and ($diff > 0)) {
					$tmpForPlotCutoffHsh{$diff} = ${$splicingSiteDiffMajorHsh{$bound53}}{$diff};
				} elsif ($diff > $plotCutoff) {
					$tmpForPlotCutoffHsh{$plotCutoff} = 0 if (not exists $tmpForPlotCutoffHsh{$plotCutoff});
					$tmpForPlotCutoffHsh{$plotCutoff} += ${$splicingSiteDiffMajorHsh{$bound53}}{$diff};
				}
			}
			
			GNUPlotBarChartNumberItem(\%tmpForPlotFullHsh, "major intron $bound53 end altSplicingSiteDistance", "$outDir/plotData/intron$bound53.major.altSplicingSiteDistance.full.dat", "$outDir/plotPdf/intron$bound53.major.altSplicingSiteDistance.full.pdf", "Distance (nt)", "Frequency", "no", "no");
			GNUPlotBarChartNumberItem(\%tmpForPlotCutoffHsh, "major intron $bound53 end altSplicingSiteDistance", "$outDir/plotData/intron$bound53.major.altSplicingSiteDistance.cutoff.dat", "$outDir/plotPdf/intron$bound53.major.altSplicingSiteDistance.cutoff.pdf", "Distance (nt)", "Frequency", "no", "no");
		}

		return (\%splicingSiteDiffAllHsh, \%splicingSiteDiffHshByClusterHsh);
	
	}
	############################################## defineJunctionsOfProminentIsoform ########################################################
	sub defineJunctionsOfProminentIsoform {#---junctions on a prominent isoform is defined as a junction that has the highest number of read among all overlapping junctions

		#---rationale:
		#---check for each superjunction, to check whether it has the highest read num among all overlapping clusters.
		#---if yes, the so_ and sc_cluster (i.e. superjunction) will be stored in the prominent cluster hash
		#---if no, the overlapping clusters will be stored in the prominent cluster hash
		#---then go through each cluster and pick up the junctions with the highest number of reads (i.e. prominent read)
		#---go through each prominent read to see if they overlap with each other, it may happen for two s_ clusters since superjunction clusters 
		
		print "Defining the junctions of the most prominent isoform.\n";
		my (%junctOnProminentIsofmHsh);
		
		my %junctClusterNameByJunctHsh = %{$_[0]};
		my %junctClusterInfoHsh = %{$_[1]};
		my %junctClusterSSOverlapHsh = %{$_[2]};
		
		#----scan the superJunctionFirst;
		foreach my $clusterName (keys %junctClusterInfoHsh) {
			if ($clusterName =~ m/^s/) {#---superJunction
				my $superClusterName = $clusterName;
				my $superSupportReadNum = ${$junctClusterInfoHsh{$superClusterName}}{"mostSupportReadNum"};
				my $superIsProminent = "yes";
				
				#----check the inferior clusters
				foreach my $overlapInferiorCluster (keys %{$junctClusterSSOverlapHsh{$superClusterName}}) {
					my $inferiorSupportReadNum = ${$junctClusterInfoHsh{$overlapInferiorCluster}}{"mostSupportReadNum"};
					if ($inferiorSupportReadNum > $superSupportReadNum) {
						$superIsProminent = "no";
					}
				}
				
				if ($superIsProminent eq "yes") {
					#---skip the inferiorClsuter if the superJunction cluster is prominent;
					foreach my $overlapInferiorCluster (keys %{$junctClusterSSOverlapHsh{$superClusterName}}) {
						${$junctClusterInfoHsh{$overlapInferiorCluster}}{"onProminentIsoform"} = "no";#---to skip the inferior junct later on since the super junction is major
					}
					
					my $clusterIsMajor = ${$junctClusterInfoHsh{$superClusterName}}{"majorCluster"};
					
					if ($clusterIsMajor eq "yes") {
						my $prominentJunctStr = ${$junctClusterInfoHsh{$superClusterName}}{"prominentJunct"};
						${$junctClusterInfoHsh{$superClusterName}}{"onProminentIsoform"} = "yes";
						$junctOnProminentIsofmHsh{$prominentJunctStr} = $superSupportReadNum;
						
					} else {#---superJunction has the major support num but not major
						${$junctClusterInfoHsh{$superClusterName}}{"onProminentIsoform"} = "no";
					}
					
				} else { #---superJunction is not prominent
					${$junctClusterInfoHsh{$superClusterName}}{"onProminentIsoform"} = "no";
				}
			}
		} #---end of foreach my $clusterName (keys %junctClusterInfoHsh)
		
		#----scan the unprocessed inferior junctions i.e. not exists ${$junctClusterInfoHsh{$clusterName}}{"onProminentIsoform"}
		foreach my $clusterName (keys %junctClusterInfoHsh) {
			if (not exists ${$junctClusterInfoHsh{$clusterName}}{"onProminentIsoform"}) {
				my $clusterIsMajor = ${$junctClusterInfoHsh{$clusterName}}{"majorCluster"};
				if ($clusterIsMajor eq "yes") {
					my $prominentJunctStr = ${$junctClusterInfoHsh{$clusterName}}{"prominentJunct"};
					${$junctClusterInfoHsh{$clusterName}}{"onProminentIsoform"} = "yes";
					my $supportReadNum = ${$junctClusterInfoHsh{$clusterName}}{"mostSupportReadNum"};
					$junctOnProminentIsofmHsh{$prominentJunctStr} = $supportReadNum;
				} else {
					${$junctClusterInfoHsh{$clusterName}}{"onProminentIsoform"} = "no";
				}
			}
		}
		
		#---Just double check to make sure the junctions on the prominent isofrom doesnt overlapp with each other
		print "Making sure the junctions on the prominent isofrom do not overlap with each other.\n";
		my %tmpCheckJunctHsh;
		foreach my $refJunctStr (keys %junctOnProminentIsofmHsh) {
			$tmpCheckJunctHsh{$refJunctStr}++;
			foreach my $qryJunctStr (keys %junctOnProminentIsofmHsh) {
				next if (exists $tmpCheckJunctHsh{$qryJunctStr});
				if ($refJunctStr ne $qryJunctStr) {
					die "Unexpected scenerio: Junctions on prominent isoform overlap with each other.\n" if (exists ${$SSOvrlpJunctIntronHsh{$refJunctStr}}{$qryJunctStr});
				}
			}
		}
		%tmpCheckJunctHsh = ();
		
		#---print the junction on prominent isoform
		open (TMPLOG, ">$outDir/log/junctOnProminentIsofm.txt");
		foreach my $junctStr (sort{$junctClusterNameByJunctHsh{$a} cmp $junctClusterNameByJunctHsh{$b}} keys %junctOnProminentIsofmHsh) {
			print TMPLOG $junctStr."\t".$junctClusterNameByJunctHsh{$junctStr}."\n";
		}
		close TMPLOG;
		
		return (\%junctClusterInfoHsh, \%junctOnProminentIsofmHsh);
		
	}
	
}
########################################################################## buildGeneModel
sub buildGeneModel {
	
	my %clusterJunctTrnsfgHsh = %{$_[0]};
	my %junctOnProminentIsofmHsh = %{$_[1]};
	my %junctClusterInfoHsh = %{$_[2]};
	my %junctClusterSSOverlapHsh = %{$_[3]};
	my %trnsfgRngHsh = %{$_[4]};
	my %junctRngHsh = %{$_[5]};
	my %trnsfgCntgHsh = %{$_[6]};
	my %trnsfgStrdHsh = %{$_[7]};
	my %junctClusterNameByJunctHsh = %{$_[8]};
	
	#---build the prominent isoform from clusters
	print "Bulding the prominent isoforms.\n";
	my (%allTrsncptIsofmBoundHsh, %trnsfgInJTClusterHsh, %geneStrandHsh, %geneCntgHsh);
	foreach my $JTClusterName (sort  {$a cmp $b} keys %clusterJunctTrnsfgHsh) {

		#-- create an empty array, essential in cases of the cluster contains only minor jucntions
		@{${${$allTrsncptIsofmBoundHsh{$JTClusterName}}{"0"}}{"intronBound"}} = ();
		@{${${$allTrsncptIsofmBoundHsh{$JTClusterName}}{"0"}}{"fragBound"}} = ();

		my $cntg = ${$clusterJunctTrnsfgHsh{$JTClusterName}}{"cntg"};
		my %tmpJunctClusterOnJTClusterHsh; 
		
		#---get the bound of the prominent junctions
		foreach my $junctStr (@{${$clusterJunctTrnsfgHsh{$JTClusterName}}{"junct"}}) {

			my $junctCluster = $junctClusterNameByJunctHsh{$junctStr};

			if (exists $junctOnProminentIsofmHsh{$junctStr}) {#---prominent junctions
				${$tmpJunctClusterOnJTClusterHsh{"prominent"}}{$junctCluster} = $junctStr; #---store the prominent junction of the prominent cluster
				my @junctStrSplt = split /:/, $junctStr;
				my $intronStart = $junctStrSplt[1] - 1;
				my $intronEnd = $junctStrSplt[-1] + 1;
				push @{${${$allTrsncptIsofmBoundHsh{$JTClusterName}}{"0"}}{"intronBound"}}, $intronStart;
				push @{${${$allTrsncptIsofmBoundHsh{$JTClusterName}}{"0"}}{"intronBound"}}, $intronEnd;
				my $bedStart = ${${$junctRngHsh{$cntg}}{$junctStr}}{"start"};
				my $bedEnd = ${${$junctRngHsh{$cntg}}{$junctStr}}{"end"};
				push @{${${$allTrsncptIsofmBoundHsh{$JTClusterName}}{"0"}}{"fragBound"}}, $bedStart;
				push @{${${$allTrsncptIsofmBoundHsh{$JTClusterName}}{"0"}}{"fragBound"}}, $bedEnd;
			}
		}
		
		#---record the minor clusters
		foreach my $junctStr (@{${$clusterJunctTrnsfgHsh{$JTClusterName}}{"junct"}}) {
			my $junctCluster = $junctClusterNameByJunctHsh{$junctStr};
			
			if (not exists ${$tmpJunctClusterOnJTClusterHsh{"prominent"}}{$junctCluster}) {
				${$tmpJunctClusterOnJTClusterHsh{"minor"}}{$junctCluster} = ${$junctClusterInfoHsh{$junctCluster}}{"prominentJunct"}; #---store the prominent junction of the minor cluster
			}
		}

		
		#---consider the minor junction clusters
		my $minorIsoformNum = 0;
		foreach my $minorJunctClusters (sort {$a cmp $b} keys %{$tmpJunctClusterOnJTClusterHsh{"minor"}}) {
			$minorIsoformNum++;
			my $minorJunctStr = ${$tmpJunctClusterOnJTClusterHsh{"minor"}}{$minorJunctClusters};
			my @junctStrSplt = split /:/, $minorJunctStr;
			my $intronStart = $junctStrSplt[1] - 1;
			my $intronEnd = $junctStrSplt[-1] + 1;
			push @{${${$allTrsncptIsofmBoundHsh{$JTClusterName}}{$minorIsoformNum}}{"intronBound"}}, $intronStart;
			push @{${${$allTrsncptIsofmBoundHsh{$JTClusterName}}{$minorIsoformNum}}{"intronBound"}}, $intronEnd;
			my $bedStart = ${${$junctRngHsh{$cntg}}{$minorJunctStr}}{"start"};
			my $bedEnd = ${${$junctRngHsh{$cntg}}{$minorJunctStr}}{"end"};
			push @{${${$allTrsncptIsofmBoundHsh{$JTClusterName}}{$minorIsoformNum}}{"fragBound"}}, $bedStart;
			push @{${${$allTrsncptIsofmBoundHsh{$JTClusterName}}{$minorIsoformNum}}{"fragBound"}}, $bedEnd;


			#---check if it overlaps with the prominent clusters (in cases of superjunction is on the prominent isoform)
			foreach my $prominentJunctCluster (keys %{$tmpJunctClusterOnJTClusterHsh{"prominent"}}) {
				if (not exists ${$junctClusterSSOverlapHsh{$minorJunctClusters}}{$prominentJunctCluster}) {#---this minor cluster is not overlapping with the prominet cluster
					my $prominentJunctStr = ${$tmpJunctClusterOnJTClusterHsh{"prominent"}}{$prominentJunctCluster};
					my @junctStrSplt = split /:/, $prominentJunctStr;
					my $intronStart = $junctStrSplt[1] - 1;
					my $intronEnd = $junctStrSplt[-1] + 1;
					push @{${${$allTrsncptIsofmBoundHsh{$JTClusterName}}{$minorIsoformNum}}{"intronBound"}}, $intronStart;
					push @{${${$allTrsncptIsofmBoundHsh{$JTClusterName}}{$minorIsoformNum}}{"intronBound"}}, $intronEnd;
					my $bedStart = ${${$junctRngHsh{$cntg}}{$prominentJunctStr}}{"start"};
					my $bedEnd = ${${$junctRngHsh{$cntg}}{$prominentJunctStr}}{"end"};
					push @{${${$allTrsncptIsofmBoundHsh{$JTClusterName}}{$minorIsoformNum}}{"fragBound"}}, $bedStart;
					push @{${${$allTrsncptIsofmBoundHsh{$JTClusterName}}{$minorIsoformNum}}{"fragBound"}}, $bedEnd;
				}
			}
		}
		
		#---add the transfrag bounds to all isoform and reconstruct the bounds
		foreach my $isoformNum (keys %{$allTrsncptIsofmBoundHsh{$JTClusterName}}) {

			#---add the trnsfg bounds
			foreach my $trnsfg (@{${$clusterJunctTrnsfgHsh{$JTClusterName}}{"trnsfg"}}) {
				my $trnfgStart = ${${$trnsfgRngHsh{$cntg}}{$trnsfg}}{"start"};
				my $trnfgEnd = ${${$trnsfgRngHsh{$cntg}}{$trnsfg}}{"end"};
				push @{${${$allTrsncptIsofmBoundHsh{$JTClusterName}}{$isoformNum}}{"fragBound"}}, $trnfgStart;		
				push @{${${$allTrsncptIsofmBoundHsh{$JTClusterName}}{$isoformNum}}{"fragBound"}}, $trnfgEnd;
				$trnsfgInJTClusterHsh{$trnsfg}++;
			}

			#---reconstruct the bounds
			my @fragBoundSortedAry = sort {$a <=> $b} @{${${$allTrsncptIsofmBoundHsh{$JTClusterName}}{$isoformNum}}{"fragBound"}};
			my $leftMostBound = $fragBoundSortedAry[0];
			my $rightMostBound = $fragBoundSortedAry[-1];
			my @allIntronBoundAry = @{${${$allTrsncptIsofmBoundHsh{$JTClusterName}}{$isoformNum}}{"intronBound"}};
			@{${${$allTrsncptIsofmBoundHsh{$JTClusterName}}{$isoformNum}}{"combineBound"}} = sort {$a <=> $b} (@allIntronBoundAry, $rightMostBound, $leftMostBound);
			$geneStrandHsh{$JTClusterName} = ${$clusterJunctTrnsfgHsh{$JTClusterName}}{"strand"};
			$geneCntgHsh{$JTClusterName} = $cntg;
			delete ${${$allTrsncptIsofmBoundHsh{$JTClusterName}}{$isoformNum}}{"intronBound"};	
			delete ${${$allTrsncptIsofmBoundHsh{$JTClusterName}}{$isoformNum}}{"fragBound"};	
		}	

	}
	
	#---add the transfrag without junctions into the allTrsncptIsofmBoundHsh;
	foreach my $trnsfg (sort {$trnsfgCntgHsh{$a} cmp $trnsfgCntgHsh{$b}} keys %trnsfgCntgHsh) {
		if (not exists $trnsfgInJTClusterHsh{$trnsfg}) {
			my $cntg = $trnsfgCntgHsh{$trnsfg};
			my $trnfgStart = ${${$trnsfgRngHsh{$cntg}}{$trnsfg}}{"start"};
			my $trnfgEnd = ${${$trnsfgRngHsh{$cntg}}{$trnsfg}}{"end"};
			push @{${${$allTrsncptIsofmBoundHsh{$trnsfg}}{"0"}}{"combineBound"}}, $trnfgStart;
			push @{${${$allTrsncptIsofmBoundHsh{$trnsfg}}{"0"}}{"combineBound"}}, $trnfgEnd;
			$geneStrandHsh{$trnsfg} = $trnsfgStrdHsh{$trnsfg};
			$geneCntgHsh{$trnsfg} = $cntg;
		}
	}
	
	#---print the all, prominent and minor isoforms to different GTF files
	open (AGTFPATH, ">$outDir/finalGTF/allIsoform.gtf");
	open (PGTFPATH, ">$outDir/finalGTF/prominentIsoform.gtf");
	open (MGTFPATH, ">$outDir/finalGTF/minorIsoform.gtf");
	open (AGFFPATH, ">$outDir/finalGTF/allIsoform.gff");
	open (PGFFPATH, ">$outDir/finalGTF/prominentIsoform.gff");
	open (MGFFPATH, ">$outDir/finalGTF/minorIsoform.gff");

	foreach my $gene (sort {$a cmp $b} keys %allTrsncptIsofmBoundHsh) {
		foreach my $trnscptNum (sort {$a <=> $b} keys %{$allTrsncptIsofmBoundHsh{$gene}}) {
			my $trnsptStart = ${${${$allTrsncptIsofmBoundHsh{$gene}}{$trnscptNum}}{"combineBound"}}[0];
			my $trnsptEnd = ${${${$allTrsncptIsofmBoundHsh{$gene}}{$trnscptNum}}{"combineBound"}}[-1];
			my $cntg = $geneCntgHsh{$gene};
			my $strand = $geneStrandHsh{$gene};
			my $boundNum = @{${${$allTrsncptIsofmBoundHsh{$gene}}{$trnscptNum}}{"combineBound"}};
		
			#---print the transcript range first
			print AGTFPATH join "\t", ($cntg, "geneModelBuilder", "transcript", $trnsptStart, $trnsptEnd, "1000", $strand, ".", "gene_id \"$gene\"; transcript_id \"$gene.$trnscptNum\";\n");
			if ($trnscptNum == 0) {
				print PGTFPATH join "\t", ($cntg, "geneModelBuilder", "transcript", $trnsptStart, $trnsptEnd, "1000", $strand, ".", "gene_id \"$gene\"; transcript_id \"$gene.$trnscptNum\";\n");
			} else {
				print MGTFPATH join "\t", ($cntg, "geneModelBuilder", "transcript", $trnsptStart, $trnsptEnd, "1000", $strand, ".", "gene_id \"$gene\"; transcript_id \"$gene.$trnscptNum\";\n");
			}

			#Gff
			#DS572165	geneModelReviser	gene	897	1582	.	+	.	ID=trnsfrg_12610_0
			#DS572165	geneModelReviser	transfrag	897	1582	.	+	.	ID=rna_trnsfrg_12610_0-1;Parent=trnsfrg_12610_0
			#DS572165	geneModelReviser	exon	897	1582	.	+	.	ID=exon_trnsfrg_12610_0-1;Parent=rna_trnsfrg_12610_0-1

			print AGFFPATH join "\t", ($cntg, "geneModelBuilder", "gene", $trnsptStart, $trnsptEnd, "1000", $strand, ".", "ID=$gene\n");
			print AGFFPATH join "\t", ($cntg, "geneModelBuilder", "transcript", $trnsptStart, $trnsptEnd, "1000", $strand, ".", "ID=rna_$gene\-1;Parent=$gene\n");
			if ($trnscptNum == 0) {
				print PGFFPATH join "\t", ($cntg, "geneModelBuilder", "gene", $trnsptStart, $trnsptEnd, "1000", $strand, ".", "ID=$gene\n");
				print PGFFPATH join "\t", ($cntg, "geneModelBuilder", "transcript", $trnsptStart, $trnsptEnd, "1000", $strand, ".", "ID=rna_$gene\-1;Parent=$gene\n");
			} else {
				print MGFFPATH join "\t", ($cntg, "geneModelBuilder", "gene", $trnsptStart, $trnsptEnd, "1000", $strand, ".", "ID=$gene\n");
				print MGFFPATH join "\t", ($cntg, "geneModelBuilder", "transcript", $trnsptStart, $trnsptEnd, "1000", $strand, ".", "ID=rna_$gene\-1;Parent=$gene\n");
			}

			#---print all the exon ranges
			my $exonNum = 0;
			for (my $i=0; $i < $boundNum; $i = $i + 2) {
				$exonNum++;
				my $exonStart = ${${${$allTrsncptIsofmBoundHsh{$gene}}{$trnscptNum}}{"combineBound"}}[$i];
				my $exonEnd = ${${${$allTrsncptIsofmBoundHsh{$gene}}{$trnscptNum}}{"combineBound"}}[$i+1];
				print AGTFPATH join "\t", ($cntg, "geneModelBuilder", "exon", $exonStart, $exonEnd, "1000", $strand, ".", "gene_id \"$gene\"; transcript_id \"$gene.$trnscptNum\"; exon_number \"$exonNum\";\n");
				print AGFFPATH join "\t", ($cntg, "geneModelBuilder", "exon", $exonStart, $exonEnd, "1000", $strand, ".", "ID=exon_$gene\-$exonNum;Parent=rna_$gene\-1\n");
				if ($trnscptNum == 0) {
					print PGTFPATH join "\t", ($cntg, "geneModelBuilder", "exon", $exonStart, $exonEnd, "1000", $strand, ".", "gene_id \"$gene\"; transcript_id \"$gene.$trnscptNum\"; exon_number \"$exonNum\";\n");
					print PGFFPATH join "\t", ($cntg, "geneModelBuilder", "exon", $exonStart, $exonEnd, "1000", $strand, ".", "ID=exon_$gene\-$exonNum;Parent=rna_$gene\-1\n");
				} else {
					print MGTFPATH join "\t", ($cntg, "geneModelBuilder", "exon", $exonStart, $exonEnd, "1000", $strand, ".", "gene_id \"$gene\"; transcript_id \"$gene.$trnscptNum\"; exon_number \"$exonNum\";\n");
					print MGFFPATH join "\t", ($cntg, "geneModelBuilder", "exon", $exonStart, $exonEnd, "1000", $strand, ".", "ID=exon_$gene\-$exonNum;Parent=rna_$gene\-1\n");
				}
			}
		}
	}
	close AGTFPATH;
	close PGTFPATH;
	close MGTFPATH;	
	close AGFFPATH;
	close PGFFPATH;
	close MGFFPATH;
}
########################################################################## GNUPlotXYScatterWithLines
sub GNUPlotXYScatterWithLines {

	my %XYHsh = %{$_[0]};
	my $plotFilePath = $_[1];
	my $plotDataPath = $_[2];
	my $xlable = $_[3];
	my $ylable = $_[4];
	my $xscale = $_[5];
	my $yscale = $_[6];
	my $title = $_[7];
	my $verbose = $_[8];
	
	my $GNULogXCmd = "";
	$GNULogXCmd = "set logscale x" if ($xscale eq "log");
	my $GNULogYCmd = "";
	$GNULogYCmd = "set logscale y" if ($yscale eq "log");

	$plotFilePath .= ".pdf" if ($plotFilePath !~ m/\.pdf$/);

	my @filePathSplt = split /\//, $plotFilePath;
	my $fileName = $filePathSplt[-1];

	print "Running GNUPlotXYScatterWithLines for $fileName.\n" if ($verbose eq "yes");
	
	#---creat a tmp file
	open (TMPFILE, ">$plotDataPath");
	for my $x (sort {$a <=> $b} keys %XYHsh) {
		print TMPFILE $x."\t".$XYHsh{$x}."\n";
	}
	close TMPFILE;
	
	#---do the GNUPLOT
	open (GNUPLOT, "|gnuplot");
	print GNUPLOT <<EOPLOT;
	set terminal postscript color solid
	set output "| ps2pdf - $plotFilePath";
	unset logscale x; 
	unset logscale y; 
	$GNULogXCmd;
	$GNULogYCmd;
	set xlabel "$xlable";
	set ylabel "$ylable";
	set title "$title";
	set nokey;
   	plot '$plotDataPath' using 1:2 with lines;
EOPLOT
	close(GNUPLOT);
	#rmtree(['tmp.dat'], 0, 1); #---non-verbose removal of tmp file
}
########################################################################## GNUPLOTAryHistogram
sub GNUPLOTAryHistogram {

	my @numAry = @{$_[0]};
	my $plotFilePath = $_[1];
	my $plotDataPath = $_[2];
	my $xlable = $_[3];
	my $ylable = $_[4];
	my $title = $_[5];
	my $verbose = $_[6];
	
	my @filePathSplt = split /\//, $plotFilePath;
	my $fileName = $filePathSplt[-1];
	
	print "Running GNUPLOTAryHistogram for $fileName.\n" if ($verbose eq "yes");
	
	#---calculate the optimal bin
	@numAry = sort {$a <=> $b} @numAry;
	my $tail5PctPos = int ($#numAry*0.05);
	my $upper95PctPos = $#numAry - $tail5PctPos;
	my $lower95PctPos = $tail5PctPos;
	my $binWidth = ($numAry[$upper95PctPos] - $numAry[$lower95PctPos])/100;
	#$binWidth = 1 if ($binWidth < 1);

	#---creat a tmp file
	open (TMPFILE, ">$plotDataPath");
	for (@numAry) {print TMPFILE $_."\n";}
	close TMPFILE;
	
	#---do the GNUPLOT
	open (GNUPLOT, "|gnuplot");
	print GNUPLOT <<EOPLOT;
	set terminal postscript color solid
	bin_width = $binWidth
	bin_number(x) = floor(x/bin_width)
	rounded(x) = bin_width * (bin_number(x) + 0.5)
	UNITY = 1
	unset logscale x;
	unset logscale y;
	set output "| ps2pdf - $plotFilePath";
	set xlabel "$xlable";
	set ylabel "$ylable";
	set title "$title";
   	plot '$plotDataPath' u (rounded(\$1)):(UNITY) t 'bin at $binWidth' smooth frequency w histeps
EOPLOT
	close(GNUPLOT);
	#rmtree(['./tmpRand.dat'], 0, 1); #---non-verbose removal of tmp file
}
########################################################################## generateCumulativePctHash
sub generateCumulativePctHash {

	#---both the keys and the values have to numbers
	my %countHsh = %{$_[0]};
	my $order = $_[1]; #---a for ascending or d for descending

	#---find the total
	my $totalCount = 0;
	foreach my $item (sort {$a <=> $b} keys %countHsh) {$totalCount += $countHsh{$item};}
	
	#---generate the cumulative values;
	my $cumulativeCount = 0;
	my %cumulativeHsh;
	if ($order eq "a") {
		foreach my $item (sort {$a <=> $b} keys %countHsh) {
			$cumulativeCount += $countHsh{$item};
			my $cumulativePct = sprintf "%.02f", ($cumulativeCount/$totalCount)*100;
			$cumulativeHsh{$item} = $cumulativePct;
		}

	} elsif ($order eq "d") {
		foreach my $item (sort {$b <=> $a} keys %countHsh) {
			$cumulativeCount += $countHsh{$item};
			my $cumulativePct = sprintf "%.02f", ($cumulativeCount/$totalCount)*100;
			$cumulativeHsh{$item} = $cumulativePct;
		}
	}	
	
	return \%cumulativeHsh;
}
########################################################################## GNUPlotBarChartNumberItem
sub GNUPlotBarChartNumberItem {

	#---assuming the keys of the hash are numbers
	my %barDataHsh = %{$_[0]};
	my $plotTitle = $_[1];
	my $dataFilePath = $_[2];
	my $plotFilePath = $_[3];
	my $xLabel = $_[4];
	my $yLabel = $_[5];
	my $removeDataFile = $_[6];
	my $verbose = $_[7];
	
	#---add extension to =plot filename if it doesnt have
	$plotFilePath .= ".pdf" if ($plotFilePath !~ m/\.pdf$/);
	
	my @filePathSplt = split /\//, $plotFilePath;
	my $fileName = $filePathSplt[-1];

	print "Running GNUPlotBarChartNumberItem for $fileName.\n" if ($verbose eq "yes");

	#---get the total count
	my $totalCount = 0;
	foreach my $item (sort {$a <=> $b} keys %barDataHsh) {$totalCount += $barDataHsh{$item};}

	#---generate a tmpData file
	my $tmpBarDatPath = $dataFilePath;
	open (TMPBARDATA, ">$tmpBarDatPath");
	my $maxPct = 0;
	foreach my $item (sort {$a <=> $b} keys %barDataHsh) {
		my $itemCount = $barDataHsh{$item};
		my $tmpPct = sprintf "%.2f", ($itemCount/$totalCount)*100;
		print TMPBARDATA $item."\t".$tmpPct."\t".$itemCount."\n";
		$maxPct = $tmpPct if ($tmpPct > $maxPct);
	}
	close TMPBARDATA;
	
	#---define the max y values
	my $yMax = ((int ($maxPct/10))*10) + 10; #---in scale of 10%
	my $ytics = int ($yMax/10);
	my $y2Max = int ($totalCount*($yMax/100));
	my $y2tics = int ($y2Max/10);

	open (GNUPLOT, "|gnuplot");
	print GNUPLOT <<EOPLOT;
	set terminal postscript color solid
	set output "| ps2pdf - $plotFilePath";
	set key left;
	set title "$plotTitle";
	set xlabel "$xLabel";
	set ylabel "% $yLabel";
	set y2label "count $yLabel";
	set yrange [0:$yMax];
	set ytics $ytics nomirror;
	set y2range [0:$y2Max];
	set y2tics $y2tics nomirror;
	set boxwidth 1 absolute;
	set style fill solid 2.00 border 0;
	set style histogram;
	set style data histograms;
	plot '$tmpBarDatPath' using 2:xtic(1);
EOPLOT
	
	system ("rm $tmpBarDatPath") if $removeDataFile eq "yes";
	
}
