#!/usr/bin/perl
#
# This script parses a HC VCF output run on a normal/tumor pair.  The 10th and 11th column of
# the VCF should be variant calls from either normal or tumor.  The first argument is the name
# of the VCF file to be parsed. Second argument is sample name of the normal sample.  Third argument
# is the sample name of the tumor sample.  
#
# parseMergeVCF.pl normal-tumor.hc.snpEff.vcf normalSampleName tumorSampleName
#
# eg. parseMergeVCF.pl MMRF_1355_2_PB_Whole_C1_KAS5U_L01938-MMRF_1355_1_BM_CD138pos_T1_KAS5U_L01927.UG.vcf MMRF_1355_2_PB_Whole_C1_KAS5U_L01938 MMRF_1355_1_BM_CD138pos_T1_KAS5U_L01927
#
#Jessica Aldrich
#TGen
#Dec 30, 2012
#Modified: May 12, 2015

use List::Util qw( sum min );
use List::MoreUtils 'first_index';

open(VCF,"$ARGV[0]");
open(OFILE,">","merged.vcf.txt");
#open(OFILE2,">","baf.txt");

print OFILE "Chromosome\tPosition\tNormalReadDepth\tNormalRefAllele\tNormalAltAllele\tTumorReadDepth\tTumorRefAllele\tTumorAltAllele\n";
#print OFILE2 "Chromosome\tPosition\tNormalReadDepth\tNormalRefAllele\tNormalAltAllele\tTumorReadDepth\tTumorRefAllele\tTumorAltAllele\tCONTROLBAF\tCONTROLMAF\tBAF\tMAF\tABSBAFDEV\n";

$normalName=$ARGV[1];
$tumorName=$ARGV[2];

LOOP:while ($line=<VCF>){

	chomp($line);

	if ($line=~/^##/){next LOOP};
	if ($line=~/^#CHROM/){	
		@tmp=split(/\t+/,$line);
		$szTMP=@tmp;
                for ( my $i=9; $i <= $szTMP; $i++){
                if ($tmp[$i] eq $normalName){
                        $norm=$i;
                }
                if ($tmp[$i] eq $tumorName){
                        $tumor=$i;
                }
                }

	}
	
	@temp=split(/\t+/,$line);

	if ($temp[9] eq '.' || $temp[10] eq '.'){next LOOP};
	
	@COL8=split(/:/,$temp[8]);	
	
	$GTind = first_index { /GT/ } @COL8;
	$ADind = first_index { /AD/ } @COL8;
	$DPind = first_index { /DP/ } @COL8;
	if ($DPind == -1){next LOOP};	

	@control = split(/:/,$temp[$norm]);
	@affected = split(/:/,$temp[$tumor]);
	
	@con_allele=split(/,/,$control[$ADind]);
	@aff_allele=split(/,/,$affected[$ADind]);	

	if ( $temp[7] =~ /GMAF=/ || $temp[7] =~ /CAF=/ ){
	
		if ($temp[7] =~ /GMAF=/) {
			$temp[7] =~ /GMAF=(.*?);/;
			$GMAF = $1;
		}

		if ($temp[7] =~ /CAF=/ ){
			$temp[7] =~ /CAF=(.*?);/;
			@CAF = split(/,/,$1);
			$GMAF = min(@CAF);
		}

		$charR = length($temp[3]);
		$charA = length($temp[4]);
	
		# Write out baf.txt
		#if ( $GMAF > 0.05 && $control[0] eq '0/1' && $charR < 2 && $charA < 2 && sum(@con_allele) > 50 && sum(@aff_allele) > 50 && $con_allele[1]/sum(@con_allele) > 0.05) {
	
		#	$baf = $aff_allele[1]/sum(@aff_allele);
		#	$maf = min(@aff_allele)/sum(@aff_allele);
		#	$control_baf = $con_allele[1]/sum(@con_allele);
		#	$control_maf = min(@con_allele)/sum(@con_allele);
		#	$bafdev = abs(0.5-$baf);		

		#	if($temp[0] eq "X"){$temp[0]=23};
        	#	if($temp[0] eq "Y"){$temp[0]=24};

                #	print OFILE2 "$temp[0]\t$temp[1]\t$control[$DPind]\t$con_allele[0]\t$con_allele[1]\t$affected[$DPind]\t$aff_allele[0]\t$aff_allele[1]\t$control_baf\t$control_maf\t$baf\t$maf\t$bafdev\n";

        	#}		
		
		# Write out merged.vcf.txt
		if ( $temp[2] =~ /^rs.*/ && $control[0] eq '0/1' && sum(@con_allele) > 10 && sum(@aff_allele) > 10) {
		
		        print OFILE "$temp[0]\t$temp[1]\t$control[$DPind]\t$con_allele[0]\t$con_allele[1]\t$affected[$DPind]\t$aff_allele[0]\t$aff_allele[1]\n";
	
	        }
	}
}

close(VCF);
close(OFILE);
#close(OFILE2);

