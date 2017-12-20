use strict;
use REST::Client;
use JSON;
use Data::Dumper;
use feature qw(switch);
use clsEpiRR;



my $url = $ARGV[0];
my $look_for_reference = $ARGV[1];

if (not defined $url) {
	die "the 'url' parameter must be indicated\n";
}
if (not defined $look_for_reference) {
	$look_for_reference = '';
}

$url.='&frame=object&format=json';

my $clientSum = REST::Client->new();
my @arrClsEpi; #array to store the clsEpiRR objects

#create columns headers
my $clsEpi = clsEpiRR->new();
$clsEpi->set_project('project');
$clsEpi->set_species('species');
$clsEpi->set_accession('accession');
$clsEpi->set_ihec('ihec');
$clsEpi->set_epigenome('epigenome');
$clsEpi->set_referenceStage('reference stage');
$clsEpi->set_ATAC_Seq ('ATAC-Seq');
$clsEpi->set_ChipS_In ('ChIP-Seq Input');
$clsEpi->set_DNase_Seq ('DNase-Seq');
$clsEpi->set_H3K27ac ('H3K27ac');
$clsEpi->set_H3K27me3 ('H3K27me3');
$clsEpi->set_H3K36me3 ('H3K36me3');
$clsEpi->set_H3K4me1 ('H3K4me1');
$clsEpi->set_H3K4me3 ('H3K4me3');
$clsEpi->set_H3K9me3 ('H3K9me3');
$clsEpi->set_CTCF ('CTCF');
$clsEpi->set_avExperiments ('Available Experiments');
$clsEpi->set_others ('other');

#add the object to the array
push @arrClsEpi, $clsEpi;


#Call the Api
my $headers = {Accept => 'application/json'};
$clientSum->GET($url, $headers);
my $responseSum = decode_json($clientSum->responseContent());
	
#count the results
my @num=@{$responseSum->{'@graph'}};
my $numTot = scalar @num;

#loop over series	
for (my $i=0; $i < $numTot; $i++){
	my $project ='ENCODE';
	my $species = $responseSum->{'@graph'}[$i]->{'organism'}[0];
	$species = (split '/', $species)[-1];
	my $accession = $responseSum->{'@graph'}[$i]->{'accession'};
	my $ihec = $responseSum->{'@graph'}[$i]->{'dbxrefs'}[0];
	$ihec = (split ':', $ihec)[-1];
	my $epigenome = $responseSum->{'@graph'}[$i]->{'biosample_term_name'}[0];
	my $refStage = '';
	my $atacSeq = 'NO';
	my $chipSeqInput = 'NO';
	my $DNaseSeq='NO';
	my $H3K27ac = 'NO';
	my $H3K27me3 = 'NO';
	my $H3K36me3 = 'NO';
	my $H3K4me1 = 'NO';
	my $H3K4me3 = 'NO';
	my $H3K9me3 = 'NO';
	my $CTCF = 'NO';
	my $others = '';
	my $availableExp = 0;
	
	#count assay terms
	my @lstAssayTerms = @{$responseSum->{'@graph'}[$i]->{'assay_term_name'}};
	my $numAssayTerms = scalar @lstAssayTerms;
	#look for DNase-Seq
	for (my $e=0; $e < $numAssayTerms; $e++){
		my $assayTerm = $responseSum->{'@graph'}[$i]->{'assay_term_name'}[$e];
		if (uc $assayTerm eq 'DNASE-SEQ') {
			$DNaseSeq = 'YES';
		}
		if (uc $assayTerm eq 'ATAC-SEQ') {
			$atacSeq = 'YES';
		}
	}
	
	
	if (uc $look_for_reference eq 'Y') {
		#check for the reference stage consistence
		$refStage = reference_stage($responseSum->{'@graph'}[$i]->{'related_datasets'});
	}
	#exit 0;
	#count targets
	my @lstTargets = @{$responseSum->{'@graph'}[$i]->{'target'}};
	my $numTargets = scalar @lstTargets;
	#find Histone modification and Chip-Seq Input and count available
	for (my $e=0; $e < $numTargets; $e++){
		my $rawTarget = (split '/', $responseSum->{'@graph'}[$i]->{'target'}[$e])[-1];
		my $target = (split '-', $rawTarget)[0];
		given ( uc $target){
			when ('H3K27AC') {
				$H3K27ac = 'YES';
				$availableExp ++;
			}
			when ('H3K27ME3') {
				$H3K27me3 = 'YES';
				$availableExp ++;
			}
			when ('H3K36ME3') {
				$H3K36me3 = 'YES';
				$availableExp ++;
			}
			when ('H3K4ME1') {
				$H3K4me1 = 'YES';
				$availableExp ++;
			}
			when ('H3K4ME3') {
				$H3K4me3 = 'YES';
				$availableExp ++;
			}
			when ('H3K9ME3') {
				$H3K9me3 = 'YES';
				$availableExp ++;
			}	
			when ('CTCF') {
				$CTCF = 'YES';
			}
			when ('CONTROL') {
				$chipSeqInput = 'YES';
				$availableExp ++;
			}
			default {
				if ($others ne '') {
					$others .= ', ';
				}
				$others .= $target;
			}
		}
		
	}
	
		#load the data in the clsEpiRR object
		$clsEpi = clsEpiRR->new();
		$clsEpi->set_project($project);
		$clsEpi->set_species($species);
		$clsEpi->set_accession($accession);
		$clsEpi->set_ihec($ihec);
		$clsEpi->set_epigenome($epigenome);
		$clsEpi->set_referenceStage($refStage);
		$clsEpi->set_ATAC_Seq ($atacSeq);
		$clsEpi->set_ChipS_In ($chipSeqInput);
		$clsEpi->set_DNase_Seq ($DNaseSeq);
		$clsEpi->set_H3K27ac ($H3K27ac);
		$clsEpi->set_H3K27me3 ($H3K27me3);
		$clsEpi->set_H3K36me3 ($H3K36me3);
		$clsEpi->set_H3K4me1 ($H3K4me1);
		$clsEpi->set_H3K4me3 ($H3K4me3);
		$clsEpi->set_H3K9me3 ($H3K9me3);
		$clsEpi->set_CTCF ($CTCF);
		$clsEpi->set_avExperiments ($availableExp);
		$clsEpi->set_others ($others);
		
		#add the object to the array
		push @arrClsEpi, $clsEpi;
}


#get the number of items in the array
my $numResults = scalar @arrClsEpi;
#print the results
for (my $i=0; $i < $numResults; $i++){
	print @arrClsEpi[$i]->csv_row."\n";
}



sub reference_stage {
   
	my $params = @_[0];
	my @experiments = @{$params};
	my $clientExp = REST::Client->new();
	my $headersExp = {Accept => 'application/json'};
	my $ret=''; #consistent
	my $baseUrl = 'https://www.encodeproject.org';
	my $acumLifeStage = '';
	my $acumAge = '';
	my $swFirst = 0;
	my $numLoopsExp = 0;
	my $numLoopsRep = 0;
	
	#loop each experiment
	LBL_Exp: {
		foreach my $exp (@experiments){
			$numLoopsExp++;
			my $expUrl = $baseUrl.$exp.'?frame=object&format=json';
			$clientExp->GET($expUrl, $headersExp);
			my $responseExp = decode_json($clientExp->responseContent());
			my @replicates = @{$responseExp->{'replicates'}};
			#loop each replicate in the experiment
			foreach my $rep (@replicates) {
				$numLoopsRep++;
				my $urlReplicate = $baseUrl.$rep.'?frame=embedded&format=json';
				my $clientRep = REST::Client->new();
				$clientRep->GET($urlReplicate, $headersExp);
				my $responseRep = decode_json($clientRep->responseContent());
				#print $responseRep->{'library'}->{'biosample'}->{'life_stage'}.", ".$responseRep->{'library'}->{'biosample'}->{'age'}."\n";
				if ($swFirst==0){
					#first time the fields life stage and age are obtained
					$acumLifeStage = $responseRep->{'library'}->{'biosample'}->{'life_stage'};
					$acumAge = $responseRep->{'library'}->{'biosample'}->{'age'};
					$ret = $acumLifeStage.', '.$acumAge;
					$swFirst = -1;
				}else {
					#compare the fields to find if there are any difference
					if ((uc $acumLifeStage ne uc $responseRep->{'library'}->{'biosample'}->{'life_stage'}) || (uc $acumAge ne uc $responseRep->{'library'}->{'biosample'}->{'age'})){
						#not consistent
						$ret='';
						last LBL_Exp;
					}
				
				}
			}
		}	
		
	}

	return $ret;		
}


