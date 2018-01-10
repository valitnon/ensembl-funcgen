use strict;
use REST::Client;
use JSON;
use Data::Dumper;
use feature qw(switch);
use Bio::EnsEMBL::Funcgen::clsEpiRR;

#----------------
#Check parameter
#----------------
my $url = $ARGV[0];

if (!$url) {
	usage();
	exit;
}

#--------------------------------------------
#Add parameters to get the JSON object format
#--------------------------------------------
$url.='&frame=object&format=json';

#-----------------------------------------------
#Init REST client and array for clsEpiRR objects
#-----------------------------------------------
my $clientSum = REST::Client->new();
my @arrClsEpi; #array to store the clsEpiRR objects


#----------------------
#create columns headers
#----------------------
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

#----------------------------
#add the object to the array
#----------------------------
push @arrClsEpi, $clsEpi;

#-------------
#Call the Api
#-------------
my $headers = {Accept => 'application/json'};
$clientSum->GET($url, $headers);
my $responseSum = decode_json($clientSum->responseContent());

#----------------
#loop over series
#----------------	
foreach my $serie (@{$responseSum->{'@graph'}}){
	my $project = get_project($serie->{'award'});
	my $species = $serie->{'organism'}[0];
	$species = (split '/', $species)[-1];
	my $accession = $serie->{'accession'};
	my $ihec = $serie->{'dbxrefs'}[0];
	$ihec = (split ':', $ihec)[-1];
	my $epigenome = $serie->{'biosample_term_name'}[0];
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
	
	#-------------------
	#look for DNase-Seq
	#-------------------
	foreach my $assayTerm (@{$serie->{'assay_term_name'}}){
		if (uc $assayTerm eq 'DNASE-SEQ') {
			$DNaseSeq = 'YES';
		}
		if (uc $assayTerm eq 'ATAC-SEQ') {
			$atacSeq = 'YES';
		}
	}
	
	#---------------------------------------------------------
	#find Histone modification, Chip-Seq Input and count them
	#Also checks for reference satge consistance
	#---------------------------------------------------------
	my $numH3K27ac = 0;
	my $numH3K27me3 = 0;
	my $numH3K36me3 = 0;
	my $numH3K4me1 = 0;
	my $numH3K4me3 = 0;
	my $numH3K9me3 = 0;
	my $numControl = 0;
	
	my ($refStage, $lstTargets) = get_avilable_targets($serie->{'related_datasets'});
	foreach my $target (@{$lstTargets}){
		given (uc $target){
			when ('H3K27AC') {
				$H3K27ac = 'YES';
				if ($numH3K27ac == 0){
					$numH3K27ac=1;
				}
			}
			when ('H3K27ME3') {
				$H3K27me3 = 'YES';
				if ($numH3K27me3 == 0){
					$numH3K27me3=1;
				}
			}
			when ('H3K36ME3') {
				$H3K36me3 = 'YES';
				if ($numH3K36me3 == 0){
					$numH3K36me3=1;
				}
			}
			when ('H3K4ME1') {
				$H3K4me1 = 'YES';
				if ($numH3K4me1 == 0){
					$numH3K4me1=1;
				}
			}
			when ('H3K4ME3') {
				$H3K4me3 = 'YES';
				if ($numH3K4me3 == 0){
					$numH3K4me3=1;
				}
			}
			when ('H3K9ME3') {
				$H3K9me3 = 'YES';
				if ($numH3K9me3 == 0){
					$numH3K9me3=1;
				}
			}	
			when ('CTCF') {
				$CTCF = 'YES';
			}
			when ('CONTROL') {
				$chipSeqInput = 'YES';
				if ($numControl == 0){
					$numControl=1;
				}
			}
			default {
				if ($others ne '') {
					$others .= ', ';
				}
				$others .= $target;
			}
		}
		
	}
	$availableExp = $numH3K27ac + $numH3K27me3 + $numH3K36me3 + $numH3K4me1 + $numH3K4me3 + $numH3K9me3 + $numControl;

	#-------------------------------------
	#load the data in the clsEpiRR object
	#-------------------------------------		
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
	
	#----------------------------
	#add the object to the array
	#----------------------------	
	push @arrClsEpi, $clsEpi;
	
}

#------------------
#print the results
#------------------
foreach my $REpi (@arrClsEpi){
	print $REpi->csv_row."\n";
}

sub get_project {
	my $url = shift;
	my $clientAwd = REST::Client->new();
	my $headersAwd = {Accept => 'application/json'};
	my $baseUrl = 'https://www.encodeproject.org';
	my $awdUrl = $baseUrl.$url.'?frame=object&format=json';
	$clientAwd->GET($awdUrl, $headersAwd);
	my $responseAwd = decode_json($clientAwd->responseContent());
	my $project = "";
	if ($responseAwd->{'project'}){
		$project = $responseAwd->{'project'};
	}
	
	return $project;
}

sub get_avilable_targets {
	my $experiments = shift;
	my @lstExperiments = @{$experiments};
	my $clientExp = REST::Client->new();
	my $headersExp = {Accept => 'application/json'};
	my $baseUrl = 'https://www.encodeproject.org';
	my @lstReleasedExps;
	my $ret='';
	my @lifeStage;
	my @age;
	foreach my $exp (@lstExperiments){
		my $expUrl = $baseUrl.$exp.'?frame=object&format=json';
		$clientExp->GET($expUrl, $headersExp);
		my $responseExp = decode_json($clientExp->responseContent());
		#-----------------------------------------------------
		#filter for ChIP-Seq, DNase-Seq and status = released
		#-----------------------------------------------------		
		if ((uc $responseExp->{'status'} eq 'RELEASED') && ((uc $responseExp->{'assay_title'} eq 'CHIP-SEQ') || (uc $responseExp->{'assay_title'} eq 'DNASE-SEQ'))){
			#check for target of chip-seq assays
			if (uc $responseExp->{'assay_title'} eq 'CHIP-SEQ'){
				my $rawTarget = (split '/', $responseExp->{'target'})[-1];
				my $target = (split '-', $rawTarget)[0];
				if ($target){
					push @lstReleasedExps, $target;
				}
			}
			#-------------------------------------
			#consistency of reference_stage 
			#-------------------------------------			

			my @replicates = @{$responseExp->{'replicates'}};

			foreach my $rep (@replicates) {
				my $urlReplicate = $baseUrl.$rep.'?frame=embedded&format=json';
				my $clientRep = REST::Client->new();
				$clientRep->GET($urlReplicate, $headersExp);
				my $responseRep = decode_json($clientRep->responseContent());
				
				#add LifeStage if is not been already added
				if ( !grep { $_ eq $responseRep->{'library'}->{'biosample'}->{'life_stage'}} @lifeStage )
				{
				  push @lifeStage, $responseRep->{'library'}->{'biosample'}->{'life_stage'};
				}

				#add Age if is not been already added
				if ( !grep { $_ eq $responseRep->{'library'}->{'biosample'}->{'age'}} @age )
				{
				  push @age, $responseRep->{'library'}->{'biosample'}->{'age'};
				}					
			}				

		}		
	}
	
	#-------------------------------------
	#check consistency of reference_stage 
	#-------------------------------------		
	my $numLifeStage = @lifeStage;
	my $numAge = @age;
	if ($numLifeStage > 1 || $numAge > 1){
		$ret='WARNING -- not consistent';
		if ($numLifeStage > 1){
			$ret .= ' -- Life Stages: ';
			foreach my $lStage (@lifeStage){
				$ret .= $lStage.' - '
			}
		}
		if ($numAge > 1){
			$ret .= ' -- Ages: ';
			foreach my $age (@age){
				$ret .= $age.' - '
			}
		}
	}else{
		$ret = @lifeStage[0].', '.@age[0];
	}
	
	return $ret, \@lstReleasedExps;
}

sub usage {
    my $usage = << 'END_USAGE';

Usage: encode_reference_epigenomes.pl url [look_for_reference_stage_consistency]

Options:
url: The encode search url
 
END_USAGE

    say $usage;

    return 1;
}


