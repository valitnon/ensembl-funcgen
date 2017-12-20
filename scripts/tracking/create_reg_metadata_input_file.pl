use strict;
use clsRegisterMetadata;
use REST::Client;
use JSON;
use Data::Dumper;
use Getopt::Long qw(GetOptions);
use feature qw(say);
use Config::Tiny;
use Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor;

##################### MAIN FUNCTION ##############################

#Get Parameters
my $pathFile;
my $localFilePath;
my @lstTargets;
my @lstAssays;
my $help;
my $cfgFile;
my @lstErrors;

GetOptions(
			'f=s' => \$pathFile, 
			'p=s' => \$localFilePath, 
			'a=s' => \@lstAssays,
			'c=s' => \$cfgFile,
			't=s' => \@lstTargets, 
			'h=s' => \$help,);

#check parameters			
if ( $help || !$pathFile || !$localFilePath || !@lstAssays || !$cfgFile) {
    usage();
    exit;
}
if (not(-f $pathFile)){
	#file not found
	die "$pathFile file can not be found\n";
}

if (not(-f $cfgFile)){
	#file not found
	die "$cfgFile file can not be found\n";
}

#read config file
my $cfg = Config::Tiny->read($cfgFile);

#get adaptors
my $adaptors = fetch_adaptors($cfg);

#create hashes to compare and normalize the data
my ($hAnalysis, $hFeatureType, $hGender) = get_compare_hashes($adaptors);
my %hshAnalysis = %{$hAnalysis};
my %hshFeatureType = %{$hFeatureType};
my %hshGender = %{$hGender};

#compare parameter values with DB values (assays and targets) to check if they exists and to normalize them
my @lstDbAssays;
foreach my $assayVal (@lstAssays){
	my $dbAssay = check_db_value($assayVal, \%hshAnalysis);
	if (!$dbAssay){
		push @lstErrors, "parameter: assay\tvalue: $assayVal\terror: Not found in analysis table";
	}else {
		push @lstDbAssays, $dbAssay;
	}
}

my @lstDbTargets;
foreach my $targetVal (@lstTargets){
	my $dbTarget = check_db_value($targetVal, \%hshFeatureType);
	if (!$dbTarget){
		push @lstErrors, "parameter: target\tvalue: $targetVal\terror: Not found in feature_type table";
	}else {
		push @lstDbTargets, $dbTarget;
	}
}

#if there are any errors in the parameters exits and show the errors
if (@lstErrors){
	foreach my $errorVal (@lstErrors){
		print $errorVal."\n";
	}
	exit;
}

#add the final bar if necessary
if ((substr $localFilePath, -1) ne '/'){
	$localFilePath.='/';
}

my @lstRegMeta;
my @lstCtlClsRegMeta;

#Create a clsRegisterMetadata object with the column headers and add it to the array
my $clsRegMetaHead = store_row('accession', 'epigenome', 'feature_type', 'biological_replicate', 'technical_replicate', 'gender', 'md5_checksum', 'local_url', 'analysis', 'experimental_group', 'assay_xrefs', 'ontology_xrefs', 'xrefs', 'epigenome_description', 'control_id', 'download_url', 'info');
push @lstRegMeta, $clsRegMetaHead;

#open file
open(my $fh, '<:', $pathFile) or die "Could not open file $pathFile\n";

#read file 
while (my $fullRow = <$fh>) {
	chomp $fullRow;
  	my @row = (split "\t", $fullRow);
  	my $epiAccess = @row[0];
	
	#init varaiables
	my $accession= '-';
	my $epigenome = '-';
	my $featureType = '-';
	my $bioReplicate = '-';
	my $techReplicate = '-';
	my $gender = '-';
	my $md5Check = '-';
	my $localUrl = '-';
	my $analysis = '-';
	my $expGroup = 'ENCODE'; #fix
	my $assayXrefs ='-';
	my $ontXrefs = '-';
	my $xrefs = '-';
	my $epiDesc = '-';
	my $controlId = '-';
	my $downUrl = '-';
	my $info = '-';
	my $fileName ='';
	
	my $epiFeature='';
	my $swBioRep = 1;
	my $swTecRep = 1;
	
	#get epigenome data
	my $rEData = get_reference_epigenome_data($epiAccess);
	if ($rEData){
		if ($rEData->{'dbxrefs'}[0]){
			$xrefs = $rEData->{'dbxrefs'}[0];
		}
		
		#Experiments
		my @experiments = @{$rEData->{'related_datasets'}};
		
		#get control ID
		my ($controlAccession, $infoControl) = get_control_id(\@experiments);
		
		my @lstControls;
		#loop experiments
		foreach my $exp (@experiments) {
			my $target;
			if ($exp->{'target'}->{'label'}){
				$target = $exp->{'target'}->{'label'};
			}

			my $swTarget =0;
			if (@lstDbTargets){
				#Target parameter has been informed
				if ($target){
					if ( not ((grep( /^$target$/, @lstDbTargets) ) || (uc $target eq 'CONTROL') || (uc $target eq 'RABBIT-IGG-CONTROL')) ) {
						#it is not a control or one of our targets
						$swTarget =1;
					} 
				}				
			}
			
			#change Control target to WCE
			if ((uc $target eq 'CONTROL') || (uc $target eq 'RABBIT-IGG-CONTROL')){
				$target = 'WCE';
			}
			
			#filter experiments by assay
			my $assay = $exp->{'assay_title'};
			if ((grep( /^$assay$/i, @lstDbAssays)) && ($swTarget==0)){
				#normalize the assay term
				$assay = check_db_value($assay, \%hshAnalysis);

				#check and normalize target
				if ($target){
					my $DBtarget = check_db_value($target, \%hshFeatureType);
					if (!$DBtarget){
						push @lstErrors, "experiment: ".$exp->{'accession'}."\ttarget: $target\terror: Not found in feature_type table";
					}else{
						$target = $DBtarget;
					}
				}else{
					#no target, check if assay is DNase-seq
					if (uc $assay eq 'DNASE-SEQ'){
						$target = 'DNase1';
					}else{
						#Error no target.
						push @lstErrors, "experiment: ".$exp->{'accession'}."\ttarget: Empty\terror: No target in ENCODE data";
					}
					
				}
				if ($target){
					$featureType = $target;
					if (uc $target ne 'WCE'){
						$analysis = $exp->{'assay_title'};
						#check and normalize analysis
						my $DBAnalysis = check_db_value($analysis, \%hshAnalysis);
						if (!$DBAnalysis){
							push @lstErrors, "experiment: ".$exp->{'accession'}."\tanlysis: $analysis\terror: Not found in analysis table";
						}
						$analysis = $DBAnalysis;
						$controlId = $controlAccession;
						$info ='-';
					}else{
						$analysis = '-';
						$controlId = '-';
						$info = $infoControl;
					}
				}
				
				$ontXrefs = $exp->{'biosample_term_id'};
				#replicates (Biosamples)
				my @replicates = @{$exp->{'replicates'}};
				my $rep = @replicates[0]; #get the information from the first replicate
				my $termName = '-';
				if (uc $exp->{'biosample_term_name'} ne 'UNKNOWN'){
					$termName = $exp->{'biosample_term_name'};
				}
				my $lifeStage = '-';
				if (uc $rep->{'library'}->{'biosample'}->{'life_stage'} ne 'UNKNOWN'){
					$lifeStage = $rep->{'library'}->{'biosample'}->{'life_stage'};
				}
				my $age = '-';
				if (uc $rep->{'library'}->{'biosample'}->{'age'} ne 'UNKNOWN'){
					$age = $rep->{'library'}->{'biosample'}->{'age'};
				}						
				$assayXrefs = $exp->{'assay_term_id'};
				$epigenome = join (':', $termName, $lifeStage, $age);	
				$epiDesc = join (' ', $termName, $lifeStage, $age);
				
				#Normalize gender
				$gender = check_db_value($rep->{'library'}->{'biosample'}->{'sex'}, \%hshGender);
				
				
				
				#files
				my @files = @{$exp->{'files'}};
				foreach my $file (@files){
					#get only fastq files
					if (uc $file->{'file_format'} eq 'FASTQ'){
						#accession
						$accession = $file->{'accession'};
						#md5 checksum
						$md5Check = $file->{'md5sum'};
						$bioReplicate = $file->{'biological_replicates'}[0];
						$techReplicate = (split '_', $file->{'technical_replicates'}[0])[-1];
						
						#check replicates numeration
						if ($epigenome.$target ne $epiFeature){
							#epigenome + target has changed
							if ($swBioRep == 0){
								#Error
								push @lstErrors, "experiment: ".$exp->{'accession'}."\tReplicates: Biological\terror: Numeration error found";
							}
							if ($swTecRep == 0){
								#Error
								push @lstErrors, "experiment: ".$exp->{'accession'}."\tReplicates: Technical\terror: Numeration error found";
							}
							$swBioRep=0;
							$swTecRep=0;
							$epiFeature = $epigenome.$target;
						}

						if ($swBioRep == 0){
							if ($bioReplicate == 1){
								#first Bioreplicate found
								$swBioRep=1;
							}
						}
						if ($swTecRep == 0){
							if ($techReplicate==1){
								#first Techreplicate found
								$swTecRep=1;
							}
						}
						
						#dowload url
						if ($file->{'href'}){
							$downUrl = 'https://www.encodeproject.org'.$file->{'href'};
							$localUrl = $localFilePath.(split '/', $downUrl)[-1];
						}
						
						if (uc $target eq 'WCE'){
							#assign values to object Register Metadata
							my $clsCtlRegMeta = store_row($accession, $epigenome, $featureType, $bioReplicate, $techReplicate, $gender, $md5Check, $localUrl, $analysis, $expGroup, $assayXrefs, $ontXrefs, $xrefs, $epiDesc, $controlId, $downUrl, $info);
							#Add the object to the list of control rows
							push @lstControls, $clsCtlRegMeta;
						}else{
							#assign values to object Register Metadata
							my $clsRegMeta = store_row($accession, $epigenome, $featureType, $bioReplicate, $techReplicate, $gender, $md5Check, $localUrl, $analysis, $expGroup, $assayXrefs, $ontXrefs, $xrefs, $epiDesc, $controlId, $downUrl, $info);
							#Add the object to the list of rows
							push @lstRegMeta, $clsRegMeta;
						}
					}
				}
					
			}

		}
		#add control rows
		push (@lstRegMeta, @lstControls);
	}
}

#close file
close($fh);

#print data or Errors
if (@lstErrors){
	foreach my $errorVal (@lstErrors){
		print $errorVal."\n";
	}
}else{
	foreach my $reg (@lstRegMeta) {
		print $reg->csv_row."\n";
	}
}

##################### END MAIN FUNCTION ##############################

sub get_control_id{
	my $experiments = shift;
	my $accession; 
	my $info;
	my $numExp = 0;
	foreach my $exp (@{$experiments}){
		if (uc $exp->{'target'}->{'label'} eq 'CONTROL'){
			$numExp ++;
			$info .='Exp:'.$exp->{'accession'};
			my @files = @{$exp->{'files'}};
			foreach my $file (@files){
				if (uc $file->{'file_format'} eq 'FASTQ'){
					$info .=' File:'.$file->{'accession'}.'; ';
					if ($accession ne ''){
						$accession.=',';
					}
					$accession.= $file->{'accession'};
				}
			} 
			
		}
		
	}
	if ($numExp < 2){
		$info = '-';
	}
	return $accession, $info;
}

sub get_reference_epigenome_data{
	my $refEpiAccession = shift;
	my $url = 'https://www.encodeproject.org/reference-epigenomes/'.$refEpiAccession.'/?frame=embedded&format=json';
	my $headers = {Accept => 'application/json'};
	my $rEClient = REST::Client->new();
	$rEClient->GET($url, $headers);
	my $rEResponse = decode_json($rEClient->responseContent());
	return $rEResponse;
}


sub store_row{
	my $accession= shift;
	my $epigenome = shift;
	my $featureType = shift;
	my $bioReplicate = shift;
	my $techReplicate = shift;
	my $gender = shift;
	my $md5Check = shift;
	my $localUrl = shift;
	my $analysis = shift;
	my $expGroup = shift; 
	my $assayXrefs = shift;
	my $ontXrefs = shift;
	my $xref = shift;
	my $epiDesc = shift;
	my $controlId = shift;
	my $downUrl = shift; 
	my $info = shift;
	my $clsRegMeta = clsRegisterMetadata->new();
	$clsRegMeta->set_accession($accession);
	$clsRegMeta->set_epigenome($epigenome);
	$clsRegMeta->set_feature_type($featureType);
	$clsRegMeta->set_biological_replicate($bioReplicate);
	$clsRegMeta->set_technical_replicate($techReplicate);
	$clsRegMeta->set_gender($gender);
	$clsRegMeta->set_md5_checksum ($md5Check);
	$clsRegMeta->set_local_url ($localUrl);
	$clsRegMeta->set_analysis ($analysis);
	$clsRegMeta->set_experimental_group ($expGroup);
	$clsRegMeta->set_assay_xrefs ($assayXrefs);
	$clsRegMeta->set_ontology_xrefs ($ontXrefs);
	$clsRegMeta->set_xrefs ($xref);
	$clsRegMeta->set_epigenome_description ($epiDesc);
	$clsRegMeta->set_control_id ($controlId);
	$clsRegMeta->set_download_url ($downUrl);
	$clsRegMeta->set_info ($info);
	return $clsRegMeta;
}

sub usage {
    my $usage = << 'END_USAGE';

Usage: createRegMetaInputFile.pl -f <source_file> -p <local_store_path> -a <assay> -c <config_file> [-t <target>]

Options:
-f source_file:		this is the file that contains the accessions of the Epigenome Summaries
-p local_store_path:	this is the path where the files will be allocated
-a assay:		this is the analysis description to filter experiments
-c config_file:		this is the configuration file that contains the database connection details
-h help:		shows this help message
[-t] target:	this is the feature_type to filter by
 
END_USAGE

    say $usage;

    return 1;
}

sub fetch_adaptors {
    my ($cfg) = @_;
    my %adaptors;

    my $dbAdaptor = Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor->new(
        -user    => $cfg->{efg_db}->{user},
        -host    => $cfg->{efg_db}->{host},
        -port    => $cfg->{efg_db}->{port},
        -dbname  => $cfg->{efg_db}->{dbname}
    );
	
	$adaptors{analysis} = $dbAdaptor->get_adaptor("Analysis");
	$adaptors{FeatureType} = $dbAdaptor->get_adaptor("FeatureType");
	
	
    #my $dba = $dbAdaptor->db();

    #$adaptors{epigenome}    = $dba->get_EpigenomeAdaptor();
    #$adaptors{feature_type} = $dba->get_FeatureTypeAdaptor();
    #$adaptors{analysis}     = $dba->get_AnalysisAdaptor();

    return \%adaptors;
}

sub get_compare_hashes(){
	my %adaptors = @_;
	my $lstAnalysis = $adaptors->{analysis}->fetch_all();
	my $lstFeatureTypes = $adaptors->{FeatureType}->fetch_all();
	
	my %hshAnalysis;
	foreach my $objAnalysis (@{$lstAnalysis}){
		my $logName = $objAnalysis->logic_name();
		$hshAnalysis{uc $logName}=$logName;
	}
	
	my %hshFeatureType;
	foreach my $objFeatureType (@{$lstFeatureTypes}){
		my $ftNamne = $objFeatureType->name;
		$hshFeatureType{uc $ftNamne}=$ftNamne;
	} 
	
	my %hshGender;
	$hshGender{'MALE'}='male';
	$hshGender{'FEMALE'}='female';
	$hshGender{'HERMAPHRODITE'}='hermaphrodite';
	$hshGender{'MIXED'}='mixed';
	$hshGender{'UNKNOWN'}='unknown';
	
	return \%hshAnalysis, \%hshFeatureType, \%hshGender;
}

sub check_db_value(){
	my $value= shift;
	my $hCompare = shift;
	my $ret;
	my %hshCompare = %{$hCompare};
	
	if ($hshCompare{uc $value}){
		#exists
		$ret = $hshCompare{uc $value};
	}
	return $ret;
	
}