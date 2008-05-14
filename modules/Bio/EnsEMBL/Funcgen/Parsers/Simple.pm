#
# EnsEMBL module for Bio::EnsEMBL::Funcgen::Parsers::Simple
#

=head1 NAME

Bio::EnsEMBL::Funcgen::Parsers::Simple

=head1 SYNOPSIS

  use vars qw(@ISA);
  @ISA = qw(Bio::EnsEMBL::Funcgen::Parser::Simple);

  
=head1 DESCRIPTION

This is a base class to support simple file format parsers.

=head1 AUTHOR

This module was created by Nathan Johnson.

=head1 CONTACT

Post questions to the EnsEMBL development list ensembl-dev@ebi.ac.uk

=head1 METHODS

=cut

package Bio::EnsEMBL::Funcgen::Parsers::Simple;

use Bio::EnsEMBL::Funcgen::ExperimentalSet;
use Bio::EnsEMBL::Funcgen::FeatureSet;
use Bio::EnsEMBL::Funcgen::AnnotatedFeature;
use Bio::EnsEMBL::Utils::Exception qw( throw warning deprecate );
use Bio::EnsEMBL::Funcgen::Utils::EFGUtils qw(species_chr_num open_file);
use Bio::EnsEMBL::Utils::Argument qw( rearrange );
use Bio::EnsEMBL::Funcgen::Utils::Helper;
use strict;


use vars qw(@ISA);
@ISA = qw(Bio::EnsEMBL::Funcgen::Utils::Helper);

=head2 new

  Example    : my $self = $class->SUPER::new(@_);
  Description: Constructor method for Bed class
  Returntype : Bio::EnsEMBL::Funcgen::Parsers::Simple
  Exceptions : throws if caller is not Importer
  Caller     : Bio::EnsEMBL::Funcgen::Parsers:Simple
  Status     : at risk

=cut


sub new{
  my $caller = shift;
  
  my $class = ref($caller) || $caller;
  my $self  = $class->SUPER::new(@_);
  

  #Could potentially take fields params directly to define a custom format
  #Take direct field mappings, plus special fields which needs parsing differently
  #i.e. default is tab delimited, and GFF would define Attrs field as compound field and provide special parsing and field mapping
  

  throw("This is a skeleton class for Bio::EnsEMBL::Importer, should not be used directly") 
	if(! $self->isa("Bio::EnsEMBL::Funcgen::Importer"));
  
  $self->{'config'} =  
	{(
	  #can we omit these?
      #array_data   => [],#['experiment'],
      #probe_data   => [],#["probe"],
      #norm_method => undef,
	  #protocols => {()},
	  'results_data' => ["and_import_simple"], 
     )};
	

  return $self;
}


=head2 set_config

  Example    : my $self->set_config;
  Description: Sets attribute dependent config
  Returntype : None
  Exceptions : None
  Caller     : Bio::EnsEMBL::Funcgen::Importer
  Status     : at risk

=cut


sub set_config{
  my $self = shift;

  throw('Must provide an ExperimentalSet name for a '.uc($self->vendor).' import') if ! defined $self->experimental_set_name();

  #Mandatory checks
  if(! defined $self->feature_analysis){
	throw('Must define a -feature_analysis parameter for '.uc($self->vendor).' imports');
  }


  #We need to undef norm method as it has been set to the env var
  $self->{'norm_method'} = undef;

  #dir are not set in config to enable generic get_dir method access

  return;
}

sub define_sets{
  my ($self) = @_;


  my $eset_adaptor = $self->db->get_ExperimentalSetAdaptor();
 
 

  #Use define_and_validate with append as we may have a pre-existing set
 
  my $dset = $self->define_and_validate_sets
	(
	 -dbadaptor    => $self->db,
	 -name         => $self->experimental_set_name,
	 -feature_type => $self->feature_type,
	 -cell_type    => $self->cell_type,
	 -analysis     => $self->feature_analysis,
	 -type         => 'annotated', 
	 -description  => $self->set_description,
	 -rollback     => $self->rollback,
	 -append       => 1,
	);

  #It may already be defined as we may be adding a subset to a previously imported set.
  
  my ($eset) = @{$dset->get_supporting_sets};#Only expecting one.
 
  if($eset->name ne $self->experimental_set_name){
	throw('Found pre-exiting DataSet('.$dset->name.') with incorrect supporting '
		  .ucfirst($eset->type).'Set('.$eset->name.')');
  }


  if(! defined $eset){
	$eset = Bio::EnsEMBL::Funcgen::ExperimentalSet->new
	  (
	   -name         => $self->experimental_set_name(),
	   -experiment   => $self->experiment(),
	   -feature_type => $self->feature_type(),
	   -cell_type    => $self->cell_type(),
	   -vendor       => $self->vendor(),
	   -format       => $self->format(),
	   -analysis     => $self->feature_analysis,
	  );
	($eset)  = @{$eset_adaptor->store($eset)};

	#add to dset here and store
	$dset->add_supporting_set($eset);
	$dset->adaptor->update_supporting_sets($dset);
  }

  #We could retrieve the DataSet by supporting set here, but define_and_validate_sets is more exhautive
  #Just need to check the ExperimentalSet hasn't already be added as a supporting set

  #We are now using IMPORTED to define wheather a FeatureSet has been imported succesfully
  #However we already have IMPORTED on the ExperimentalSubSet
  #We shoiuld add it to FeatureSet to remain consistent.
  #Do we not need to remove IMPORTED from FeatureSet whilst we are doing the import
  #Yes as we may try importing more, but not point at the old files
  #so this will never get noticed if we don't have set wide status
  
  return $dset;

}

 
sub read_and_import_simple_data{
  my $self = shift;
    
  $self->log("Reading and importing ".$self->vendor()." data");
  my (@header, @data, @design_ids, @lines);
  my ($anal, $fh, $file);
  
 
  my $af_adaptor = $self->db->get_AnnotatedFeatureAdaptor();
  my $new_data = 0;

  my $dset   = $self->define_sets;
  my ($eset) = @{$dset->get_supporting_sets};   

  #Get file
  if (! @{$self->result_files()}) {
	my $list = "ls ".$self->input_dir().'/'.$self->name().'*.'.lc($self->vendor);#could use vendor here?
	my @rfiles = `$list`;
	$self->result_files(\@rfiles);
  }
    
  if (scalar(@{$self->result_files()}) >1) {
	warn("Found more than one bed file:\n".
		 join("\n", @{$self->result_files()})."\nBed does not yet handle replicates.".
		 "  We need to resolve how we are going handle replicates with random cluster IDs");
	#do we even need to?
  }
  
  #Here were are tracking the import of individual bed files by adding them as ExperimentalSubSets
  #Recovery would never know what to delete
  #So would need to delete all, Hence no point in setting status?
  
  
  
  ### VALIDATE FILES ###
  #We need validate all the files first, so the import doesn't fall over half way through
  #Or if we come across a rollback halfway through
  my (%new_data, %feature_params);
  my $roll_back = 0;
  
  foreach my $filepath( @{$self->result_files} ) {
	chomp $filepath;
	my $filename;
	($filename = $filepath) =~ s/.*\///;
	my $sub_set;
	$self->log('Found '.$self->vendor." file\t$filename");
	
	if( $sub_set = $eset->get_subset_by_name($filename) ){
	  
	  if( $sub_set->has_status('IMPORTED') ){
		$self->log("ExperimentalSubset(${filename}) has already been imported");
	  } 
	  else{
		$self->log("Found partially imported ExperimentalSubset(${filename})");
		$roll_back = 1;
		
		if ( $self->recovery && $roll_back ) {
		  $self->log("Rolling back results for ExperimentalSubset:\t".$filename);
		  warn "Cannot yet rollback for just an ExperimentalSubset, rolling back entire set\n";
		  
		  $self->rollback_FeatureSet($fset);
		  $self->rollback_ExperimentalSet($eset);
		  last;
		}
		elsif( $roll_back ){
		  throw("Found partially imported ExperimentalSubSet:\t$filepath\n".
				"You must specify -recover  to perform a full roll back for this ExperimentalSet:\t".$eset->name);
		}
	  }
	}
	else{
	  $self->log("Found new ExperimentalSubset(${filename})");
	  $new_data{$filepath} = undef;
	  $sub_set = $eset->add_new_subset($filename);
	  $eset_adaptor->store_ExperimentalSubsets([$sub_set]);
	}
  }

  
  ### READ AND IMPORT FILES ###
  my ($filename);


  ###WAS HERE, MOVING MY VARS outside of loop where possible




  foreach my $filepath(@{$self->result_files()}) {
	chomp $filepath;
	my $count = 0;
	($filename = $filepath) =~ s/.*\///;
	
	if( $roll_back || exists $new_data{$filepath} ){
	  
	  $self->log('Reading '.$se;f->vendor." file:\t".$filename);
	  my $fh = open_file($filepath);
	  my @lines = <$fh>;
	  close($fh);
		 
	  my ($line, $f_out);
	  my $fasta = '';
	  
		 #warn "we need to either dump the pid rather than the dbID or dump the fasta in the DB dir";
		 my $fasta_file = $ENV{'EFG_DATA'}."/fastas/".$self->experiment->name().'.'.$filename.'.fasta';

		 if($self->dump_fasta()){

		   #we need to check whether we have a seq field accesible here

		   $self->backup_file($fasta_file);
		   $f_out = open_file($fasta_file, '>');
		 }
		
		 foreach my $line (@lines) {
		   $line =~ s/\r*\n//o;

		   #Genericise these
		   next if $line =~ /^\#/;	
		   next if $line =~ /^$/;


		   next $line !~ /^chr/i;
		   #next if $line =~ /^chr/i;#Mikkelson hack
		   
		   my ($chr, $start, $end, $pid, $score) = split/\t/o, $line;				  
		   #my ($chr, $start, $end, $score) = split/\t/o, $line;#Mikkelson hack				  
		   
		   my @fields = split/\t/o, $line;


		   #start building parameters hash
		   foreach my $field_index(@{$self->get_field_indices}){

			 

			 if($self->get_field_by_index($field_index)


		   if($self->ucsc_coords){
			 $start +=1;
			 $end +=1;
		   }
		   
		   if(!  $self->cache_slice($chr)){
			 warn "Skipping AnnotatedFeature import, cound non standard chromosome: $chr";
		   }
		   else{
			 
			 #this is throwing away the encode region which could be used for the probeset/family?	
			 my $feature = Bio::EnsEMBL::Funcgen::AnnotatedFeature->new
			   (
				-START         => $start,
				-END           => $end,
				-STRAND        => 1,
				-SLICE         => $self->cache_slice($chr),
				-ANALYSIS      => $fset->analysis,
				-SCORE         => $score,
				-DISPLAY_LABEL => $pid,
				-FEATURE_SET   => $fset,
			   );
			 
			 $af_adaptor->store($feature);
			 
			 $count++;

			 #dump fasta here
			 if ($self->dump_fasta()){
			   $fasta .= '>'.$pid."\n".$self->cache_slice($chr)->sub_Slice($start, $end, 1)->seq()."\n";
			 }
		   }
		 }
		 
		 
		 if ($self->dump_fasta()){
		   print $f_out $fasta;
		   close($f_out);
		 }


		 $self->log("Finished importing $count features:\t$filepath");
		 my $sub_set = $eset->get_subset_by_name($filename);
		 $sub_set->adaptor->store_status('IMPORTED', $sub_set);
	   }
    }

    $self->log("No new data, skipping result parse") if ! keys %new_data && ! $roll_back;   
    $self->log("Finished parsing and importing results");


	#Retrieve DataSet here as sanity check it is valid
	#We had one which had a mismatched feature_type between the feature_set and the experimental_set
	#How was this possible, surely this should have failed on generation/storage?

    
    return;
}
  


1;
