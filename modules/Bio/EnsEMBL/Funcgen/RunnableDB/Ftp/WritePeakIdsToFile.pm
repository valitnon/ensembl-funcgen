package Bio::EnsEMBL::Funcgen::RunnableDB::Ftp::WritePeakIdsToFile;

use strict;
use Data::Dumper;
use Bio::EnsEMBL::Registry;
use base ('Bio::EnsEMBL::Hive::Process');

use constant {
  OUTPUT_BRANCH => 2,
};

sub run {
  my $self = shift;

  my $species                   = $self->param('species');
  my $epigenome_production_name = $self->param('epigenome_production_name');
  my $feature_type_name         = $self->param('feature_type_name');
  my $output_file               = $self->param('output_file');

  my $funcgen_adaptor = Bio::EnsEMBL::Registry->get_DBAdaptor( $species, 'funcgen' );

  my $helper = Bio::EnsEMBL::Utils::SqlHelper->new(
    -DB_CONNECTION => $funcgen_adaptor->dbc
  );

  open my $fh, '>', $output_file;
  
  my $number_of_ids_written;
  
  $helper->execute_no_return(
    -SQL      => "
      select 
        peak.peak_id 
      from 
        peak 
        join peak_calling using (peak_calling_id) 
        join epigenome using (epigenome_id) 
        join experiment using (experiment_id) 
        join feature_type on (feature_type.feature_type_id = experiment.feature_type_id) 
      where 
        epigenome.production_name = ? 
        and feature_type.name     = ?
    ",
    -PARAMS => [ 
      $epigenome_production_name, 
      $feature_type_name 
     ],
    -CALLBACK => sub {
      my $row = shift;
      my $annotated_feature_id = $row->[0];
      $fh->print($annotated_feature_id);
      $fh->print("\n");
      $number_of_ids_written++;
      return;
    },
  );
  $fh->close;
  
  my $job = $self->get_job_hash_ref;
  my $new_data_flow = {
    %$job,
    'number_of_ids_in_file' => $number_of_ids_written,
  };

  $self->dataflow_output_id(
    $new_data_flow, 
    OUTPUT_BRANCH
  );
}

sub get_job_hash_ref {

  my $self = shift;

  my $input_job = $self->input_job;
  use Bio::EnsEMBL::Hive::Utils qw( destringify );
  my $input_id = destringify($input_job->input_id);
  return $input_id;
}

1;
