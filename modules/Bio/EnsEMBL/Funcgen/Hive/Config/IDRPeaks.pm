=head1 LICENSE

Copyright [1999-2016] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

     http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

=head1 CONTACT

  Please email comments or questions to the public Ensembl
  developers list at <http://lists.ensembl.org/mailman/listinfo/dev>.

  Questions may also be sent to the Ensembl help desk at
  <http://www.ensembl.org/Help/Contact>.

=head1 NAME

    Bio::EnsEMBL::Funcgen::Hive::Config::IDRPeaks;

=head1 SYNOPSIS


=head1 DESCRIPTION


=head1 CONTACT

    Please contact http://lists.ensembl.org/mailman/listinfo/dev mailing list with questions/suggestions.

=cut
package Bio::EnsEMBL::Funcgen::Hive::Config::IDRPeaks;

use strict;
use warnings;

use base ('Bio::EnsEMBL::Funcgen::Hive::Config::BaseSequenceAnalysis');

sub pipeline_wide_parameters {
  my $self = shift;
  return {
    %{$self->SUPER::pipeline_wide_parameters},
    
    can_run_SWEmbl_R0005_replicate => 1,#'IDRPeaks',
    can_PreprocessIDR              => 1,#'IDRPeaks',
    can_DefineMergedDataSet        => 0, 
   };
}

sub pipeline_analyses {
  my $self = shift;

  return [
   @{$self->SUPER::pipeline_analyses}, #To pick up BaseSequenceAnalysis-DefineMergedOutputSet
    {
     -logic_name    => 'run_SWEmbl_R0005_replicate',  #SWEmbl permissive
     -module        => 'Bio::EnsEMBL::Funcgen::Hive::RunPeaks',
     -rc_name => 'normal_5GB_2cpu_monitored', # Better safe than sorry... size of datasets tends to increase...       
    },
     {
     -logic_name    => 'PreprocessIDR',
     -module        => 'Bio::EnsEMBL::Funcgen::Hive::PreprocessIDR',
     -analysis_capacity => 100,#Unlikely to get anywhere near this
     -rc_name    => 'default',
     -batch_size => 30, #Should really take ~1min to process each set of replicates
     -parameters => { permissive_peaks => $self->o('permissive_peaks') },
     -flow_into => {
       '2->A' => [ 'RunIDR' ],
       'A->3' => [ 'PostProcessIDRReplicates' ], 
      }, 
    },
    {
     -logic_name    => 'RunIDR',
     -module        => 'Bio::EnsEMBL::Funcgen::Hive::RunIDR',
     -analysis_capacity => 100,
     -rc_name    => 'normal_2GB',
     -batch_size => 6,
    -flow_into => {
      2 => [ ':////accu?idr_peak_counts=[accu_idx]' ],
     }
    },
    {
     -logic_name    => 'PostProcessIDRReplicates',
     -module        => 'Bio::EnsEMBL::Funcgen::Hive::PostprocessIDR',
     -analysis_capacity => 100,
     -rc_name    => 'default', #<1mins
     -batch_size => 10,#?
     -flow_into => {
       '2' => [ 'DefineMergedReplicateResultSet' ],
      }, 
    },
    {
     -logic_name    => 'DefineMergedReplicateResultSet',
     -module        => 'Bio::EnsEMBL::Funcgen::Hive::DefineMergedReplicateResultSet',
     -analysis_capacity => 100,
     -rc_name => 'default',
     -flow_into => { '2' => [ 'DefineMergedDataSet' ] },
    },
  ];
}

1;
