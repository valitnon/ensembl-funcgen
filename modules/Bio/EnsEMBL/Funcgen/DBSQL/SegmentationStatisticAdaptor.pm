=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2017] EMBL-European Bioinformatics Institute

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

=head1 SYNOPSIS

=head1 DESCRIPTION

=cut

package Bio::EnsEMBL::Funcgen::DBSQL::SegmentationStatisticAdaptor;

use strict;
use Carp;
use base 'Bio::EnsEMBL::Funcgen::DBSQL::GenericAdaptor';

sub object_class {
  return 'Bio::EnsEMBL::Funcgen::SegmentationStatistic';
}

sub _tables {
  return ['segmentation_statistic', 'ss']
}

# sub fetch_by_statistic_Epigenome {
# 
#   my $self         = shift;
#   my $statistic    = shift;
#   my $segmentation = shift;
#   my $epigenome    = shift;
#   
#   return $self->fetch_single_object(
#     'epigenome_id = ? and statistic = ? and segmentation_id = ?', 
#     [ 
#       $epigenome->dbID,
#       $statistic,
#       $segmentation->dbID,
#     ]
#   );
# }


sub fetch_num_epigenomes_by_Segmentation {
  my $self = shift;
  my $segmentation = shift;
  
  if (! defined $segmentation) {
    confess("Segmentation parameter is missing!");
  }
  if (! $segmentation->isa('Bio::EnsEMBL::Funcgen::Segmentation')) {
    confess("Segmentation parameter is not a segmentation!");
  }
  
  my $statistic = $self->fetch_single_object(
    'segmentation_id = ? and statistic = ?',
    [ 
      $segmentation->dbID,
      'num_epigenomes',
    ]
  );
  
  if (! defined $statistic) {
    $statistic = Bio::EnsEMBL::Funcgen::SegmentationStatistic->new(
      -statistic => 'num_epigenomes',
      -value     => 'not defined',
    );
  }
  
  return $statistic;
}

sub _delete_possibly_existing_statistics {

  my $self = shift;
  my $statistics = shift;
  
  foreach my $statistic (@$statistics) {
    $self->_delete_possibly_existing_statistic($statistic);
  }
  return;
}

sub _delete_possibly_existing_statistic {

  my $self = shift;
  my $statistic = shift;
  
  my $segmentation_clause = defined $statistic->segmentation_id ? 'segmentation_id = ?' : 'segmentation_id is null';
  my $label_clause        = defined $statistic->label           ? 'label           = ?' : 'label           is null';
  my $statistic_clause    = defined $statistic->statistic       ? 'statistic       = ?' : 'statistic       is null';
  my $epigenome_id_clause = defined $statistic->epigenome_id    ? 'epigenome_id    = ?' : 'epigenome_id    is null';

  my $sql = "
    delete from segmentation_statistic 
      where 
            $segmentation_clause
        and $label_clause
        and $statistic_clause
        and $epigenome_id_clause
  ";
  my $sth = $self->prepare( $sql );
  
  use DBI qw(:sql_types);
  my $position = 1;
  
  if (defined $statistic->segmentation_id) { $sth->bind_param($position++, $statistic->segmentation_id, SQL_INTEGER); }
  if (defined $statistic->label          ) { $sth->bind_param($position++, $statistic->label,           SQL_VARCHAR); }
  if (defined $statistic->statistic      ) { $sth->bind_param($position++, $statistic->statistic,       SQL_VARCHAR); }
  if (defined $statistic->epigenome_id   ) { $sth->bind_param($position++, $statistic->epigenome_id,    SQL_INTEGER); }

  $sth->execute;
  return;
}

sub _fetch_something_by_label {
  my $self = shift;
  my $the_thing    = shift;
  my $label        = shift;
  
  my $statistic_name = $the_thing;
  
  my $statistic = $self->fetch_single_object(
    'segmentation_id is null and label = ? and statistic = ?',
    [ 
      $label,
      $statistic_name,
    ]
  );
  if (! defined $statistic) {
    $statistic = Bio::EnsEMBL::Funcgen::SegmentationStatistic->new(
      -statistic => $statistic_name,
      -value     => 'not defined',
    );
  }
  return $statistic;
}

sub _fetch_something_by_Segmentation_label {
  my $self = shift;
  my $the_thing    = shift;
  my $segmentation = shift;
  my $label        = shift;
  
  if (! defined $segmentation) {
    confess("Segmentation parameter is missing!");
  }
  if (! $segmentation->isa('Bio::EnsEMBL::Funcgen::Segmentation')) {
    confess("Segmentation parameter is not a segmentation!");
  }

  my $statistic_name = $the_thing;
  
  my $statistic = $self->fetch_single_object(
    'segmentation_id = ? and label = ? and statistic = ?',
    [ 
      $segmentation->dbID,
      $label,
      $statistic_name,
    ]
  );

  if (! defined $statistic) {
    $statistic = Bio::EnsEMBL::Funcgen::SegmentationStatistic->new(
      -statistic => $statistic_name,
      -value     => 'not defined',
    );
  }
  
  return $statistic;
}

sub fetch_average_length_by_label {
  my $self = shift;
  return $self->_fetch_something_by_label('average_length', @_);
}

sub fetch_q0_by_label {
  my $self = shift;
  return $self->_fetch_something_by_label('length_q0', @_);
}

sub fetch_q1_by_label {
  my $self = shift;
  return $self->_fetch_something_by_label('length_q1', @_);
}

sub fetch_q2_by_label {
  my $self = shift;
  return $self->_fetch_something_by_label('length_q2', @_);
}

sub fetch_q3_by_label {
  my $self = shift;
  return $self->_fetch_something_by_label('length_q3', @_);
}

sub fetch_q4_by_label {
  my $self = shift;
  return $self->_fetch_something_by_label('length_q4', @_);
}

sub fetch_average_length_by_Segmentation_label {
  my $self = shift;
  return $self->_fetch_something_by_Segmentation_label('average_length', @_);
}

sub fetch_num_states_by_Segmentation_label {
  my $self = shift;
  return $self->_fetch_something_by_Segmentation_label('num_states', @_);
}

sub fetch_q0_by_Segmentation_label {
  my $self = shift;
  return $self->_fetch_something_by_Segmentation_label('length_q0', @_);
}

sub fetch_q1_by_Segmentation_label {
  my $self = shift;
  return $self->_fetch_something_by_Segmentation_label('length_q1', @_);
}

sub fetch_q2_by_Segmentation_label {
  my $self = shift;
  return $self->_fetch_something_by_Segmentation_label('length_q2', @_);
}

sub fetch_q3_by_Segmentation_label {
  my $self = shift;
  return $self->_fetch_something_by_Segmentation_label('length_q3', @_);
}

sub fetch_q4_by_Segmentation_label {
  my $self = shift;
  return $self->_fetch_something_by_Segmentation_label('length_q4', @_);
}

1;