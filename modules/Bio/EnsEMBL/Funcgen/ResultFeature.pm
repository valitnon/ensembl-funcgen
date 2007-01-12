#
# Ensembl module for Bio::EnsEMBL::Funcgen::ResultFeature
#
# You may distribute this module under the same terms as Perl itself

=head1 NAME

Bio::EnsEMBL::Funcgen::ResultFeature - A module to represent a lightweight ResultFeature object

=head1 SYNOPSIS

use Bio::EnsEMBL::Funcgen::ResultFeature;

my $rfeature = Bio::EnsEMBL::Funcgen::ResultFeature->new_fast
			  (	{
				 score => $score,
				 start => $start,
				 end   => $end,
				});

my @rfeatures = @{$rset->get_displayable_ResultFeature_by_Slice($slice)};

foreach my $rfeature (@rfeatures){
    my $score = $rfeature->score();
    my $rf_start = $rfeature->start();
    my $rf_end = $rfeature->end();
}

=head1 DESCRIPTION

This is a very sparse class designed to be as lightweight as possible to enable fast rendering in the web browser.
As such only the information absolutely required is contained.  Any a piori information is omitted e.g. seq_region_id, 
this will already be known as ResultFeatures are retrieved via a Slice method in ResultSet via the ResultSetAdaptor, 
likewise with analysis and experimental_chip information.  ResultFeatures are transient objects, in that they are not 
stored in the DB, but are a very small subset of information from the result and oligo_feature tables. ResultFeatures 
should only be generated by the ResultSetAdaptorast here is no parameter checking in place.


=head1 AUTHOR

This module was written by Nathan Johnson.

=head1 CONTACT

Post comments or questions to the Ensembl development list: ensembl-dev@ebi.ac.uk

=head1 METHODS

=cut

use strict;
use warnings;

package Bio::EnsEMBL::Funcgen::ResultFeature;

use Bio::EnsEMBL::Utils::Argument qw( rearrange ) ;
use Bio::EnsEMBL::Utils::Exception qw( throw warning );

#use vars qw(@ISA);
#@ISA = qw(Bio::EnsEMBL::Storable);???????????????????????????????????????????


=head2 new_fast

  Args       : Hashref with all internal attributes set
  Example    : none
  Description: Quick and dirty version of new. Only works if the code is very
               disciplined.
  Returntype : Bio::EnsEMBL::Funcgen::ResultFeature
  Exceptions : None
  Caller     : ResultSetAdaptor
  Status     : At Risk

=cut

sub new_fast {
   my ($class, $hashref)  = @_;
   return bless ($hashref, $class);
}



=head2 score

  Example    : my $score = $rf->score();
  Description: Getter of the score attribute for ResultFeature
               objects
  Returntype : string
  Exceptions : None
  Caller     : General
  Status     : At Risk

=cut

sub score {
    my $self = shift;
    return $self->{'score'};
}

=head2 start

  Example    : my $start = $rf->start();
  Description: Getter of the start attribute for ResultFeature
               objects.
  Returntype : int
  Exceptions : None
  Caller     : General
  Status     : At Risk

=cut

sub start {
    my $self = shift;
	return $self->{'start'};
}


=head2 end

  Example    : my $start = $rf->end();
  Description: Getter of the end attribute for ResultFeature
               objects.
  Returntype : int
  Exceptions : None
  Caller     : General
  Status     : At Risk

=cut

sub end {
    my $self = shift;
	return $self->{'end'};
}

1;

