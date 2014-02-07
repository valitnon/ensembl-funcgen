#
# Ensembl module for Bio::EnsEMBL::Funcgen::DBSQL::BaseAdaptor
#


=head1 LICENSE

Copyright [1999-2013] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute

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
  developers list at <ensembl-dev@ebi.ac.uk>.

  Questions may also be sent to the Ensembl help desk at
  <http://www.ensembl.org/Help/Contact>.


=head1 NAME

Bio::EnsEMBL::Funcgen::DBSQL::BaseAdaptor - Simple wrapper class for Funcgen Storable Adaptors

=head1 SYNOPSIS

$adaptor->store_states($storable);

=head1 DESCRIPTION

This is a simple wrapper class to hold common methods to all Funcgen StorableAdaptors.
Includes status methods.

=head1 SEE ALSO

Bio::EnsEMBL::DBSQL::BaseAdaptor

=cut

package Bio::EnsEMBL::Funcgen::DBSQL::BaseAdaptor;

use strict;
use warnings;
use Bio::EnsEMBL::Utils::Exception qw( throw deprecate );
use Bio::EnsEMBL::Utils::Scalar    qw( assert_ref );
use DBI qw(:sql_types);
#have to re-import them, as we lose :sql_types $DBI::EXPORT_TAGS in core BaseAdaptor
#we could just re-export everything

use base qw(Bio::EnsEMBL::DBSQL::BaseAdaptor Exporter);

require Exporter; #Still required to use vars @EXPORT
use vars qw( @EXPORT );
@EXPORT = (@{$DBI::EXPORT_TAGS{'sql_types'}});


=head2 new

  Example    : my $adaptor = 
  Description: 
  Returntype : Bio::EnsEMBL::Funcgen::BaseAdaptor
  Exceptions : None
  Caller     : Bio::EnsEMBL::Funcgen::BaseAdaptor subclass constructor methods
  Status     : At Risk

=cut

#This is currently adding a requirement for file based adaptor
#to have a _true_tables method

sub new {
  my ($class, @args) = @_;
  my $self = $class->SUPER::new(@args);
  
  $self->reset_true_tables; #Set the _tables values
  
  #Set the main table attrs to avoid having to 
  #repeatedly handle them in the contrain methods
  my @tables = $self->_tables;
  ($self->{_table_name}, $self->{_table_syn}) = @{$tables[0]};
    
  return $self;
}


=head2 compose_constraint_query

  Arg [1]    : Hash - Params hash containing 'constraints' hash key value pairs
  Example    : my @objs = @{$self->generic_fetch($self->compose_contraint_query($params_hash))};
  Description: Given a params Hashref containing a constraints hash e.g.
                {
                 constraints =>
                  {
                   names    => [ 'name1', 'name2' ],
                   analyses => [ $analysis_obj1, $analysis_obj2 ],
                  }
                }
               This method will iterate through the cosntraints keys calling
               _constrain_${constraints_key} e.g. _constrain_analyses.
               SQL is built from all the specified constraint, and tables 
               added for use when a constraint uses dynamic query composition across 
               tables outside the normal specification for this adaptor. 
               Invalid constraints are caught and a helpful list of valid constraints 
               are listed.

               NOTE: Most constraints operation with OR logic implicitly via an
               SQL join or IN clause. There is currently one exception to this
               is the 'states' constraint which uses AND logic, as this provide
               more appropriate functionality for this constraint.

  Returntype : Scalar - Constraint SQL string
  Exceptions : Throws if params hash argument is not valid.
               Throws if specific constraint name if not valid.
  Caller     : fetch methods
  Status     : At Risk

=cut

#This approach originated in the need for a fully flexible method call
#for the FeatureSetAdaptor to support the Experiment view, but is also
#now used to support flexible filtering of pipeline input sets
#TODO Add support for final_clause (order/group) and default_where

#Currently the API methods only take a constraints hash as an argument
#not a full params hash (i.e. constraints, optional_constraints or other config

#optional_parameters doesn't actually work that well
#as [undef] is not caught as it still has a size of 1
#e.g when passing optional params which are not defined e.g. [$logic_name] 

sub compose_constraint_query{
  my ($self, $params) = @_;

  if($params &&
	 (ref($params) ne 'HASH') ){
	throw('You must pass a valid params HASHREF to compose_constraint_query');
  }

  my @constraints;

  for my $con_type (qw(constraints optional_constraints)){ 

    if( exists ${$params}{$con_type} ){
  
      foreach my $constraint_key(keys (%{$params->{$con_type}})){
        #warn "$con_type $constraint_key = ".$params->{$con_type}{$constraint_key};
        my $constrain_method = '_constrain_'.$constraint_key;
          
        if(! $self->can($constrain_method)){
          throw($constraint_key." is not a valid constraint type.");
          # Valid constraints for ".ref($self)." are:\n\t".join("\n\t", @{$self->list_valid_constraints}));
        }
    
        #Only call constraint method if we have data or it is not optional
        #The following test allows empty arrayrefs to be passed for optional_cosntraints
  
        if((( (ref($params->{$con_type}{$constraint_key}) eq 'ARRAY') &&
               @{$params->{$con_type}{$constraint_key}}) || #arrayref and populated or
              defined $params->{$con_type}{$constraint_key} ) # otherwise defined (likely scalar or object)
          || $con_type ne 'optional_constraints'){ #or not optional (implicitly undefined)
  
          my ($constraint, $constraint_conf) =
            $self->$constrain_method($params->{$con_type}{$constraint_key}, $params);
          push @constraints, $constraint;
  
          #Currently only handle tables here but could also
          #set other dynamic config e.g. final_clause etc.
          if (exists ${$constraint_conf}{tables}) {    
              $self->_tables($constraint_conf->{tables});
          }
        }#else do nothing
      } # END OF constraint_keys
    }    
  } # END OF con_types
 
  return join(' AND ', @constraints) || '';
}


#=head2 list_valid_constraints
#
#  Example    : print "Valid constraints are:\t".join("\t", @{$adaptor->list_valid_constraints});
#  Description: This method simply returns a list of valid constraint hash keys for use with
#               fetch method which support the compose_query_constraint method.
#  Returntype : Listref of scalar contraint keys
#  Exceptions : None
#  Caller     : compose_query_constraint
#  Status     : At Risk
#
#=cut

#Grab all the keys from the symbol table and see if they have a coderef assigned
#no strict 'refs';
#my @methods = grep { defined &{$_} && ($_ =~ /^_constrain_.*/)} keys %{ref($self).'::'};
#use strict 'refs';
#This doesn't recurse down @ISA and don't seem to have a coderef associated???
#although it is supposed to

#sub list_valid_constraints{
#  my $self = $_[0];

#  #Minor helpful voodoo
#  return [ map {/_constrain_([a-zA-Z_]+$)/ ? $1 : () }
#               @{Class::Inspector->methods(ref($self), 'private')} ];
#}


#Class::Inspector is no a core perl module


=head2 _tables

  Args       : None
  Example    : my @adaptor_tables = @{$adaptor->_tables};
  Description: PROTECTED implementation of superclass abstract method.
               Returns the names and aliases of the tables as defined by the
               _true_tables method in the subclass and adds tables as required
               for query composition
  Returntype : List of listrefs of strings
  Exceptions : None
  Caller     : Internal
  Status     : At Risk

=cut

#todo this currently doesn't account for tables being added more than once between constraints
#This could be handled quite easily with a hash
#but would not solved the problem of duplicating the join clause in the constraint method
#would be sensible to be key on alias rather than table name, 
#to allow multiple joins to the same table

sub _tables {
  my ($self, $new_tables) = @_;
  
  if(defined $new_tables){
    push @{$_[0]->{_tables}}, @{$new_tables};  
  }
  
  return @{$_[0]->{_tables}};
}


=head2 reset_true_tables

  Args       : None
  Example    : $adaptor->reset_true_tables;
  Description: Resets the tables attributes to the default or 'true' tables
               for a given composable adaptor.
  Returntype : None
  Exceptions : None
  Caller     : Internal
  Status     : At Risk

=cut

#change to _reset_true_tables?

sub reset_true_tables{
  my $self = $_[0];
  @{$self->{_tables}} = $self->_true_tables;
  return;
}


=head2 _main_table

  Example    : my $syn = $adaptor->_main_table->[1];
  Description: Convenience method to retrieve the main table or main table synonym for this adaptor
               Entirely dependent on ensembl convention of always having main table as first element
               of tables array.
  Returntype : Array ref
  Exceptions : None
  Caller     : General
  Status     : At Risk

=cut

#set this in new and return attr

sub _main_table{
  my $self = $_[0];

  #Need to do this to put it in list context to avoid just returning the last value
  my @tables = $self->_tables;
  return $tables[0];
}

#why is this calling SUPER::_tables rather than the _tables method in this module?
#This is due to the order of @ISA in the BaseFeatureAdaptor


sub _table_syn {
  return $_[0]->{_table_syn};  
}

sub _table_name {
 return $_[0]->{_table_name};   
}



=head2 store_states

  Arg [1]    : Bio::EnsEMBL::Funcgen::Storable
  Example    : $rset_adaptor->store_states($result_set);
  Description: Stores states of Storable in status table.
  Returntype : None
  Exceptions : Throws if Storable is not stored
  Caller     : General
  Status     : At Risk

=cut

sub store_states{
  my ($self, $storable) = @_;
  assert_ref($storable, 'Bio::EnsEMBL::Funcgen::Storable');

  foreach my $state(@{$storable->get_all_states()}){
    
    if (! $self->has_stored_status($state, $storable)){
      $self->store_status($state, $storable) 
    }
  }

  return;
}


=head2 fetch_all

  Arg[1]     : Hashref - optional parameters e.g. {constraints => {states => ['DISPLAYABLE']}}
               or (for backwards compatibility) a status string e.g. 'DISPLAYABLE'.
  Example    : my @dsets = @{$dsa->fetch_all()};
  Description: Gets all available objects from the DB modulated by the optional parameters argument.
  Returntype : Arrayref
  Exceptions : None
  Caller     : General
  Status     : At Risk

=cut

sub fetch_all {
  my ($self, $params) = @_;

  if($params){  # Still supporting passing a status String arg for now

    if(ref($params) eq ''){ #This is a scalar status arg
      $params = {constraints => {states => [$params]}};
    }
  }

  my $objs = $self->generic_fetch($self->compose_constraint_query($params));
  $self->reset_true_tables;
   
  return $objs;
}


=head2 fetch_all_displayable

  Example    : my @displayable_dset = @{$dsa->fetch_all_displayable()};
  Description: Gets all displayable DataSets
  Returntype : ARRAYREF
  Exceptions : None
  Caller     : General
  Status     : At Risk - can we just reimplement fetch_all with optional status arg

=cut

sub fetch_all_displayable{
  my $self = shift;
  return $self->fetch_all({states => ['DISPLAYABLE']});
}


=head2 status_to_constraint

  Arg [1]    : string - status e.g. 'DISPLAYABLE'
  Arg [2]    : string - Constraint
  Example    : $sql = $self->status_to_constraint($self, $constraint, $status);
  Description: Appends the appropriate status constraint dependant on the BaseAdaptor sub class.
  Returntype : string - constraint
  Exceptions : None
  Caller     : Bio::EnsEMBL::Funcgen::DBSQL::"BaseAdaptors"
  Status     : At risk - to be deprecated

=cut

#TODO remove this in favour of compose_constraint

sub status_to_constraint{
  my ($self, $status) = @_;

  #This is now supported in compose_constraint_query
  #Which avoid some of the problems below

  my $constraint;

  #This will throw if status not valid, but still may be absent
  my $status_id = $self->_get_status_name_id($status);


  return if (! $status_id);

  my @tables = $self->_tables;
  my ($table_name, $syn) = @{$tables[0]};

  my @status_ids;

  my $sql = "SELECT table_id from status where table_name='$table_name' and status_name_id='$status_id'";
  @status_ids = map $_ = "@$_", @{$self->db->dbc->db_handle->selectall_arrayref($sql)};

  #This is causing problems as we might get none, which will invalidate the sql
  #Hence we return nothing

  $constraint = " $syn.${table_name}_id IN (".join(',', @status_ids).")" if @status_ids;
  return $constraint;
}


=head2 _test_funcgen_table

  Arg [1]    : Bio::EnsEMBL::"OBJECT"
  Example    : $status_a->_is_funcgen_object($experimental_chip)};
  Description: Tests if the object is a valid funcgen object with an identifiable table_name
  Returntype : string - table_name
  Exceptions : Throws if argument if a Bio::EnsEMBL::Funcgen::"OBJECT" not supplied
               Throws if not table name identified
  Caller     : general
  Status     : At risk

=cut

sub _test_funcgen_table{
  my ($self, $obj) = @_;

  #Can we change this to is_stored_and_valid for a Storable?
  $self->db->is_stored_and_valid('Bio::EnsEMBL::Funcgen::Storable', $obj);
  #Does this test for ad

  my @tables = $obj->adaptor->_tables;

  my ($table) = @{$tables[0]};
  #InputSubSet fix, as doesn't have own adaptor
  $table = 'input_subset' if $obj->isa('Bio::EnsEMBL::Funcgen::InputSubset');

  return $table || $self->throw("Cannot identify table name from $obj adaptor");
}


=head2 fetch_all_states

  Arg [1]    : Bio::EnsEMBL::"OBJECT"
  Arg [2...] : listref of states
  Example    : my @ec_states = @{$status_a->fetch_all_states($experimental_chip)};
  Description: Retrieves all states associated with the given "OBJECT"
  Returntype : ARRAYREF
  Exceptions : None
  Caller     : general
  Status     : At risk

=cut

sub fetch_all_states{
  my ($self, $obj) = @_;

  my $table = $self->_test_funcgen_table($obj);
  my $sql = "SELECT name FROM status_name sn, status s WHERE s.table_name='$table' ".
    "AND s.table_id='".$obj->dbID()."' and s.status_name_id=sn.status_name_id";
  my @states = map $_ = "@$_", @{$self->db->dbc->db_handle->selectall_arrayref($sql)};

  return \@states;
}


=head2 has_stored_status

  Arg [1]    : string - status e.g. IMPORTED, DISPLAYABLE
  Arg [2]    : Bio::EnsEMBL::Storable
  Example    : if($status_a->has_stored_status('IMPORTED', $array){ ... skip import ... };
  Description: Tests wether a given object has a given state
  Returntype : BOOLEAN
  Exceptions : None
  Caller     : Bio::EnsEMBL::Funcgen::BaseAdaptor
  Status     : At risk

=cut

#Only used for set_status, merge with set_status?
 
sub has_stored_status{
  my ($self, $state, $obj) = @_;

  my (@row);
  my $status_id = $self->_get_status_name_id($state);

  if ($status_id){

    my $table = $self->_test_funcgen_table($obj);

    if($status_id){
      my $sql = "SELECT status_name_id FROM status WHERE table_name=\"$table\"".
        " AND table_id=\"".$obj->dbID()."\" AND status_name_id=\"$status_id\"";

      #could just return the call directly?
      @row = $self->db->dbc->db_handle->selectrow_array($sql);
    }
  }

  return (@row) ? 1 : 0;
}


=head2 store_status

  Arg [1]    : string - status e.g. IMPORTED, DISPLAYABLE
  Arg [2]    : Bio::EnsEMBL::"OBJECT"
  Example    : $status_a->store_status('IMPORTED', $array_chip);
  Description: Sets a state for a given object
  Returntype : None
  Exceptions : None
  Caller     : general
  Status     : At risk - Move to Status.pm?

=cut

sub store_status{
  my ($self, $state, $obj) = @_;

  my $sql;
  my $table = $self->_test_funcgen_table($obj);

  if(! $self->has_stored_status($state, $obj)){
    my $status_id = $self->_get_status_name_id($state);

    if(! $status_id){
      throw("$state is not a valid status_name for $obj:\t".$obj->dbID);
    }

	$sql = "INSERT into status(table_id, table_name, status_name_id) VALUES('".$obj->dbID()."', '$table', '$status_id')";
	$self->db->dbc->do($sql);

	#Setting it in the obj if it is not already present.
	$obj->add_status($state) if(! $obj->has_status($state, $obj));
  }

  return;
}


=head2 revoke_status

  Arg [1]    : string - status name e.g. 'IMPORTED'
  Arg [2]    : Bio::EnsEMBL::Funcgen::Storable
  Example    : $rset_adaptor->revoke_status('DAS DISPLAYABLE', $result_set);
  Description: Revokes the given state of Storable in status table.
  Returntype : None
  Exceptions : Warns if storable does not have state
               Throws is status name is not valid
               Throws if not state passed
  Caller     : General
  Status     : At Risk

=cut

sub revoke_status{
  my ($self, $state, $storable) = @_;

  throw('Must provide a status name') if(! defined $state);
  my $table_name = $self->_test_funcgen_table($storable);
  my $status_id = $self->_get_status_name_id($state);

  if ($status_id){

    #hardcode for InputSubset
    $table_name = 'input_subset' if $storable->isa('Bio::Ensembl::Funcgen:InputSubset');


    if(! $self->has_stored_status($state, $storable)){
  	warn $storable.' '.$storable->dbID()." does not have status $state to revoke\n";
  	return;
    }

    #do sanity checks on table to ensure that IMPORTED does not get revoke before data deleted?
    #how would we test this easily?

    my $sql = "delete from status where table_name='${table_name}'".
  	 " and status_name_id=${status_id} and table_id=".$storable->dbID();

    $self->db->dbc->db_handle->do($sql);

    #now splice from status array;
    #splice in loop should work as we will only see 1
    #Just hash this?

    for my $i(0..$#{$storable->{'states'}}){

  	  if($storable->{'states'}->[0] eq $state){
  	    splice @{$storable->{'states'}}, $i, 1;
  	    last;
  	  }
    }
  }
  #throw here?

  return;
}


=head2 revoke_states

  Arg [1]    : Bio::EnsEMBL::Funcgen::Storable
  Example    : $rset_adaptor->revoke_status($result_set);
  Description: Deletes all status records associated with the passed Storable.
  Returntype : Bio::EnsEMBL::Funcgen::Storable
  Exceptions : None
  Caller     : General + Helper rollback methods
  Status     : At Risk

=cut

sub revoke_states{
  my ($self, $storable) = @_;

  my $table_name = $self->_test_funcgen_table($storable);
  #add support for InputSubset which doesn't currently have an adaptor
  $table_name = 'input_subset' if $storable->isa('Bio::Ensembl::Funcgen::InputSubset');
  my $sql = "delete from status where table_name='${table_name}'".
	" and table_id=".$storable->dbID();

  #Do delete and clear stored states
  $self->db->dbc->db_handle->do($sql);
  undef $storable->{'states'};
  return $storable;
}


=head2 set_imported_states_by_Set

  Arg [1]    : Bio::EnsEMBL::Funcgen::Set e.g. a FeatureSet or ResultSet
  Example    : $self->set_imported_states_by_Set($set);
  Description: Sets default states for imported Feature|ResultSets
  Returntype : None
  Exceptions : None
  Caller     : Import parsers and RunnableDBs
  Status     : At risk - move to BaseImporter

=cut

#This needs to be used by RunnableDBs too!
#All state stuff is handled by BaseAdaptor?
#Can we put this in the SetAdaptor?

sub set_imported_states_by_Set{
  my ($self, $set) = @_;

  $self->db->is_stored_and_valid('Bio::EnsEMBL::Funcgen::Set', $set);
  #This should really be restricted to FeatureSet and ResultSet

  #Store default states for FeatureSets
  #DAS_DISPLAYABLE IMPORTED_'CSVERSION'
  #These need to insert ignore as may already be present?
  #Insert ignore may not catch an invalid status
  #So add states and store states as this checks
  $set->adaptor->store_status('DAS_DISPLAYABLE', $set);


  #To get assembly version here we need to
  # 1 get the current default chromosome version
  # or/and
  # 2 Use the assembly param to guess it from the coord_sys table
  # #This may pose problems for DB names which use numbers in their genebuild version
  # Would need to set this as a flag and/or specify the genebuild version too
  # Currently dnadb is set to last dnadb with 'assembly' as default version
  # We should match as test, just to make sure

  #Get default chromosome version for this dnadb
  my $cs_version = $self->db->dnadb->get_CoordSystemAdaptor->fetch_by_name('chromosome')->version;

  #Sanity check that it matches the assembly param?
  #Woould need to do this if ever we loaded on a non-default cs version

  $set->adaptor->store_status("IMPORTED_${cs_version}", $set);
}


=head2 status_filter

  Arg [1]    : string - status e.g. IMPORTED, DISPLAYABLE
  Arg [2]    : string - table name e.g. experimental_chip
  Arg [3]    : list   - table dbIDs
  Exmaple    : my @displayable_ec_ids = @{$ec_adaptor->status_filter('DISPLAYABLE',
                                                                     'experimental_chip',
                                                                     (map $_->dbID, @echips))};
  Description: Quick method for filtering dbIDs based on their table and and status
  Returntype : ARRAYREF
  Exceptions : Warns if state already set
               Throws is status name is not already stored.
  Caller     : general - ResultSetAdaptor->fetch_Resultfeatures_by_Slice
  Status     : At risk - Move to Status?

=cut

sub status_filter{
  my ($self, $status, $table_name, @table_ids) = @_;

  my @status_ids;
  my $status_id = $self->_get_status_name_id($status);

  if($status_id){

    throw("Must provide a table_name and table_ids to filter non-displayable ids")
      if(! ($table_name && @table_ids));

    my $sql = "SELECT table_id from status where table_name='$table_name' and ".
               "table_id in (".join(", ", @table_ids).") and status.status_name_id".
               "='$status_id'";
    @status_ids = map $_ = "@$_",
                   @{$self->db->dbc->db_handle->selectall_arrayref($sql)};
  }

  return \@status_ids;
}


=head2 _get_status_name_id

  Arg [1]    : String - status_name e.g. IMPORTED, DISPLAYABLE
  Arg [2]    : Boolean - optional flag to throw error if status_name is not
               present in the DB.
  Example    : my $status_id = $self->_get_status_name_id('IMPORTED');
  Description: Retrieves the dbID of a given status_name record
  Returntype : Int
  Exceptions : Throws if status name argument not defined or if throw flag is
               set and status_name is not in the DB.
  Caller     : Bio::EnsEMBL::Funcgen::BaseAdaptor
  Status     : At risk - Move to Status?

=cut

sub _get_status_name_id{
  my ($self, $status, $throw) = @_;

  if(! defined $status){
    throw('You must provide a status_name string argument');
  }

  my $sql = "SELECT status_name_id from status_name where name='$status'";
  my ($status_id) = $self->db->dbc->db_handle->selectrow_array($sql);

  if (! $status_id){
    if($throw){
      throw("Status name $status is not valid. ".
        'Maybe you need to add it to the status_name table?');
    }
    else{
      warn("Status name $status is not valid. ".
        'Maybe you need to add it to the status_name table?');
    }
  }

  return $status_id;
}


=head2 fetch_all_by_external_name

  Arg [1]    : String $external_name
               An external identifier of the feature to be obtained
  Arg [2]    : (optional) String $external_db_name
               The name of the external database from which the
               identifier originates.
  Example    : my @features =
                  @{ $adaptor->fetch_all_by_external_name( 'NP_065811.1') };
  Description: Retrieves all features which are associated with
               an external identifier such as an Ensembl Gene or Transcript
               stable ID etc.  Usually there will only be a single
               feature returned in the list reference, but not
               always.  Features are returned in their native
               coordinate system, i.e. the coordinate system in which
               they are stored in the database.  If they are required
               in another coordinate system the Feature::transfer or
               Feature::transform method can be used to convert them.
               If no features with the external identifier are found,
               a reference to an empty list is returned.
  Returntype : arrayref of Bio::EnsEMBL::Funcgen::Storable objects
               Maybe any Feature, FeatureType, Probe or ProbeSet
  Exceptions : Warns if method not available for given object adaptor
  Caller     : general
  Status     : at risk

=cut

sub fetch_all_by_external_name {
  my ( $self, $external_name, $external_db_name ) = @_;

  my $entryAdaptor = $self->db->get_DBEntryAdaptor();
  my (@ids, $type, $type_name);
  ($type = ref($self)) =~ s/.*:://;
  $type =~ s/Adaptor$//;
  ($type_name = $type) =~ s/Feature$/_feature/;
  my $xref_method = 'list_'.lc($type_name).'_ids_by_extid';

  if(! $entryAdaptor->can($xref_method)){
	warn "Does not yet accomodate $type external names";
	return [];
  }

  #Would be better if _list_ method returned and arrayref
  return $self->fetch_all_by_dbID_list([$entryAdaptor->$xref_method($external_name, $external_db_name)]);
}


=head2 fetch_all_by_external_names

  Arg [1]    : ARRAYREF of strings. External identifiers of the features to be obtained
  Arg [2]    : (optional) String $external_db_name
               The name of the external database from which the
               identifier originates.
  Example    : my @features =
                  @{ $adaptor->fetch_all_by_external_names(['ENST00003548913', ...])};
  Description: Retrieves all features which are associated with
               the external identifiers such as a Ensembl gene or transcript
               stable IDs, etc.  Features are returned in their native
               coordinate system, i.e. the coordinate system in which
               they are stored in the database.  If they are required
               in another coordinate system the Feature::transfer or
               Feature::transform method can be used to convert them.
               If no features with the external identifier are found,
               a reference to an empty list is returned.
  Returntype : arrayref of Bio::EnsEMBL::Funcgen::Storable objects
               Maybe any Feature, FeatureType, Probe or ProbeSet
  Exceptions : Warns if xref method not available for given object adaptor
  Caller     : general
  Status     : at risk

=cut

sub fetch_all_by_external_names{
  my ( $self, $external_names, $external_db_name ) = @_;

  my $entryAdaptor = $self->db->get_DBEntryAdaptor();
  my ($type, $type_name);
  ($type = ref($self)) =~ s/.*:://;
  $type =~ s/Adaptor$//;
  ($type_name = $type) =~ s/Feature$/_feature/;
  my $xref_method = 'list_'.lc($type_name).'_ids_by_extids';


  if(! $entryAdaptor->can($xref_method)){
	warn "Does not yet accomodate $type external names";
	return [];
  }

  #Would be better if _list_ method returned and arrayref
  my @ids = $entryAdaptor->$xref_method($external_names, $external_db_name);

  return $self->fetch_all_by_dbID_list(\@ids);
}


=head2 fetch_all_by_linked_Transcript

  Arg [1]    : Bio::EnsEMBL::Transcript
  Example    : my @psets =
                  @{ $probe_set_adaptor->fetch_all_by_linked_Transcript($tx_obj);
  Description: Retrieves all features which are associated with
               the given Ensembl Transcript.
  Returntype : ARRAYREF of Bio::EnsEMBL::Funcgen::Storable objects
               Maybe any Feature, FeatureType, Probe or ProbeSet
  Exceptions : Throws if arguments not valid
  Caller     : general
  Status     : at risk

=cut

sub fetch_all_by_linked_Transcript{
  my ($self, $tx) = @_;

  if(! $tx ||
	 ! (ref($tx) && $tx->isa('Bio::EnsEMBL::Transcript') && $tx->dbID)){
	throw('You must provide a valid stored Bio::EnsEMBL:Transcript object');
  }

  return $self->fetch_all_by_external_name($tx->stable_id, $self->db->species.'_core_Transcript')
}

=head2 fetch_all_by_linked_transcript_Gene

  Arg [1]    : Bio::EnsEMBL::Gene
  Example    : my @psets =
                  @{ $probe_set_adaptor->fetch_all_by_linked_transcript_Gene($gene_obj);
  Description: Retrieves all features which are indirectly associated with
               the given Ensembl Gene, through it's Transcripts.
  Returntype : ARRAYREF of Bio::EnsEMBL::Funcgen::Storable objects
               Maybe any Feature, FeatureType, Probe or ProbeSet
  Exceptions : Throws if arguments not valid
  Caller     : general
  Status     : at risk

=cut


sub fetch_all_by_linked_transcript_Gene{
   my ( $self, $gene ) = @_;

   if(! $gene ||
	  ! (ref($gene) && $gene->isa('Bio::EnsEMBL::Gene') && $gene->dbID)){
	 throw('You must provide a valid stored Bio::EnsEMBL:Gene object');
   }
   #No need to quote param here as this is a known int from the DB.
   my $tx_sids = $gene->adaptor->db->dbc->db_handle->selectcol_arrayref('select tsid.stable_id from transcript_stable_id tsid, transcript t where t.gene_id='.$gene->dbID.' and t.transcript_id=tsid.transcript_id');

   return $self->fetch_all_by_external_names($tx_sids, $self->db->species.'_core_Transcript');
}



=head2 store_associated_feature_types

  Arg [1]     : Bio::EnsEMBL::Funcgen::Storable
  Example     : $ext_feat_adaptor->store_associated_feature_type($ext_feat);
  Description : Store FeatureTypes assoicated with a given Storable
  Returntype  : None
  Exceptions  : Throws if FeatureTypes are not valid or stored
  Caller      : Adaptors
  Status      : At risk

=cut

#This currently does not validation
#of existing types
#we may want to rollback the associated ftypes from a specific analysis
#is rollback function safe here
#do we want to split this into update/add and store methods?


sub store_associated_feature_types {
  my ($self, $storable, $rollback) = @_;

  #Direct access to avoid lazy loading with an unstored SetFeature
  my $assoc_ftypes = $storable->{'associated_feature_types'};

  #Could be undef or empty
  return if ! defined $assoc_ftypes || scalar(@$assoc_ftypes) == 0;

  my $table_name = $storable->adaptor->_main_table->[0];
  my $dbid = $storable->dbID;
  my $sql;

  if($rollback){
    $sql = 'DELETE from associated_feature_table where '.
      "table_name='$table_name' and table_id=$dbid";
    $self->db->dbc->db_handle->do($sql);
  }

  $sql = 'INSERT into associated_feature_type(table_id, table_name, feature_type_id) values (?,?,?)';

  foreach my $ftype(@$assoc_ftypes){

	#We have already tested the class but not whether it is stored
	$self->db->is_stored_and_valid('Bio::EnsEMBL::Funcgen::FeatureType', $ftype);

	my $sth = $self->prepare($sql);
	$sth->bind_param(1, $dbid,        SQL_INTEGER);
	$sth->bind_param(2, $table_name,  SQL_VARCHAR);
	$sth->bind_param(3, $ftype->dbID, SQL_INTEGER);
	$sth->execute();
  }

  return;
}


=head2 fetch_all_by_associated_FeatureType

  Arg [1]    : Bio::EnsEMBL::Funcgen::FeatureType
  Example    : my $assoc_ftypes = $ft_adaptor->fetch_all_by_associated_FeatureType($ftype);
  Description: Fetches all objects which have associated FeatureType.
               Note this is not the main FeatureType for this object.
               Likely only relevant for FeatureTypeAdaptor and ExternalFeatureAdaptor.
  Returntype : ARRAYREF of Bio::EnsEMBL::Funcgen::Storable objects
  Exceptions : Throws if FeatureType not valid or stored
  Caller     : General
  Status     : At risk

=cut

#todo add _constrain_associated_FeatureType method

sub fetch_all_by_associated_FeatureType{
  my ($self, $ftype) = @_;

  $self->db->is_stored_and_valid('Bio::EnsEMBL::Funcgen::FeatureType', $ftype);
  my ($table_name, $table_syn) = @{$self->_main_table};

  $self->_tables([['associated_feature_type', 'aft']]);
  my $constraint = "aft.feature_type_id=? AND aft.table_name='${table_name}' ".
    "AND aft.table_id=${table_syn}.${table_name}_id";

  $self->bind_param_generic_fetch($ftype->dbID,  SQL_INTEGER);
  my $objs = $self->generic_fetch($constraint);
  $self->reset_true_tables;

  return $objs;
}


=head2 _list_dbIDs

  Example    : my @table_ids = @{$adaptor->_list_dbIDs()};
  Description: Wrapper for parent method, dynamically passes correct table name to query.
               Gets an array of internal IDs for all objects in the main table of this class.
  Returntype : List of Ints
  Exceptions : None
  Caller     : General
  Status     : Stable

=cut

sub _list_dbIDs{
  my $self = shift;
  return $self->SUPER::_list_dbIDs($self->_table_name);
}


=head2 build_feature_class_name

  Arg[1]     : String - feature_class e.g. annotated, dna_methylation, regulatory etc.
  Example    : my $fclass_name = $adaptor->build_feature_class_name;
  Description: Builds the full feature class name for a given feature class.
  Returntype : String
  Exceptions : None
  Caller     : FeatureSet::get_FeatureAdaptor and Set::feature_class_name
  Status     : Stable

=cut

sub build_feature_class_name{
  my ($self, $fclass) = @_;

  if(! defined $fclass){
    throw('You must pass a feature class argument to build the feature class name');
  }

  my @words = split('_', $fclass);

  foreach my $word(@words){
    $word = ucfirst($word);
    $word = 'DNA' if $word eq 'Dna';
  }

  return join('', (@words, 'Feature') );
}



sub stable_id_prefix{
  my $self = shift;

  if (! defined $self->{'stable_id_prefix'}) {
    ($self->{'stable_id_prefix'}) = @{$self->db->dnadb->get_MetaContainer->list_value_by_key
                                        (
                                         'species.stable_id_prefix'
                                        )};

    #Only add R if it is defined
    $self->{'stable_id_prefix'} .= 'R' if $self->{'stable_id_prefix'};
  }

  return $self->{'stable_id_prefix'}
}



### GENERIC CONSTRAINT METHODS ###
#   Used by compose_constraint_query

=head2 _constraint_states

  Arg[1]     : Arrayref - Status name strings.
  Arg[2]     : Hashref  - Constraint parameters e.g.
                {string_param_exists => 1} # This validates the status name exists
  Example    : my ($constraint, $constraint_conf) =
                 $self->$constrain_method($params->{constraints}{$constraint_key}, 
                                          $params);
  Description: Builds the constraint sql for resitricting to Storables with all the
               status names specified (AND logic).
  Returntype : List constaining a constraint sql string, and a hashref of constraint config.
  Exceptions : Throws if args are not valid.
  Caller     : BaseAdaptor::compose_constraint_query
  Status     : At risk 

=cut

#This uses AND logic, rather than the OR logic of the other constrain methods
#todo allow OR logic via a different method?
#todo move string_params_exists to caller or just use directly?
#and never assume if exists is true?

sub _constrain_states {
  my ($self, $states, $params) = @_;

  if(! (defined $states &&
        (ref($states) eq 'ARRAY') &&
        (scalar(@$states) > 0) )){
    throw('Must pass an Arrayref of states (strings) to contrain by');
  }

  if(defined $params &&
     (ref($params) ne 'HASH') ){
    throw('Params argument must be a Hashref');     
  }


  my @tables = $self->_tables;
  my ($table_name, $syn) = @{$self->_main_table};
  my ($status_table, $sn_ids_clause);


  my @sn_ids = sort {$a<=>$b}
                (map $self->_get_status_name_id($_, $params->{string_param_exists}) ||
                 'NULL', @$states);
  #|| NULL here accounts for absent status_names
  #i.e. $sn_ids_clause will never be true

  if(scalar(@$states) != 1){
    #add in table_name to make it faster
    #can't put in table_id as this would be a join between select and subselect
    $status_table = '(SELECT table_id, table_name, group_concat(status_name_id) ids '.
                    'FROM status WHERE table_name="'.$table_name.'" and ('.
                    join(' OR ', (map "status_name_id=$_", @sn_ids)).
                    ') group by table_id order by status_name_id)';

    #This enforces AND logic, whilst allowing for records with a superset of states
    $sn_ids_clause = ' s.ids like "%'.join('%,', @sn_ids).'%"';
  }
  else{
    $status_table  = 'status';
    $sn_ids_clause = 's.status_name_id='.$sn_ids[0];
  }


  my $constraint_conf = { tables => [[$status_table, 's']]};  #,['status_name', 'sn']),

  my $constraint = " $syn.${table_name}_id=s.table_id AND ".
    "s.table_name='$table_name' AND ".$sn_ids_clause;

  return ($constraint, $constraint_conf);
}



### DEPRECATED ###

sub list_dbIDs { #Deprecated in v69
	my $self = shift;
  deprecate('Please use _list_dbIDs.');
	return $self->_list_dbIDs($self->_main_table->[0]);
}

sub _constrain_status { #Deprecated in v73
  my ($self, $state) = @_;

  deprecate("The 'state' contraint key is deprecated, please use the following instead:\n\t".
    "'states' => ['state1', ...]");

  return $self->_constrain_states([$state]);
}

sub fetch_all_by_status{ #deprecated in v51
  my ($self, $status) = @_;

  deprecate('Use fetch_all($status) instead');
  return $self->fetch_all($status);

}


1;

