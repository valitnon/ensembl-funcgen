use strict;
use Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor;
 
Bio::EnsEMBL::DBSQL::DBAdaptor->new(
  -species => 'homo_sapiens',
  -group   => 'core',
  -host    => 'mysql-ens-sta-1.ebi.ac.uk',
  -port    => 4519,
  -user    => 'ensro',
  -dbname  => 'homo_sapiens_funcgen_91_38',
);
 
  
Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor->new(
  -species => 'homo_sapiens',
  -group   => 'funcgen',
  -host    => 'mysql-ens-reg-prod-1.ebi.ac.uk',
  -port    => 4526,
  -user    => 'ensadmin',
  -pass    => 'ensembl',
  -dbname  => 'jcmarca_homo_sapiens_funcgen_91_38',
);
 
# Avoids unnecessary open connections. Important when running on large numbers of species and the registry connects to many databases upon startup.
Bio::EnsEMBL::Registry->disconnect_all();
 
1;