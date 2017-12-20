package clsEpiRR;

    sub new {
	
        my ($class, $args) = @_;

    	my $self = {
        	project => $args->{project} || '',
        	species  => $args->{species} || '',
			accession  => $args->{accession} || '',
			ihec  => $args->{ihec} || '',
			epigenome  => $args->{epigenome} || '',
			refStage  => $args->{refStage} || '',
			ATAC_Seq  => $args->{ATAC_Seq} || '',
			ChipS_In  => $args->{ChipS_In} || '',
			DNase_Seq  => $args->{DNase_Seq} || '',
			H3K27ac  => $args->{H3K27ac} || '',
			H3K27me3  => $args->{H3K27me3} || '',
			H3K36me3  => $args->{H3K36me3} || '',
			H3K4me1  => $args->{H3K4me1} || '',
			H3K4me3  => $args->{H3K4me3} || '',
			H3K9me3  => $args->{H3K9me3} || '',
			CTCF  => $args->{CTCF} || '',
			avExperiments  => $args->{avExperiments} || '',
			others  => $args->{others} || ''
    	};
    	return bless $self, $class;

	}
	
	sub set_project {
		my ($self, $project) = @_;
	    $self->{project} = $project;	
	}
	
	sub set_species {
		my ($self, $species) = @_;
	    $self->{species} = $species;	
	}
	
	sub set_accession {
		my ($self, $accession) = @_;
	    $self->{accession} = $accession;	
	}
	
	sub set_ihec {
		my ($self, $ihec) = @_;
	    $self->{ihec} = $ihec;	
	}
	
	sub set_epigenome {
		my ($self, $epigenome) = @_;
	    $self->{epigenome} = $epigenome;	
	}
	
	sub set_referenceStage {
		my ($self, $refStage) = @_;
	    $self->{refStage} = $refStage;	
	}
	
	sub set_ATAC_Seq {
		my ($self, $ATAC_Seq) = @_;
	    $self->{ATAC_Seq} = $ATAC_Seq;	
	}
	
	sub set_ChipS_In {
		my ($self, $ChipS_In) = @_;
	    $self->{ChipS_In} = $ChipS_In;	
	}
	
	sub set_DNase_Seq {
		my ($self, $DNase_Seq) = @_;
	    $self->{DNase_Seq} = $DNase_Seq;	
	}
	
	sub set_H3K27ac {
		my ($self, $H3K27ac) = @_;
	    $self->{H3K27ac} = $H3K27ac;	
	}
	
	sub set_H3K27me3 {
		my ($self, $H3K27me3) = @_;
	    $self->{H3K27me3} = $H3K27me3;	
	}
	
	sub set_H3K36me3 {
		my ($self, $H3K36me3) = @_;
	    $self->{H3K36me3} = $H3K36me3;	
	}
	
	sub set_H3K4me1 {
		my ($self, $H3K4me1) = @_;
	    $self->{H3K4me1} = $H3K4me1;	
	}
	
	sub set_H3K4me3 {
		my ($self, $H3K4me3) = @_;
	    $self->{H3K4me3} = $H3K4me3;	
	}
	
	sub set_H3K9me3 {
		my ($self, $H3K9me3) = @_;
	    $self->{H3K9me3} = $H3K9me3;	
	}
	
	sub set_CTCF {
		my ($self, $CTCF) = @_;
	    $self->{CTCF} = $CTCF;	
	}
	
	sub set_avExperiments {
		my ($self, $avExperiments) = @_;
	    $self->{avExperiments} = $avExperiments;	
	}
	
	sub set_others {
		my ($self, $others) = @_;
	    $self->{others} = $others;	
	}
	
	sub csv_row {
		my $self = shift;
		my $csv_row = join ("\t", $self->{project}, $self->{species}, $self->{accession}, $self->{ihec}, $self->{epigenome}, $self->{refStage}, 
		$self->{ATAC_Seq}, $self->{ChipS_In}, $self->{DNase_Seq}, $self->{H3K27ac}, $self->{H3K27me3}, $self->{H3K36me3}, $self->{H3K4me1}, 
		$self->{H3K4me3}, $self->{H3K9me3}, $self->{CTCF}, $self->{avExperiments}, $self->{others});
		
		return $csv_row;
	}	
1;