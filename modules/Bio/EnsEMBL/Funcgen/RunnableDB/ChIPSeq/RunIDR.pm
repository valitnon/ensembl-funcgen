
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

Bio::EnsEMBL::Funcgen::RunnableDB::ChIPSeq::RunIDR

=head1 DESCRIPTION

=cut

package Bio::EnsEMBL::Funcgen::RunnableDB::ChIPSeq::RunIDR;

use warnings;
use strict;
use base 'Bio::EnsEMBL::Hive::Process';

use Bio::EnsEMBL::Utils::Scalar                 qw( assert_ref ); 
use Bio::EnsEMBL::Utils::Exception              qw( throw );
use Bio::EnsEMBL::Funcgen::Sequencing::SeqTools qw( run_IDR );
use Data::Dumper;

sub run {
  my $self = shift;
  
  my $pair    = $self->param_required('pair');
  my $tempdir = $self->param_required('tempdir');
  
  # --------------------------------------------------------------------------
  # Compute thresholds
  
  my %peaks_per_swembl_output_file;
  foreach my $permissive_peak_call_result (@$pair) {
  
    my $swembl_file = $permissive_peak_call_result->{peak_file};
  
    my $cmd = "grep -vE '(#|(^Region[[:space:]]+Start))' $swembl_file | wc -l | awk '{print \$1}'";
    
    use Bio::EnsEMBL::Funcgen::Utils::EFGUtils qw( run_backtick_cmd run_system_cmd );
    my $num_peaks = run_backtick_cmd($cmd);
    
    $self->say_with_header("Found " . $num_peaks . " peaks in bed file ${swembl_file}.", 1);
    
    if ($num_peaks == 0) {
      $self->throw("Found bed file with 0 peaks. This will cause IDR to fail, please remove or fix:\n\t$swembl_file");
    }
    $peaks_per_swembl_output_file{$swembl_file} = $num_peaks;
  }
  
  my @all_peak_numbers = values %peaks_per_swembl_output_file;
  
  $self->say_with_header("Number of peaks found in each bed file:\n" . Dumper(\%peaks_per_swembl_output_file), 1);
  
  use List::Util qw( min );
  my $min_num_peaks = min( @all_peak_numbers );
  
  $self->say_with_header("The minumum number of peaks is: " . $min_num_peaks, 1);
  
  my $idr_threshold = $self->idr_threshold($min_num_peaks);
  my $number_of_peaks_considered = $self->number_of_peaks_considered($min_num_peaks);

  $self->say_with_header("min_num_peaks              = $min_num_peaks", 1);
  $self->say_with_header("idr_threshold              = $idr_threshold", 1);
  $self->say_with_header("number_of_peaks_considered = $number_of_peaks_considered", 1);
  
  # --------------------------------------------------------------------------
  # Convert swembl output to bed files.
  
  my @bed_files;
  foreach my $permissive_peak_call_result (@$pair) {
  
    #die(Dumper($permissive_peak_call_result));
  
    my $swembl_file = $permissive_peak_call_result->{peak_file};
  
    my $bed_file = $swembl_file . '.bed';
    my $clean_bed_file = $swembl_file . '.clean.bed';

    my $cmd = 'awk \'BEGIN {OFS="\t"} { if($0 !~ /^(#|(Region[[:space:]]+Start))/) {print $1,$2,$3,".",$7,".",$7,-1,-1,int($9-$1)} }\' '.
      "$swembl_file | sort -k 7nr,7nr | head -n $min_num_peaks | sort -k 1,2n > " . $bed_file;
    
    $self->say_with_header("cmd = $cmd", 1);
    run_system_cmd($cmd);
    
    $self->write_bed_file_without_swembl_issues(
      $bed_file, 
      $clean_bed_file
    );
    push @bed_files, $clean_bed_file;
  }
  

  # --------------------------------------------------------------------------
  # Run idr analysis
  
  my $idr_output_file = $self->param_required('idr_output_file');
  
  my $cmd = "idr --idr-threshold $idr_threshold --output-file $idr_output_file --plot --use-old-output-format --samples " . join(' ', @bed_files);
  $self->say_with_header("cmd = $cmd", 1);
  
  use Bio::EnsEMBL::Funcgen::Utils::GoodUtils qw( run_cmd );
  #my $output = run_cmd($cmd);
  
  use Capture::Tiny ':all';
  (
    my $stdout, 
    my $stderr, 
    my $failed
  ) = capture {
    system( $cmd );
  };
  
  my $insufficient_merged_peaks_error;
  
  if ($failed) {
    $insufficient_merged_peaks_error 
      = $stderr =~ /ValueError: Peak files must contain at least 20 peaks post-merge/;
  }

  if ($failed && ! $insufficient_merged_peaks_error) {
    $self->throw(
      "Failed running command:\n\n"
      . $cmd . "\n\n"
      . "Stdout:\n\n"
      . $stdout . "\n\n"
      . "Stderr:\n\n"
      . $stderr . "\n\n"
    );
  }

  if ($insufficient_merged_peaks_error) {
    $self->dataflow_output_id(
      {
        idr_result => {
          idr_num_peak_threshold    => undef,
          insufficient_merged_peaks => $insufficient_merged_peaks_error,
          peak_calling_pair         => $pair,
        }
      },
      2
    );
    return;
  }
  
  $cmd = "cat $idr_output_file | wc -l";
  my $num_peaks = run_cmd($cmd);
  chomp($num_peaks);
  
  $self->dataflow_output_id({
      idr_result => {
        idr_num_peak_threshold    => $num_peaks,
        insufficient_merged_peaks => undef,
        peak_calling_pair         => $pair,
      }
    }, 
    2
  );
  return;
}

sub write_bed_file_without_swembl_issues {

  my $self           = shift;
  my $bed_file       = shift;
  my $clean_bed_file = shift;
  
  open my $bed_fh, '<', $bed_file;
  open my $clean_bed_fh, '>', $clean_bed_file;
  
  my $print_new_line = undef;
  
  LINE:
  while (my $line = <$bed_fh>) {
    chomp($line);
    (
      my $region,
      my $start_pos,
      my $end_pos,
      my $count,
      my $length,
      my $unique_pos,
      my $score,
      my $ref_count,
      my $max_coverage,
      my $summit,
      my $pvalue
    ) = split "\t", $line;
    
    use Scalar::Util qw( looks_like_number );
    next LINE if (! looks_like_number($summit));
    
    #my $summit_ok = ($summit > $start_pos) && ($summit < $end_pos);
    my $summit_ok = $summit > 0;
    next LINE if (! $summit_ok);
    
    # Prevent flurry of warnings like:
    # Use of uninitialized value $pvalue in join or string at
    # ...
    no warnings 'uninitialized';
    
    if ($print_new_line) {
      $clean_bed_fh->print("\n");
    }
    
    $clean_bed_fh->print(
      join "\t",
        $region,
        $start_pos,
        $end_pos,
        $count,
        $length,
        $unique_pos,
        $score,
        $ref_count,
        $max_coverage,
        $summit,
        $pvalue
    );
    $print_new_line = 1;
  }
  close($clean_bed_fh);
  close($bed_fh);
  return;
}

sub idr_threshold {
  my $self = shift;
  my $min_num_peaks = shift;
  
  my $idr_threshold;
  
  my $number_of_peaks_considered = $self->number_of_peaks_considered($min_num_peaks);
  
  if ($number_of_peaks_considered < 100_000) {
    $idr_threshold = 0.05;
  } else {
    $idr_threshold = 0.01;
  }
  return $idr_threshold;
}

sub number_of_peaks_considered {

  my $self = shift;
  my $min_num_peaks = shift;

  my $number_of_peaks_considered = $self->max_peaks_for_this_peak_caller('swembl');

  if($min_num_peaks < $number_of_peaks_considered) {
    # We take the lowest number of peaks, as we need comparable numbers of peaks across all inputs
    $number_of_peaks_considered = $min_num_peaks;
  }
  return $number_of_peaks_considered;
}

sub max_peaks_for_this_peak_caller {

  my $self        = shift;
  my $peak_caller = shift;
  
  my $max_peaks_for_this_peak_caller;
  
  if ($peak_caller eq 'swembl') {
    $max_peaks_for_this_peak_caller = 300_000;
  }
  if ($peak_caller eq 'macs') {
  
    # Copied threshold over from the old pipeline where it was never
    # used.
    #
    $max_peaks_for_this_peak_caller = 100_000;
  }
  if (! defined $max_peaks_for_this_peak_caller) {
    throw("Unknown peak caller $peak_caller");
  }
  return $max_peaks_for_this_peak_caller;
}

1;


