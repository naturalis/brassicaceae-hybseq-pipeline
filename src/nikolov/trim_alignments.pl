#!/usr/bin/env perl

use warnings;
use strict;
use Bio::SeqIO;
use Bio::Seq;
use Bio::SimpleAlign;
use Bio::AlignIO;

my $pop_thresh=30;

if($#ARGV != 2) {
 print "Usage: ./trim_alignment.pl input_dir output_dir reference_sequence\n";
 die;
}

my $input_dir=$ARGV[0];
my $output_dir=$ARGV[1];
my $ref_seq=$ARGV[2];

print "Input directory: $input_dir\n";
print "Output directory: $output_dir\n";
print "Reference sequence: $ref_seq\n";

my @alignments=<$input_dir/align_*.fasta>;

foreach my $falign (@alignments) {
  chomp($falign);
  print "$falign\n";

  # read the alignment and create the alignment object
  my $str = Bio::AlignIO->new(-file => $falign , -format => 'fasta');
  my $aln = $str->next_aln();

  # length of alignment
  my $al = $aln->length;

  # start and end positions of references
  my $start=1;
  my $end=$al;
  print "$start $end\n";

  foreach my $seq ($aln->each_seq) {
    if($seq->display_id eq $ref_seq) {
      my $init=0;
      for (my $i = 1; $i <= $al; $i++) {
        my $res = $seq->subseq($i,$i);
        if ($res =~ m/-/ && $init == 0) {
          $start++;
        };
        if ($res =~ m/[aAtTcCgG]/) {
          $init=1;
        };
      };
      my $final=0;
      for (my $i = $al; $i >= 1; $i--) {
        my $res = $seq->subseq($i,$i);
        if ($res =~ m/-/ && $final == 0) {
          $end--;
        };
        if ($res =~ m/[aAtTcCgG]/) {
          $final=1;
        };
      };
    };
  }; 

  print "$start $end\n";

  # slice keeping start and end positions and gap-only
  my $trim_aln = $aln->slice($start,$end,1);

  # erase alignment site if the population is less than threshold amount of species
  my $trim_len = $trim_aln->length;

  my @del_sites=();

  for (my $pos=1; $pos <= $trim_len; $pos++) {
    my $pop=0;
    foreach my $seq ($trim_aln->each_seq) {
      my $res = $seq->subseq($pos,$pos);
      if ($res =~ m/[aAtTcCgG]/) {
        $pop++;
      };
      last if ($pop > $pop_thresh);
    };
    if ($pop <= $pop_thresh) {
      push(@del_sites,$pos);
    };
  };

  my $trim_n = $#del_sites+1;
  print "Trimming $trim_n sites: ";
  for (my $i=$#del_sites; $i>=0; $i--) {
    my $pos = $del_sites[$i];
    print "$pos ";
    $trim_aln = $trim_aln->remove_columns([$pos-1,$pos-1]);
  };
  print "\n";

  # write the trimmed alignment
  my @w = split /\//, $falign;
  my ($first,$froot) = split /_/,$w[$#w];
  my $ftrim = '>'.$output_dir.'/trim_'.$froot;
  print "$ftrim\n";
  my $trim = Bio::AlignIO->new(-file => $ftrim, -format => 'fasta');
  $trim->write_aln($trim_aln);
};

