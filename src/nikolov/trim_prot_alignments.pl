#!/usr/bin/env perl

use warnings;
use strict;
use Bio::SeqIO;
use Bio::Seq;
use Bio::LocatableSeq;
use Bio::SimpleAlign;
use Bio::AlignIO;

# percetage length present
my $percent_length=0.35;
# Xs present cutfoff
my $nx_percent_thresh=0.20;

if($#ARGV != 1) {
 print "Usage: ./trim_prot_alignment.pl input_dir output_dir\n";
 die;
}

my $input_dir=$ARGV[0];
my $output_dir=$ARGV[1];

print "Input directory: $input_dir\n";
print "Output directory: $output_dir\n";

print "Percent length: $percent_length\n";
print "X present cutoff: $nx_percent_thresh\n";

# total number of stop codons
my $nstop_tot=0;
# exon coverage array
my %exons=();

my $maxexons=0;

open(MYS,"species.txt");
 my @species=<MYS>;
close(MYS);

my @alignments=<$input_dir/trim_*NT.fasta>;

open(MYSTATS,">$output_dir/trim_macse.out");

foreach my $falign (@alignments) {
  chomp($falign);

  $maxexons++;

  open(MYAL,$falign);
  my @lines=<MYAL>;
  close(MYAL);

  open(MYAL,">tmp.fasta");
  foreach my $line (@lines) {
   chomp($line);
   if($line =~ m/>/) {
    my ($spec,$region) = split /\//,$line;
    $line = $spec;
   } else {
    $line =~ s/\!/N/g;
   };
   print MYAL "$line\n";
  }
  close(MYAL);
 
  my @ws = split(/\//,$falign);
  my @w = split(/_/,$ws[$#ws]);
  my $exon_name = $w[1];
  print "Exon name: $exon_name\n";

  # initialize the exon list
  foreach my $sp (@species){
   chomp($sp);
   $exons{$sp}{$exon_name}=0.0;
  }

  # read the alignment and create the alignment object
  my $str = Bio::AlignIO->new(-file => "tmp.fasta" , -format => 'fasta');
  my $aln = $str->next_aln();

  # length of alignment
  my $al = $aln->length;
  # minimum number of bases in alignment
  my $min_nnucl = $percent_length * $al;
  # seq counter
  my $nseq = 0;

  # count nucleotides in alignment and remove if less than threshold
  foreach my $seq ($aln->each_seq) {
    $nseq++;
    my $nnucl=0;
    for (my $i = 1; $i <= $al; $i++) {
      my $c = $seq->subseq($i,$i);
      if ($c =~ m/[aAtTcCgG]/) {
        $nnucl++;
      };
    };
    if($nnucl <= $min_nnucl) {
      my $did = $seq->display_id();
      print "Removing $did from alignment.\n";
      print "Too few nucleotides - $nnucl out of $al.\n";
      $aln->remove_seq($seq); 
    };
  };

  $nseq=0;

  # substitute gaps for Ns (including initial and tailing) and translate sequence
  my $paln = Bio::SimpleAlign->new();
  foreach my $seq ($aln->each_seq) {

    $nseq++;
    my $ncseq = $seq->seq;
    $ncseq =~ s/-/N/g;
    $seq->seq($ncseq);
    my $did = $seq->display_id();
    my $pseq = $seq->translate(-frame => 0);
    my $pcseq = $pseq->seq;
    my $pcseq_leng = $pseq->length;

    my $plseq = Bio::LocatableSeq->new(-seq => $pcseq, -id  => $did, -start => 1, -end   => $pcseq_leng);
    $paln->add_seq($plseq);
  };

  print "Number of sequences: $nseq\n";

  # evaluate protein alignment 
  my @del_sites=();

  my $paln_leng = $paln->length;

  my $nx_thresh = $nx_percent_thresh * $nseq;

  my $nstop_exon=0;

  for (my $pos=1; $pos <= $paln_leng; $pos++) {
    my $nx=0;
    my $nstop=0;
    foreach my $seq ($paln->each_seq) {
      my $res = $seq->subseq($pos,$pos);
      if ($res =~ m/[X]/) {
        $nx++;
      };
#     last if ($nx > $nx_thresh);
      if ($res =~ m/\*/) {
        $nstop++;
#       last;
      }
    };
    # erase codon column if number of ambiguous AA greater than threshold and/or there is a stop codon 
    if ($nx > $nx_thresh || ($nstop > 0 && $pos < $paln_leng)) {
      push(@del_sites,$pos);
    };
    if ($pos < $paln_leng) {
     $nstop_exon=$nstop_exon+$nstop;
    };
  };

  # trim the NT alignment
  my $trim_nts = ($#del_sites+1)*3;
  print "Trimming $trim_nts sites: ";
  for (my $i=$#del_sites; $i>=0; $i--) {
    my $pos = $del_sites[$i];
    my $ntstart=3*$pos-3;
    my $ntend=3*$pos-1;
    print "[$ntstart,$ntend] ";
    $aln = $aln->remove_columns([$ntstart,$ntend]);
  };
  print "\n";

  $nseq=0;
  $al=$aln->length;
  my $nts=0;

  # translate the trimmed alignment
  my $trim_paln = Bio::SimpleAlign->new();
  foreach my $seq ($aln->each_seq) {
    # collect statistics
    $nseq++;
    my $nnucl=0;
    for (my $i = 1; $i <= $al; $i++) {
      my $c = $seq->subseq($i,$i);
      if ($c =~ m/[aAtTcCgG]/) {
        $nnucl++;
      };
    }; 
    my $did = $seq->display_id();
    $exons{$did}{$exon_name}=$nnucl/$al;
    $nts=$nts+$nnucl;
    # translate
    my $pseq = $seq->translate(-frame => 0);
    my $pcseq = $pseq->seq;
    my $pcseq_leng = $pseq->length;
    my $plseq = Bio::LocatableSeq->new(-seq => $pcseq, -id  => $did, -start => 1, -end   => $pcseq_leng);
    $trim_paln->add_seq($plseq);
  };

  $nstop_tot=$nstop_tot+$nstop_exon;

  my $ntot = $nseq*$al;
  my $frac_nts = 100.*$nts/$ntot;
  print "Total number of characters: $ntot\n";
  print "Percent of nucleotides: $frac_nts\n";
  print "Number of stop codons: $nstop_exon\n";
  print MYSTATS "$exon_name\t$nseq\t$ntot\t$frac_nts\n";

  # write the trimmed NT alignment
  my $ftrim = $output_dir.'/trimal_'.$exon_name.'.fasta';
  my $nt_aln = Bio::AlignIO->new(-file => ">$ftrim", -format => 'fasta');
  $aln->sort_alphabetically;
  $nt_aln->write_aln($aln);

  # remove the region
  @lines=();
  open(MYAL,$ftrim);
   @lines=<MYAL>;
  close(MYAL);

  open(MYAL,">$ftrim");
  foreach my $line (@lines) {
   chomp($line);
   if($line =~ m/>/) {
    my ($name,$region) = split /\//,$line;
    $line = $name;
   };
   print MYAL "$line\n";
  }
  close(MYAL);

  # write the trimmed protein alignment
  #my $fprot = '>'.'prot_trimal/prot_'.$exon_name.'.fasta';
  #my $prot_aln = Bio::AlignIO->new(-file => $fprot, -format => 'fasta');
  #$trim_paln->sort_alphabetically;
  #$prot_aln->write_aln($trim_paln);
 
 }; 
close(MYSTATS);

print "Number of stop codons in the concatenated alignment: $nstop_tot\n";

open(MYOUT,">$output_dir/exon_coverage.txt");
  print MYOUT "   ";
  for (my $i=1;$i<=$maxexons;$i++) { print MYOUT "\t$i";};
  print MYOUT "\n";

  foreach my $spec (sort keys %exons) {
    print MYOUT "$spec";
    foreach my $ex (sort keys %{$exons{$spec}}) {
      my $excover = $exons{$spec}{$ex};
      if ($excover > 1.0) {
        print MYOUT "\t1.00";
      } else {
        my $outstr = sprintf("%.2f", $excover);
        print MYOUT "\t$outstr";
      };
    }
    print MYOUT "\n";
  }
close(MYOUT);
