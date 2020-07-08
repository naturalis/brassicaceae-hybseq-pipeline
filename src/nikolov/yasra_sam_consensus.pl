#!/usr/bin/env perl

use warnings;
use strict;

sub read_var_seq($) {
  my $fvar = $_[0];
  my $seq = "";
  my $nsnps=0;
  my $avcov=0;
  open(MYVAR,$fvar);
  foreach my $line (<MYVAR>) {
   chomp($line);
   if(index($line,"Position") < 0) {
     my ($chrom,$pos,$ref,$cons,$read1,$read2,$varfreq,$strands1,$strands2,$qual1,$qual2,$pval,$mapqual1,$mapqual2,$r1plus,$r1minus,$r2plus,$r2minus,$varallel) = split(m/\s+/,$line);
     if($cons =~ m/(.).*/) {
       $seq = $seq . $1;
       if($1 =~ m/[^N]/) {
        $nsnps++;
        $avcov+=$read2;
       }
     };
   }
  }
  if($nsnps > 0) {
   $avcov=$avcov/$nsnps;
  };
  close(MYVAR);
  return ($seq,$nsnps,$avcov);
}

sub read_out_stats($) {
  my $fout = $_[0];
  my $stats = "";
  open(MYOUT,$fout);
  foreach my $line (<MYOUT>) {
   chomp($line);
   if(index($line,"ERROR") >= 0) {
     print "$line\n";
   }
   if(index($line,"positions were called") > 0) {
     if($line =~ m/([0-9]+)/) {
       $stats = $stats . "\t" . $1 . ' total';
     };
   }
   if(index($line,"called Reference") > 0) {
     if($line =~ m/([0-9]+)/) {
       $stats = $stats . "\t" . $1 . ' Refs';
     };
   }
   if(index($line,"called SNP") > 0) {
     if($line =~ m/([0-9]+)/) {
       $stats = $stats . "\t" . $1 . ' SNPs';
     };
   }
   if(index($line,"called indel") > 0) {
     if($line =~ m/([0-9]+)/) {
       $stats = $stats . "\t" . $1 . ' indels';
     };
   }
  }
  close(MYOUT);
  return $stats;
};

my $fsam=$ARGV[0];

my $contig_name="";
my $ncontig=0;
my $nreads=0;

open(MYSAM,"$fsam");
foreach my $line (<MYSAM>) {
 chomp($line);
 if(index($line,'@') < 0) {
  my ($qname,$flag,$contig,$pos,$mapq,$cigar,$rnext,$pnext,$tlen,$read,$phred) = split(" ",$line);
  if($contig ne $contig_name) {
    if($contig_name ne "") {
       close(MYOUTSAM);
       # make the new consensus sequence here
       my @comp = split(/_/,$contig_name);
       my $cmd = 'makeconsensus.sh ' . $comp[0] . '_' . $comp[1] . '.sam > out.log 2>&1';
       system("$cmd");

       my $stats = read_out_stats("out.log");

       print "$contig_name:\t$stats\n";

       my $fvar = $comp[0] . '_' . $comp[1] . '.var';
       my ($seq,$nsnps,$avcov) = read_var_seq($fvar);

       print "Stats_coverage $nsnps $avcov\n";

       my $seq_len = length($seq);

       if($seq_len > 0) {
        
        my $ffa = '>' . $contig_name . '.fasta';
        open(MYFASTA,$ffa);
        print MYFASTA ">$contig_name\n";
        print MYFASTA "$seq\n";
        close(MYFASTA);
       } else {
        print "Fasta file is not printed: $contig_name\n";
       };
    };
    $contig_name = $contig;
    $ncontig++;
    my @comp = split(/_/,$contig_name);
    my $len = $comp[$#comp]-$comp[$#comp-1]+1;

    my $fname = '>' . $comp[0] . '_' . $comp[1] . '.sam';
    open(MYOUTSAM,$fname);

    print "Processing: $contig_name ...\n";
    print MYOUTSAM "\@HD\tVN:1.3\n";
    print MYOUTSAM "\@SQ\tSN:$contig_name\tLN:$len\n";
  };
  print MYOUTSAM "$line\n";
  $nreads++;
 };
};
close(MYOUTSAM);
close(MYSAM);

# make the new consensus sequence here
  my @comp = split(/_/,$contig_name);
  my $cmd = 'makeconsensus.sh ' . $comp[0] . '_' . $comp[1] . '.sam > out.log 2>&1';
  system("$cmd");

  my $stats = read_out_stats("out.log");

  print "$contig_name:\t$stats\n";

  my $fvar = $comp[0] . '_' . $comp[1] . '.var';
  my ($seq,$nsnps,$avcov) = read_var_seq($fvar);

  print "Stats_coverage $nsnps $avcov\n";

  my $seq_len = length($seq);

  if($seq_len > 0) {

   my $ffa = '>' . $contig_name . '.fasta';
   open(MYFASTA,$ffa);
   print MYFASTA ">$contig_name\n";
   print MYFASTA "$seq\n";
   close(MYFASTA);
  } else {
   print "Fasta file is not printed: $contig_name\n";
  };

print "Number of contigs is: $ncontig\n";
print "Number of reads: $nreads\n";
