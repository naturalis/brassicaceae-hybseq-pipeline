#!/usr/bin/env perl

use warnings;
use strict;
use Bio::SeqIO;
use Bio::Seq;
use Bio::SimpleAlign;
use Bio::AlignIO;
use Cwd;

################# COMMAND-LINE ARGUMENTS #################
# 1: Database Filename (2bit format)
# 2: Species name
# 3: Step: 1 - blatx, 2 - mafft, 3 - consensus
##########################################################

# At least half of the target and query have to be matched
my $tCutoff = 0.0;
my $qCutoff = 0.4;

my $blatoptions = "-t=dnax -q=dnax -stepSize=5 -repMatch=2253 -minScore=20 -minIdentity=0";
my $mafftoptions = "--maxiterate 1000 --genafpair";

#make an output directory
my $out_dir='final_exons';
if( !( -d $out_dir ) ) {
 print "$out_dir does not exist\n";
 mkdir $out_dir;
};

# make a blat dir
my $blat_dir = 'blats';
if( !( -d $blat_dir ) ) {
 print "$blat_dir does not exist\n";
 mkdir $blat_dir;
};

# make an alignment directory
my $mafft_dir = 'maffts';
if( !( -d $mafft_dir ) ) {
 print "$mafft_dir does not exist\n";
 mkdir $mafft_dir;
};

# current working directory
my $dir_root= cwd();

# identifiers of references
my @identifiers = ("at","aa","si");
my @identifiers_caps = ("AT","SI","AA");

#command-line arguments
my $file=$ARGV[0];
my $spec=$ARGV[1];
my $step=$ARGV[2];

## common data
# temporary variable
my @pslString = ();

#### Subroutines

# All file names include the absolute path to the file

# BLAT a file with a database and print the result in directory
# Use: DoBlat(path-to-database,seqfile,contig name,suffix,printing directory)
sub DoBlat($$$) {
  my ($database,$file,$fout) = @_;

  #blat the contig 
  my $comm=join " ",("blat",$blatoptions,$database,$file,$fout,' > blat.out 2>&1');

  return $comm;
  #system($comm);
};

sub DoMafft($$) {
  my ($filein,$fileout) = @_;

  #mafft filein 
  my $comm=join " ",('mafft',$mafftoptions,$filein,'>',$fileout,' 2> mafft.out');
  
  return $comm;
  #system($comm);
};

#########################################################################################################from pslScore.pl
sub pslIsProtein($$$$$$$) {
  my ($blockCount, $strand, $tStart, $tEnd, $tSize, $tStarts, $blockSizes) = @_;
  my @starts = split(',', $tStarts);
  my @sizes = split(',', $blockSizes);
  my $lastBlock = $blockCount - 1;
  my $answer = 1;
  if ($strand =~ m/^\+/) {
    my $test = $starts[$lastBlock] + (3 * $sizes[$lastBlock]);
    $answer = 3 if ($tEnd == $test);
  } elsif ($strand =~ m/^-/) {
    my $test = $tSize - ($starts[$lastBlock] + (3 * $sizes[$lastBlock]));
    $answer = 3 if ($tStart == $test);
  }
  return $answer;
} # pslIsProtein()

sub pslCalcMilliBad($$$$$$$$$$$) {
 my ($sizeMul, $qEnd, $qStart, $tEnd, $tStart, $qNumInsert, $tNumInsert, $matches, $repMatches, $misMatches, $isMrna) = @_;
 my $milliBad = 0;
 my $qAliSize = $sizeMul * ($qEnd - $qStart);
 my $tAliSize = $tEnd - $tStart;
 my $aliSize = $qAliSize;
 $aliSize = $tAliSize if ($tAliSize < $qAliSize);
 if ($aliSize <= 0) {
  return $milliBad;
 }
 my $sizeDif = $qAliSize - $tAliSize;
 if ($sizeDif < 0) {
  if ($isMrna) {
      $sizeDif = 0;
  } else {
      $sizeDif = -$sizeDif;
  }
 }
 my $insertFactor = $qNumInsert;
 if (0 == $isMrna) {
  $insertFactor += $tNumInsert;
 }
 my $total = ($sizeMul * ($matches + $repMatches + $misMatches));
 if ($total != 0) {
  my $roundAwayFromZero = 3*log(1+$sizeDif);
  if ($roundAwayFromZero < 0) {
    $roundAwayFromZero = int($roundAwayFromZero - 0.5);
  } else {
    $roundAwayFromZero = int($roundAwayFromZero + 0.5);
  }
  $milliBad = (1000 * ($misMatches*$sizeMul + $insertFactor + $roundAwayFromZero)) / $total;
 }
 return $milliBad;
} # sub pslCalcMilliBad()
#####################################################################################################end from pslScore.pl

# Provided the absolute path to the psl file get the contig and exon with highest score
# and write the contig name and exon name in pslString
sub ReadPsl($) {
 my ($fpsl) = @_;

 # read the psl file
 open(MYPSL,$fpsl);
  my @lpsl=<MYPSL>;
 close(MYPSL);

 # take the last line
 my $nline = $#lpsl;
 if($nline >= 5) {
  # check if multiple matches
  my $indmax=5;
  my $scoremax=0;
  for(my $ind=5; $ind<@lpsl; $ind++) {
   my $line=$lpsl[$ind];
   chomp($line);
   # parse the line
   my ($matches, $misMatches, $repMatches, $nCount, $qNumInsert, $qBaseInsert, $tNumInsert, $tBaseInsert, $strand, $qName, $qSize, $qStart, $qEnd, $tName, $tSize, $tStart, $tEnd, $blockCount, $blockSizes, $qStarts, $tStarts) = split('\t', $line);
   # compute score
   my $sizeMul = pslIsProtein($blockCount, $strand, $tStart, $tEnd, $tSize, $tStarts, $blockSizes);
   my $pslScore = $sizeMul * ($matches + ( $repMatches >> 1) ) - $sizeMul * $misMatches - $qNumInsert - $tNumInsert;
   my $milliBad = int(pslCalcMilliBad($sizeMul, $qEnd, $qStart, $tEnd, $tStart, $qNumInsert, $tNumInsert, $matches, $repMatches, $misMatches, 1));
   my $percentIdentity = 100.0 - $milliBad * 0.1;
   if($pslScore>$scoremax) {
    $scoremax=$pslScore;
    $indmax=$ind;
   };
  };
  my $line=$lpsl[$indmax];
  chomp($line);
  # parse the line
  my ($matches, $misMatches, $repMatches, $nCount, $qNumInsert, $qBaseInsert, $tNumInsert, $tBaseInsert, $strand, $qName, $qSize, $qStart, $qEnd, $tName, $tSize, $tStart, $tEnd, $blockCount, $blockSizes, $qStarts, $tStarts) = split('\t', $line);
  # compute score
  my $sizeMul = pslIsProtein($blockCount, $strand, $tStart, $tEnd, $tSize, $tStarts, $blockSizes);
  my $pslScore = $sizeMul * ($matches + ( $repMatches >> 1) ) - $sizeMul * $misMatches - $qNumInsert - $tNumInsert;
  my $milliBad = int(pslCalcMilliBad($sizeMul, $qEnd, $qStart, $tEnd, $tStart, $qNumInsert, $tNumInsert, $matches, $repMatches, $misMatches, 1));
  my $percentIdentity = 100.0 - $milliBad * 0.1;
  my $qRatio = $pslScore/$qSize;
  my $tRatio = $pslScore/$tSize;
  my $outstr = sprintf("%s\t%d\t%d\t%d\t%s\t%d\t%d\t%d\t%d\t%.2f\t%.3f\t%.3f\n", $tName, $tSize, $tStart, $tEnd, $qName, $qSize, $qStart, $qEnd, $pslScore, $percentIdentity,$tRatio,$qRatio);

  # add the contigs to the pslString
  @pslString = ();
  push(@pslString,$tName);
  push(@pslString,$qName);
  push(@pslString,$percentIdentity);
  push(@pslString,$pslScore);
  push(@pslString,$tRatio);
  push(@pslString,$qRatio);
  push(@pslString,$tStart);
  push(@pslString,$tEnd);
  push(@pslString,$tBaseInsert);
  push(@pslString,$qSize);
  push(@pslString,$strand);  

  return $outstr;
 } else {
  my $outstr = "  WARNING: No sequence match found :(\n";
  return $outstr;
 };
};

# make a consensus sequence out of a mafft alignment
sub MakeCons($$) {
my ($fmafft,$fadir)=@_;

# read the alignment and create the alignment object
my $str = Bio::AlignIO->new(-file => $fmafft , -format => 'fasta');
my $aln = $str->next_aln();

$aln->uppercase();
my $al = $aln->length;
my $nseq = 0;
# consensus sequence 
my %consensus = ();

foreach my $seq ($aln->each_seq) {
 for (my $i = 1; $i <= $al; $i++) {
  my $res = $seq->subseq($i,$i);
  if(exists $consensus{$i}{$res}) {
   $consensus{$i}{$res}++;   
  }else{
   $consensus{$i}{$res}=1;
  };
 };
};

for (my $i = 1; $i <= $al; $i++) {
 my %col = %{$consensus{$i}};
 my @key = sort keys %col;
 my $size = $#key + 1;
 # consensus criteria
 if($size == 1) {
  $consensus{$i}{'cons'}=$key[0];
 } elsif($size == 2 && (exists $col{'-'})) {
  $consensus{$i}{'cons'}=$key[1];
 } else {
  $consensus{$i}{'cons'}='N';
 }
}

my $cons_seq="";
for (my $i = 1; $i <= $al; $i++) {
 $cons_seq = $cons_seq . $consensus{$i}{'cons'};
}
# remove gaps
$cons_seq =~ tr/-//d;

# make a new file object
my @a = split /\//,$fmafft;
my @b = split /\.mafft\.fa/,$a[$#a];
my $fn = join("",(">",$fadir,'/',$b[0],'.fasta'));
my $out_file = Bio::SeqIO->new(-file => $fn, -format => 'Fasta', -alphabet => 'dna');

# make a new sequence object
my $did = $b[0];
my $out_seq  = Bio::Seq->new(-seq => $cons_seq, -display_id => $did, -alphabet => 'dna');

# write to file
$out_file->write_seq($out_seq);
};

sub FragmentOvlp($$$$$$) {
 my ($chr1,$s1,$e1,$chr2,$s2,$e2)=@_;

 if($chr1 eq $chr2) {
  my $ovl=0;
  if ($s1 >= $s2 && $s1 <= $e2) {$ovl++;};
  if ($e1 >= $s2 && $e1 <= $e2) {$ovl++;};
  if ($s2 >= $s1 && $s2 <= $e1) {$ovl++;};
  if ($e2 >= $s1 && $e2 <= $e1) {$ovl++;};
  if ($ovl > 0) {$ovl=1;};
  return $ovl;
 } else {
  return 0;
 }
};

#### End of subroutines section

# Switch to analysis after running the doblat.sh on the queue
  goto PSL if ($step == 2);
# Switch to consensus calling after mafft alignment
  goto CONS if ($step == 3);

BLAT: {
 print "Working directory:\n";
 print "$dir_root ...\n";
 chomp($file);
 print "Database (2bit format):\n";
 print "$file ...\n";

 # full path to database
 my $file_data = $dir_root.'/'.$file;

 my $fdatabase = '../'.$file;

 open(MYBLATSCR,">doblat.sh");
 print MYBLATSCR "#!/bin/bash\n";
 print MYBLATSCR "cd $blat_dir\n";

 foreach my $ref (@identifiers) {
  chomp($ref);
  my $ref_name = 'ref-'.$ref.'_split.fasta';
  print "Reference: $ref_name\n";  

  my $fin = '../'.$ref_name;
  my $in = Bio::SeqIO->new(-file => $fin, -format => 'Fasta', -alphabet => 'dna');

  while(my $rec=$in->next_seq()) {

    # make a new file object
    my $singleseq = $rec->seq;
    my $did = $rec->display_id();
    my $fname = '>'.$blat_dir.'/'.$did.'.fa';
    my $fname_full = $dir_root.'/'.$blat_dir.'/'.$did.'.fa';
    my $out_file = Bio::SeqIO->new(-file => $fname, -format => 'Fasta', -alphabet => 'dna');
    my $out_seq  = Bio::Seq->new(-seq => $singleseq, -display_id => $did, -alphabet => 'dna');
    $out_file->write_seq($out_seq);
 
    # BLAST to reference (first argument)
    my $fname_in = $did.'.fa';
    my $fpsl = $did.'.psl';
    my $str = DoBlat($fdatabase,$fname_in,$fpsl);

    print MYBLATSCR "$str\n";
  };
 }; 

 close(MYBLATSCR);
 exit 0;
};

PSL: {
 open(MYPSLOUT,">psl_stats.out");
 open(MYPSLERR,">psl_err.out");

 chomp($file);
 print "$file ...\n";

 my $file_data = $file;

 my %exon_target_all=();

 # repeat analysis foreach target reference
 foreach my $tRef (@identifiers_caps) {

 my $nexons=0;
 my $nexcov=0;
 my $nexf=0;

 my %exon_target=();
 my %exon_score=();
 my %exon_pid=();
 my %exon_tRatio=();
 my %exon_qRatio=();
 my %exon_fname=();
 my %exon_tStart=();
 my %exon_tEnd=();
 my %exon_strand=();

 my @psl=<$blat_dir/*.psl>;

 foreach my $fpsl (@psl) {
    chomp($fpsl);

    $nexons++;

    # read psl and find the highest score exon
    my $str = ReadPsl($fpsl);

    # check psl output and print stats
    if(index($str,'WARNING') > 0) {
      print MYPSLERR "$fpsl\n";
    } else {
   
      # analyze data
      my $refexon = $pslString[1];# qName (second argument)
      my $contig = $pslString[0]; # tName (first argument)
      my $pid = $pslString[2];    # percent identity
      my $score = $pslString[3];  # Blast score
      my $tRatio = $pslString[4]; # tSize/pslScore
      my $qRatio = $pslString[5]; # qSize/pslScore
      my $tStart = $pslString[6]; # Start nucleotide in target
      my $tEnd = $pslString[7];   # End nucleotide in target
      my $tBaseInsert = $pslString[8];   # #Inserted pairs in target
      my $qSize = $pslString[9];   # query size
      my $strand = $pslString[10]; # strandness query-target

      # analyze the list of existing exon(query)-contig(target) homolog pairs
      # keep only the pairs with highest score
      if($tRatio > $tCutoff && $qRatio > $qCutoff && $tBaseInsert < $qSize) {

       $nexf++;
       print MYPSLOUT "$str";

       my ($refname,$exon) = split /\_/,$refexon;

       if ($refname eq $tRef) {

        if(exists $exon_score{$exon}) {
         if( $score > $exon_score{$exon}) {
          $exon_target{$exon}=$contig; 
          $exon_score{$exon}=$score;
          $exon_pid{$exon}=$pid;
          $exon_tRatio{$exon}=$tRatio;
          $exon_qRatio{$exon}=$qRatio;
          $exon_fname{$exon}=$refexon;  
          $exon_tStart{$exon}=$tStart;
          $exon_tEnd{$exon}=$tEnd;    
          $exon_strand{$exon}=$strand;
         };
        } else {
          $exon_target{$exon}=$contig;
          $exon_score{$exon}=$score;
          $exon_pid{$exon}=$pid;
          $exon_tRatio{$exon}=$tRatio;
          $exon_qRatio{$exon}=$qRatio;
          $exon_fname{$exon}=$refexon;
          $exon_tStart{$exon}=$tStart;
          $exon_tEnd{$exon}=$tEnd;
          $exon_strand{$exon}=$strand;
        };

       };

     } else {
      print MYPSLERR "$fpsl $tRatio $qRatio DID NOT SATISFY CUTOFF\n";
     };
   };
 }; 

 # write sequences in maffts dir for alignment
 print "Printing exons in fasta files ...\n";
 foreach my $exon (keys %exon_target) {
  $nexcov++;
  $exon_target_all{$exon}{$tRef}=$exon_fname{$exon};
  
  # extract the exon from the genome
  my $tStart=$exon_tStart{$exon};
  my $tEnd=$exon_tEnd{$exon};
  my $tName=$exon_target{$exon};
  my $cmd = 'twoBitToFa -seq='.$tName.' -start='.$tStart.' -end='.$tEnd.' '.$file_data.' exon.fasta';
  system($cmd); 

  # change the identifier and append
  my $in = Bio::SeqIO->new(-file => 'exon.fasta', -format => 'Fasta', -alphabet => 'dna');
  my $rec = $in->next_seq();
  my $did = $tRef.'_'.$rec->display_id();
  # reverse complement the sequence if needed
  my $strand = $exon_strand{$exon};
  if($strand eq "+-" || $strand eq "-+") {
   $rec=$rec->revcom();
  };
  my $outf = '>>'.$mafft_dir.'/'.$exon.'.fa';
  my $out = Bio::SeqIO->new(-file => $outf, -format => 'Fasta', -alphabet => 'dna');
  my $outseq  = Bio::Seq->new(-seq => $rec->seq, -display_id => $did, -alphabet => 'dna');
  $out->write_seq($outseq);
 };

 # statistics for each reference
 print "Reference identifier: $tRef:\n";
 print "$nexons processed\n";
 print "$nexf found\n";
 print "$nexcov retained\n";

 };

 close(MYPSLOUT);
 close(MYPSLERR);

 # print out statistics 
 my $nexcov=0;
 my $flist = '>' . $out_dir . '/' . $spec . '_cov_ref.txt';
 open(MYLIST,$flist);
 foreach my $exon (sort keys %exon_target_all) {
  $nexcov++;
  my $outstr=$exon;
  foreach my $id (sort keys %{ $exon_target_all{$exon}}) {
   $outstr=$outstr.' '.$id;
  };
  print MYLIST "$outstr\n";
 };
 close(MYLIST);
 print "$nexcov found\n";

 # prepare a mafft script
 my @mafftf=<$mafft_dir/*.fa>;

 open(MYMAFFTSCR,">domafft.sh");
 print MYMAFFTSCR "#!/bin/bash\n";

 foreach my $ffasta (@mafftf) {
    # count number of sequences in the input file
    # check that the fragments are overlapping - throw away non-overlapping fragments and exons altogether, if not found
    my $nseq=0;
    my %did_name=();
    my %chr_name=();
    my %start_n=();
    my %end_n=();
    my $in = Bio::SeqIO->new(-file => $ffasta, -format => 'Fasta', -alphabet => 'dna');
    while(my $rec = $in->next_seq) { 
      $nseq++;
      # get the genomic region
      my $did = $rec->display_id();
      $did_name{$nseq}=$did;
      my @a = split /[\:-]/,$did;
      $chr_name{$nseq}=substr($a[0],3);
      $start_n{$nseq}=$a[1];
      $end_n{$nseq}=$a[2];
    };

    if($nseq>1) {
      # case 2 sequences
      if($nseq == 2) {
        # compute overlap
        my $ovl=FragmentOvlp($chr_name{1},$start_n{1},$end_n{1},$chr_name{2},$start_n{2},$end_n{2});
        # non-overlapping remove
        if($ovl == 0) {
         print "$ffasta contains non-overlapping seqs - removing\n";
         my $cmd = 'rm '.$ffasta;
         system($cmd);
         next;
        };
      };
      if($nseq == 3) {
        # compute overlaps
        my $o12 = FragmentOvlp($chr_name{1},$start_n{1},$end_n{1},$chr_name{2},$start_n{2},$end_n{2});
        my $o13 = FragmentOvlp($chr_name{1},$start_n{1},$end_n{1},$chr_name{3},$start_n{3},$end_n{3});
        my $o23 = FragmentOvlp($chr_name{3},$start_n{3},$end_n{3},$chr_name{2},$start_n{2},$end_n{2});
        # sum of overlaps
        my $sumo = $o12 + $o13 + $o23;

        # the three are non-overlapping - remove
        if($sumo == 0) {
         print "$ffasta contains non-overlapping seqs - removing\n";
         my $cmd = 'rm '.$ffasta;
         system($cmd);
         next;
        };

        if($sumo == 1) {
         # find non-overlapping
         my $did_no;
         if($o12 == 0 && $o13 == 0) {$did_no=$did_name{1};};
         if($o12 == 0 && $o23 == 0) {$did_no=$did_name{2};};
         if($o23 == 0 && $o13 == 0) {$did_no=$did_name{3};};

         print "$ffasta $did_no is non-overlapping - removing\n";

         # remove the sequence from the file
         my $nseq=0;
         my $in = Bio::SeqIO->new(-file => $ffasta, -format => 'Fasta', -alphabet => 'dna');
         my $out = Bio::SeqIO->new(-file => '>exons.fasta', -format => 'Fasta', -alphabet => 'dna');
         while(my $rec = $in->next_seq) {
          my $did = $rec->display_id();
          if(!($did eq $did_no)) {
            $nseq++;
            my $outseq  = Bio::Seq->new(-seq => $rec->seq, -display_id => $did, -alphabet => 'dna');
            $out->write_seq($outseq);
          };
         };
         my $cmd = 'mv exons.fasta '.$ffasta;
         system($cmd);
        };

      }
    };

    # prepare output name
    my @b=split /\.fa/,$ffasta;
    my $fout = "$b[$#b]" . '.mafft.fa'; 

    if($nseq > 1) {
     # align only if more than one
     my @b=split /\.fa/,$ffasta;
     my $fout = "$b[$#b]" . '.mafft.fa';
     my $str=DoMafft($ffasta,$fout);
     print MYMAFFTSCR "$str\n";
    } else {
     print "Only one sequence - copy: $ffasta\n";
     # just copy input in output
     my $cmd='mv '.$ffasta.' '.$fout;
     system($cmd);
    };
 }; 
 close(MYMAFFTSCR);

 exit 0;
};

# reads mafft aligned files and writes a consensus seq
CONS: {
  my @mafftf=<$mafft_dir/*.mafft.fa>;

  my $nexcov=0;
  my $flist = '>' . $out_dir . '/' . $spec . '_cov.txt';
  open(MYLIST,$flist);
  foreach my $file (@mafftf) {
    $nexcov++;
    my @a = split /\//,$file;
    my @b = split /\.mafft\.fa/,$a[$#a];
    print MYLIST "$b[0]\n";
  };
  close(MYLIST);
  print "$nexcov found\n";

  foreach my $ffasta (@mafftf) {
    MakeCons($ffasta,$out_dir);
  };
}
