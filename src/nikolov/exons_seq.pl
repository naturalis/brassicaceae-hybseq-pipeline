#!/usr/bin/env perl

use warnings;
use strict;
use Bio::SeqIO;
use Bio::Seq;
use Bio::SimpleAlign;
use Bio::AlignIO;

# The pipeline works in steps (check after the subroutines section) and at each step writes a bash script for more efficient utilization of the computational resources. This script carries out all processing after the YASRA alignment to the sequences of the final exons for each species separately. It assumes YASRA is run for each of three references (AT, AA, SI) separately and the sam files are provided as command line arguments. The script also requires a directory, called exons, that contains the sequences of the reference exons, their location in the ref-at.fasta, ref-aa.fasta, ref-si.fasta files and 4 files (see below) that contain the full path + name of each exon for each reference.  

# cutoff values - check the paper for definition
my $pid_cutoff=75.0;
my $score_cutoff=20;
# options for blat
my $blatoptions = "-t=dnax -q=dnax -stepSize=5 -repMatch=2253 -minScore=0 -minIdentity=0";
# options for mafft
my $mafftoptions = "--maxiterate 1000 --genafpair";
# the reference identifiers
my @identifiers = ("AT","AA","SI");
# command-line arguments that are set up below
my $dir_at="";
my $dir_aa="";
my $dir_si="";
my $spec_suffix="";
my $dir_root = "";
my $step;
#directory with the exon sequences
my $exon_db_dir = '/biodata/dep_tsiantis/grp_tsiantis/nikolov/exons/';
#file that contains the full path + exon-file-name - adjust when you change $exon_db_dir
#similarly for exons_AT.txt, exons_AA.txt, exons_SI.txt - this is a bit inconvinient and I meant to improve it in the future
my $exon_db = $exon_db_dir . 'exons.txt';

# the program writes contigs_AT, contigs_AA, contigs_SI, blats, maffts, and final_fasta directories in CWD 
my $blat_dir = 'blats';
my $mafft_dir = 'maffts';
my $final_fasta = 'final_fasta';

if( !( -d $blat_dir ) ) {
 print "$blat_dir does not exist\n";
 mkdir $blat_dir;
};

if( !( -d $mafft_dir ) ) {
 print "$mafft_dir does not exist\n";
 mkdir $mafft_dir;
};

if( !( -d $final_fasta ) ) {
 print "$final_fasta does not exist\n";
 mkdir $final_fasta;
};

# extract command-line arguments 
if($#ARGV != 5) {
 print STDERR "Provide paths and file names of contigs aligned to AT, AA, and SI in this order, and species abbreviation.\n";
 exit 0;
} else {
 # First 3 arguments are the path and name of the YASRA sam files in the YASRA_related_files directory
 $dir_at = $ARGV[0];
 $dir_aa = $ARGV[1];
 $dir_si = $ARGV[2];
 # The path with the directory with the species folders
 $dir_root = $ARGV[3];
 # The species name which is the name of the folder
 $spec_suffix = $ARGV[4];
 # The step of the pipeline
 $step = $ARGV[5];
};

# If you provide the complete path to YASRA sam files, comment these 3 lines out.
$dir_at = $dir_root . '/' . $dir_at;
$dir_aa = $dir_root . '/' . $dir_aa;
$dir_si = $dir_root . '/' . $dir_si;

print "$dir_at\n";
print "$dir_aa\n";
print "$dir_si\n";

my %alignedfile = ($identifiers[0] => $dir_at, $identifiers[1] => $dir_aa, $identifiers[2] => $dir_si);

# hash with exon number and identifiers
my %dict_exon = ();

# hash with sequenced exons
my %seq_exon = ();

# temporary variable
my @pslString = ();

#### Subroutines used in the main part of the pipeline

# All file names include the absolute path to the file

# Do the consensus calling using the Samtools/VarScan pipeline from the Yasra Sam file
# The script yasra_sam_consensus.pl must be in the path
sub RunYasraSam2Contigs($$) {
   my ($fsam,$suffix) = @_;
   my $print_dir = join("",("contigs_",$suffix));
   my $fcontignames = join("",("contigs_",$suffix,".txt"));
 
   # create directory if it does not exist
   if( !( -d $print_dir ) ) {
     print STDERR "$print_dir does not exist\n";
     mkdir $print_dir;
   };

   # run yasra_sam_consensus.pl in printdir
   chdir $print_dir;
   my $cmd = 'yasra_sam_consensus.pl ' . $fsam . ' > yasra_sam_'.$suffix.'.out 2>&1';
   system($cmd);
   my @lcontigs = <*.fasta>;
   print "Number of contigs_$suffix: $#lcontigs\n";
   chdir "..";

   my $fout = join("",(">",$fcontignames));
   open(MYCONTIGS,$fout);
   foreach my $ex (@lcontigs) {
     chomp($ex);
     print MYCONTIGS "$ex\n";
   };
   close(MYCONTIGS);
};

sub MakeDbfile($) {
  my ($dbfile) = @_;
  my @a = split /\./,$dbfile;
  my $froot = $a[0];
  open(MYF,$dbfile);
   my @exons=<MYF>;
  close(MYF);

  foreach my $id (@identifiers) {
    my $fname = join("",($froot,"_",$id,'.txt'));
    open(FOUT,">$fname");
    foreach my $ex (@exons) {
      chomp($ex);
      my $exout = join("",($ex,"_",$id,'.fasta'));
      print FOUT "$exout\n";
    };
    close(FOUT);
  };
};

# BLAT a file with a database and print the result in directory
# Use: DoBlat(path-to-database,seqfile,contig name,suffix,printing directory)
sub DoBlat($$$$$) {
  my ($database,$file,$contig,$suffix,$print_dir) = @_;

  #blat the contig 
  my $fout=join "",($print_dir,"/",$contig,"_",$suffix,'.psl');
  my $comm=join " ",("blat",$blatoptions,$database,$file,$fout,' > blat.out 2>&1');

  # due to weird lqueue stuff just print the command in a shell script
  #system($comm);

  return $comm;
};

sub DoMafft($$) {
  my ($filein,$fileout) = @_;

  #mafft filein 
  my $comm=join " ",('mafft',$mafftoptions,$filein,'>',$fileout,' 2> mafft.out');
  
  # due to weird lqueue stuff just print the command in a shell script
  #system($comm);

  return $comm;
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
  my $outstr = sprintf("%s\t%d\t%d\t%s\t%d\t%d\t%d\t%.2f\n", $tName, $tStart, $tEnd, $qName, $qStart, $qEnd, $pslScore, $percentIdentity);

  # add the contigs to the pslString
  @pslString = ();
  push(@pslString,$tName);
  push(@pslString,$qName);
  push(@pslString,$percentIdentity);
  push(@pslString,$pslScore);
  return $outstr;
 } else {
  my $outstr = "  WARNING: No sequence match found :(\n";
  return $outstr;
 };
};

sub MakeCons($$$) {
  my ($fmafft,$suffix,$fadir)=@_;

  # read the alignment and create the alignment object
  my $str = Bio::AlignIO->new(-file => $fmafft , -format => 'fasta');
  my $aln = $str->next_aln();

  $aln->uppercase();
  my $al = $aln->length;
  my $nseq = 0;
  # consensus sequence 
  my %consensus = ();

  foreach my $seq ($aln->each_seq) {
    $nseq++;
    if($nseq > 3) {
      my $did = $seq->display_id();
      my $ifcontg = substr($did,0,6);
      if ($ifcontg eq "Contig") {
       for (my $i = 1; $i <= $al; $i++) {
         my $res = $seq->subseq($i,$i);
         if(exists $consensus{$i}{$res}) {
           $consensus{$i}{$res}++;   
         } else {
           $consensus{$i}{$res}=1;
         };
        };
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
    } elsif($size == 2) {
      if (exists $col{'-'}) {
       $consensus{$i}{'cons'}=$key[1];
      } elsif($col{$key[0]} == $col{$key[1]}) {
       $consensus{$i}{'cons'}='N';
      } else {
       if($col{$key[0]} > $col{$key[1]}) {
        $consensus{$i}{'cons'}=$key[0];
       } else {
        $consensus{$i}{'cons'}=$key[1];
       }
      }
    } else {
      $consensus{$i}{'cons'}='N';
    };
  };

  my $cons_seq="";
  for (my $i = 1; $i <= $al; $i++) {
    $cons_seq = $cons_seq . $consensus{$i}{'cons'};
  }
  # remove gaps
  $cons_seq =~ tr/-//d;

  # make a new file object
  my @a = split /\//,$fmafft;
  my @b = split /\./,$a[$#a];
  my $fn = join("",(">",$fadir,'/',$b[0],"_",$suffix,'.fasta'));
  my $out_file = Bio::SeqIO->new(-file => $fn, -format => 'Fasta', -alphabet => 'dna');

  # make a new sequence object
  my $did = $b[0];
  my $out_seq  = Bio::Seq->new(-seq => $cons_seq, -display_id => $did, -alphabet => 'dna');

  # write to file
  $out_file->write_seq($out_seq);
};

sub FragmentOvlp($$$$) {
 my ($s1,$e1,$s2,$e2)=@_;

 my $ovl=0;
 if ($s1 >= $s2 && $s1 <= $e2) {$ovl++;};
 if ($e1 >= $s2 && $e1 <= $e2) {$ovl++;};
 if ($s2 >= $s1 && $s2 <= $e1) {$ovl++;};
 if ($e2 >= $s1 && $e2 <= $e1) {$ovl++;};
 if ($ovl > 0) {$ovl=1;};
 return $ovl;
};

#### End of subroutines section

# STEP 1: extract the contigs consensus sequence from the YASRA sam file (uses in addition yasra_sam_consensus.pl and makeconsensus.sh scripts) and prepare a blat script for homology assessment.
# STEP 2: analize the blat results and prepare a mafft script and files for contig alignment to the references
# STEP 3: analyze the mafft results and infer the final reference independent consensus sequence

# if BLAT is done use this switch to do PSL analysis and write the mafft script
 goto DICT_EX if ($step == 2);
# if statistics and consensus calling is to be done use this switch
 goto CONS if ($step == 3);

# Extract contigs in files and separate folders
foreach my $id (@identifiers) {
   RunYasraSam2Contigs($alignedfile{$id},$id);

   # sort and trim contig names after RunYasraSam2Contigs consensus calling
   my $fname='contigs_' . $id . '.txt';
   open(MYCONTIG,$fname);
    my @fcontig = <MYCONTIG>;
   close(MYCONTIG);

   my @sort_fcontig = sort { ($a =~ /(\d+)/)[0] <=> ($b =~ /(\d+)/)[0] } @fcontig;

   open(MYCONTIG,">$fname");
   foreach my $cont (@sort_fcontig) {
     chomp($cont);
     my @nname=split(/\./,$cont);
     print MYCONTIG "$nname[0]\n";
   };
   close(MYCONTIG);
};

BLAT: {
  # Make the database files for BLAT
  MakeDbfile($exon_db);

  open(MYBLATSCR,">doblat.sh");
  print MYBLATSCR "#!/bin/bash\n";

  # Start the BLAT alignment section
  foreach my $id (@identifiers) {
    my @a = split /\./,$exon_db;
    my $db = "$a[0]" . "_" . "$id" . '.txt'; 

    my $fcontgs="contigs_" . "$id" . '.txt';
    open FCONTIGS, "$fcontgs";
     my @fcont=<FCONTIGS>;
    close FCONTIGS;

    foreach my $contg (@fcont) {
      chomp($contg);
      my $fc = "contigs_" . "$id" . '/' . "$contg" . '.fasta';
      my $str = DoBlat($db,$fc,$contg,$id,$blat_dir);
      print MYBLATSCR "$str\n";
    };
  };

  close(MYBLATSCR);
};

# after the blat block exit
exit 0;

DICT_EX: {
  # make dictionary
  open(FEX,$exon_db);
   my @fexons=<FEX>;
  close(FEX);

  # loop over exon homologues
  foreach my $fex (@fexons) {
    chomp($fex);
    # open fasta file
    my $fin = "$fex" . '.fasta';
    my $in = Bio::SeqIO->new(-file => $fin, -format => 'Fasta', -alphabet => 'dna');
    # make exon ID
    my @a = split /\//,$fex;
    my $exid = $a[$#a];

    # go over the sequences in the exon file
    while(my $rec = $in->next_seq) {
      my $did = $rec->display_id();
      # add to dictionary    
      $dict_exon{$did}=$exid;    
    }; 
  }; 
};

my %contig_exon=();

ASSIGNMENT_BY_BOUNDARY: {
 
 # loop over references
 foreach my $ref (@identifiers) {

  # input file names
  my $db_fname=$exon_db_dir.'/'.$ref.'_exon_enum.txt';
  my $contigs_fname='contigs_'.$ref.'.txt';
  
  # read exons
  my %db_st=();
  my %db_en=();
  open(MYDB,$db_fname);
  my @db_lines=<MYDB>;
  close(MYDB);
  foreach my $l (@db_lines) {
   chomp($l);
   my @a=split(/\s/,$l);
   $db_st{$a[0]}=$a[1];
   $db_en{$a[0]}=$a[2];
  };

  # read contigs
  my %ct_st=();
  my %ct_en=();
  open(MYCT,$contigs_fname);
  my @ct_lines=<MYCT>;
  close(MYCT);
  foreach my $l (@ct_lines) {
   chomp($l);
   my @a=split(/_/,$l);
   $ct_st{$l}=$a[$#a-1];
   $ct_en{$l}=$a[$#a];
  }

  # search all for a match
  foreach my $cont (keys %ct_st) {
   my $start = $ct_st{$cont};
   my $end   = $ct_en{$cont};
   foreach my $dbex (keys %db_st) {
    my $ovlp = FragmentOvlp($start,$end,$db_st{$dbex},$db_en{$dbex});
    if ($ovlp == 1) {
     $contig_exon{$cont}{$dbex}=1;
    }
   }
  };

 };

 # print out for debug
 open(MYCTEX,">contigs_exon_list.txt");
 foreach my $cont (sort keys %contig_exon) {
  print MYCTEX "$cont ";
  foreach my $name (keys %{ $contig_exon{$cont} }) {
   print MYCTEX "$name ";
  }
  print MYCTEX "\n";
 }
 close(MYCTEX);

};

PSL: {
  my @psl=<$blat_dir/*>;

  open(MYPSLOUT,">psl_stats.out");
  open(MYPSLERR,">psl_err.out");
  open(MYEXCONTOUT,">exon_contig_pair.txt");

  foreach my $fpsl (@psl) {

    # read psl and find the highest score exon
    my $str = ReadPsl($fpsl);

    # check psl output and print stats
    if(index($str,'WARNING') > 0) {
      chomp($fpsl);
      print MYPSLERR "$fpsl\n";
    } else {

      # analyze data and pileup for mafft
      my $exon = $pslString[0];
      my $contig = $pslString[1];
      my $pid = $pslString[2];
      my $score = $pslString[3];

      # PID cutoff
      if($pid >= $pid_cutoff && $score > $score_cutoff) {
       # YASRA comparison cutoff
       if(exists $contig_exon{$contig}{$exon}) {
        print MYPSLOUT "$str";

        my $exonid = $dict_exon{$exon};

        if(exists $seq_exon{$exonid}) {
         $seq_exon{$exonid}=join " ",($seq_exon{$exonid} , $contig);
        } else {
         $seq_exon{$exonid}=$contig; 
         # copy exon fasta file
         my $fex1 = "$exon_db_dir" . "$exonid" . '.fasta'; 
         my $cmd1 = "cp " . "$fex1" . " " . "$mafft_dir" . '/' . "$exonid" . '.fasta';
         system($cmd1); 
        };

        # append the contig seq to the three refs
        my $fex = "$mafft_dir" . '/' . "$exonid" . '.fasta';
        my @a = split /\_/, $contig;
        my $id = $a[1];
        my $fcontig = "contigs_" . "$id" . '/' . "$contig" . '.fasta';
        my $cmd = 'cat ' . "$fcontig" . ' >> ' . "$fex";
        system($cmd);
   
        # print the exon contig pair
        print MYEXCONTOUT "$exonid $contig\n";   

       } else {
        print MYPSLERR "$contig $exon pair NOT FOUND in contig_exon list\n";
       };
      } else {
       chomp($fpsl);
       print MYPSLERR "$fpsl - PID or SCORE less than CUTOFF\n";
      };
    };
  };

  # print final statistics
  my $nseq=0;
  for (sort keys %seq_exon) {$nseq++};
  print MYPSLOUT "Number of sequenced exons: $nseq\n";

  close(MYEXCONTOUT);
  close(MYPSLOUT);
  close(MYPSLERR);
};

MAFFT: {
  my @mafftf=<$mafft_dir/*>;

  open(MYMAFFTSCR,">domafft.sh");
  print MYMAFFTSCR "#!/bin/bash\n";

  foreach my $ffasta (@mafftf) {
    my @b=split /\./,$ffasta;
    my $fout = "$b[$#b-1]" . '.mafft.fa';  
    my $str = DoMafft($ffasta,$fout);
    print MYMAFFTSCR "$str\n";
  };

  close(MYMAFFTSCR);
};

# exit after writing the script
exit 0;

# reads mafft aligned files and writes a consensus seq
CONS: {
  my @mafftf=<$mafft_dir/*.mafft.fa>;

  foreach my $ffasta (@mafftf) {
    MakeCons($ffasta,$spec_suffix,$final_fasta);
  };
}

