#!/usr/bin/perl

#use Git::Repository;
## start from an existing repository
#$r = Git::Repository->new( git_dir => $ENV{'GITDIR'} );


#use String::Approx 'amatch';
#use String::Approx 'adistr';
#use String::Approx 'aindex';

use File::Basename; # Parse file paths into directory, filename and suffix CA


my %codons_all;
my %codons_all = (
'GCT' => 'A', 'GCC' => 'A', 'GCA' => 'A', 'GCG' => 'A', 'TTA' => 'L',
'TTG' => 'L', 'CTT' => 'L', 'CTC' => 'L', 'CTA' => 'L', 'CTG' => 'L',
'CGT' => 'R', 'CGC' => 'R', 'CGA' => 'R', 'CGG' => 'R', 'AGA' => 'R',
'AGG' => 'R', 'AAA' => 'K', 'AAG' => 'K', 'AAT' => 'N', 'AAC' => 'N',
'ATG' => 'M', 'GAT' => 'D', 'GAC' => 'D', 'TTT' => 'F', 'TTC' => 'F',
'TGT' => 'C', 'TGC' => 'C', 'CCT' => 'P', 'CCC' => 'P', 'CCA' => 'P',
'CCG' => 'P', 'CAA' => 'Q', 'CAG' => 'Q', 'TCT' => 'S', 'TCC' => 'S',
'TCA' => 'S', 'TCG' => 'S', 'AGT' => 'S', 'AGC' => 'S', 'GAA' => 'E',
'GAG' => 'E', 'ACT' => 'T', 'ACC' => 'T', 'ACA' => 'T', 'ACG' => 'T',
'GGT' => 'G', 'GGC' => 'G', 'GGA' => 'G', 'GGG' => 'G', 'TGG' => 'W',
'CAT' => 'H', 'CAC' => 'H', 'TAT' => 'Y', 'TAC' => 'Y', 'ATT' => 'I',
'ATC' => 'I', 'ATA' => 'I', 'GTT' => 'V', 'GTC' => 'V', 'GTA' => 'V',
'GTG' => 'V',
'TAA' => 'STOP',
'TGA' => 'STOP',
'TAG' => 'STOP',);

my %codons;
$codons{"AAT"} = 1;
$codons{"AAG"} = 1;
$codons{"ATT"} = 1;
$codons{"ATG"} = 1;
$codons{"ACT"} = 1;
$codons{"ACG"} = 1;
$codons{"AGT"} = 1;
$codons{"AGG"} = 1;
$codons{"TAT"} = 1;
$codons{"TAG"} = 1;
$codons{"TTT"} = 1;
$codons{"TTG"} = 1;
$codons{"TCT"} = 1;
$codons{"TCG"} = 1;
$codons{"TGT"} = 1;
$codons{"TGG"} = 1;
$codons{"CAT"} = 1;
$codons{"CAG"} = 1;
$codons{"CTT"} = 1;
$codons{"CTG"} = 1;
$codons{"CCT"} = 1;
$codons{"CCG"} = 1;
$codons{"CGT"} = 1;
$codons{"CGG"} = 1;
$codons{"GAT"} = 1;
$codons{"GAG"} = 1;
$codons{"GTT"} = 1;
$codons{"GTG"} = 1;
$codons{"GCT"} = 1;
$codons{"GCG"} = 1;
$codons{"GGT"} = 1;
$codons{"GGG"} = 1;

# Read in barcodes
#   FORMAT fasta
# Store fasta in hash
open(TMP, "$ARGV[0]") or die "Couldn't open $ARGV[0]: $?\n"; # open library fasta CA
my %barcode_lines;
my $header = "";
my $seq = "";
while ($line = <TMP>){ # read fasta line by line CA
    $line =~ s/\n//; #remove \n CA
    if ($line =~ /\>(\S+)/){ # if header CA

	if ($header){
	    $barcode_lines{$header} = $seq;
#	    printf("%-50s %s %6d %4d\n",$header, $seq);
	}

	$seq = "";
	$header = $1; # store header in header CA
    } else {
	$seq .= $line; # else, store seq in seq CA
    }
}
close(TMP);
$barcode_lines{$header} = $seq;

my $barcode_length = 10;
$barcode_length = $ARGV[3] if (@ARGV >= 4);

my %barcodes;
foreach $key (sort keys %barcode_lines){ # iterate over sorted fasta headers CA

    my @toks = ($key,$barcode_lines{$key});  # creates array with header and seq CA

TRYAGAIN: # label CA
    if ($toks[1] =~ /([GCAT]{$barcode_length})NNK([GCAT]{$barcode_length})/){ # checks if seq is 10 characters of GCAT then NNK then 10 characters of GCAT CA
	$barcodes{$toks[0]}{'bar1'} = $1; # store 10 char before NNK CA
	$barcodes{$toks[0]}{'bar2'} = $2; # store 10 char after NNK CA

        # reverse the DNA sequence
        my $revcomp1 = reverse($barcodes{$toks[0]}{'bar1'});

	# complement the reversed DNA sequence
        $revcomp1 =~ tr/ACGTacgt/TGCAtgca/;

        # reverse the DNA sequence
        my $revcomp2 = reverse($barcodes{$toks[0]}{'bar2'});

	# complement the reversed DNA sequence
        $revcomp2 =~ tr/ACGTacgt/TGCAtgca/;

	# stores reverse comp CA
	$barcodes{$toks[0]}{'revbar1'} = $revcomp1; 
	$barcodes{$toks[0]}{'revbar2'} = $revcomp2;

	# Check if we find the barcode in base sequence elsewhere..
	# Not sure that this block of code runs? CA
	my $b1 = $barcodes{$toks[0]}{'bar1'};
	my $b2 = $barcodes{$toks[0]}{'bar2'};
	my $r1 = $barcodes{$toks[0]}{'revbar1'};
	my $r2 = $barcodes{$toks[0]}{'revbar2'};
	if ($ref_seq =~ /$b1/ || $ref_seq =~ /$b2/ || $ref_seq =~ /$r1/ || $ref_seq =~ /$r2/){
	  print "Barcode: ".$toks[0]." requires longer barcode length than: $barcode_length\n";
	  $barcode_length += 3;
	  goto TRYAGAIN;
	}

	print "Barcode: ".$toks[0]." ".$barcodes{$toks[0]}{'bar1'}." NNK ".$barcodes{$toks[0]}{'bar2'} ." RC: ".$barcodes{$toks[0]}{'revbar2'}." NKK ".$barcodes{$toks[0]}{'revbar1'} ." wt_codon: ".$wt_codon." wt_aa:".$codons_all{$wt_codon}."\n";
    } else {

	$ref_name = $toks[0];
	$ref_seq  = $toks[1];

	$ref_seq_rc = reverse($ref_seq);
	$ref_seq_rc =~ tr/ACGTacgt/TGCAtgca/;

	print "No NNK this is the reference: $toks[0]\n\t$ref_seq\n\n$ref_seq_rc\n\n";
    }

}
close(TMP);


# PROCESS FASTA FILE
my $header = "";
my $seq    = "";
my $num_total_empty = 0;
my $total_num_reads_with_partial_barcodes_size = 0;
my $fname =fileparse($ARGV[1]);

if ($ARGV[1] =~ /fasta$/){ # changes file extension to data, if necessary CA
    $fname =~ s/.fasta/.data/;
} else {
    $fname = $fname.".data";
}
print "Writing to $fname\n"; # writing to example reads data CA
open(OUT, ">$fname");

print OUT "FNAME READ SEQLEN NUMBAR NUMNT CODON POS AA\n";

open(TMP, "$ARGV[1]") or die "Couldn't open $ARGV[1]: $?\n"; #open example reads fasta CA

my $allow_partial=0;
$allow_partial = $ARGV[2] if (@ARGV >= 3);

while ($line = <TMP>){ # read fasta line by line CA
    $line =~ s/\n//; # remove \n CA

    if ($line =~ />(\S+)/){ # if header CA
	my $h = $1; # stores header in h CA
	
	if ($header ne ""){ # if header not empty (header and seq are assigned) CA

	    # Sequence length
	    my $seq_length = length($seq);

	    # Count number of barcodes found in this sequence
	    my %partial_barcodes;
	    my %which_barcodes;
	    my $num_barcodes = 0;
	    my @barcodes;
	    my @barcode_positions;
	    foreach $barcode (keys %barcodes){ # iterates over keys in barcodes{header}{b1/b2/r1/r2} CA

		# Look for matched barcode bar1 and bar2
		my $b1 = $barcodes{$barcode}{'bar1'};
		my $b2 = $barcodes{$barcode}{'bar2'};
		my $found_barcode = 0;
		if ($seq =~ /$b1([AGCT]{3})$b2/){ # if b1 codon b2 CA
		    $num_barcodes++;
		    push @barcodes, $1; # push codon to barcodes array CA
		    push @barcode_positions, $barcode; # push header to barcode_positions array CA
		    $found_barcode = 1;
		}

		# Look for matched barcode revbar1 and revbar2
		# Repeat for r1 and r2 CA
		my $r1 = $barcodes{$barcode}{'revbar1'};
		my $r2 = $barcodes{$barcode}{'revbar2'};
		if ($seq =~ /$r2([AGCT]{3})$r1/){
		    $num_barcodes++;
		    my $nnk = $1;
		    my $revcomp = reverse($nnk);
		    $revcomp  =~ tr/ACGTacgt/TGCAtgca/;

		    push @barcodes, $revcomp;
		    push @barcode_positions, $barcode;
		    $found_barcode = 1;
		}

		
		#  ERROR CODES FOR NOT FINDING a matched barcode
		# accounts for partial barcodes (contains upstream but not downstream, or vis versa) CA
		if (!$found_barcode){
		  if ($seq =~ /$b1([AGCT]{3})/ && $seq !~ /$b2/){
		    $partial_barcodes{1}++;
		    $which_barcodes{$barcode."-1"}++;
		    
		    if ($allow_partial){
			$num_barcodes++;
			push @barcodes, $1; # push codon to barcodes array CA
			push @barcode_positions, $barcode; # push header to barcode_positions array CA
			$found_barcode = 1;
		    }


		  }

		  if ($seq =~ /([AGCT]{3})$b2/ && $seq !~ /$b1/){
		    $partial_barcodes{2}++;
		    $which_barcodes{$barcode."-2"}++;

		    if ($allow_partial){
			$num_barcodes++;
			push @barcodes, $1;
			push @barcode_positions, $barcode;
			$found_barcode = 1;
		    }

		  }

		  if ($seq =~ /([AGCT]{3})$r1/ && $seq !~ /$r2/){
		    $partial_barcodes{3}++;
		    $which_barcodes{$barcode."-3"}++;

		    if ($allow_partial){
			$num_barcodes++;
			my $nnk = $1;
			my $revcomp = reverse($nnk);
			$revcomp  =~ tr/ACGTacgt/TGCAtgca/;

			push @barcodes, $revcomp;
			push @barcode_positions, $barcode;
			$found_barcode = 1;
		    }
		  }

		  if ($seq =~ /$r2([AGCT]{3})/ && $seq !~ /$r1/){
		    $partial_barcodes{4}++;
		    $which_barcodes{$barcode."-4"}++;

		    if ($allow_partial){
			$num_barcodes++;
			my $nnk = $1;
			my $revcomp = reverse($nnk);
			$revcomp  =~ tr/ACGTacgt/TGCAtgca/;

			push @barcodes, $revcomp;
			push @barcode_positions, $barcode;
			$found_barcode = 1;
		    }

		  }

		}

		
	      }

	    # Count number of non-barcodes(number of mutations from ref - 5bp*(number of barcodes))
	    my $num_nonbarcodes = -1;
	    for ($i = 0; $i <= 60; $i+=10){ # loops that selects smaller and smaller sub seqs CA
		my $start_seq = substr($seq,$i,10); # start of seq CA
		my $end_seq   = substr($seq,length($seq)-($i+11),10); # end of seq CA
		my $sub_seq   = "";
		if ($seq =~ /$start_seq(\S+)$end_seq/){
		    $sub_seq = $start_seq.$1.$end_seq; # if start of seq and end of seq do not overlap, sub seq = start seq - end seq CA
		} 


		$num_nonbarcodes = -1;
#		print "Ref seq: \n$ref_seq\nsearch string:\n$start_seq\n$end_seq\n";
		if ($ref_seq =~ /$start_seq(\S+)$end_seq/){ # if sub seq start and end found in ref seq CA

		    my $ref_match_substr = $start_seq.$1.$end_seq;
#		    print "YES! ".length($ref_match_substr)." - ".length($sub_seq)." START: $start_seq END: $end_seq\n";
		    if (length($ref_match_substr) != length($sub_seq)){ # if length of sub seq != length of ref seq match CA
			$num_nonbarcodes = -3 + -abs(length($ref_match_substr)-length($sub_seq)); # num_nonbarcodes is determined by the difference in seq length CA
			# what if there is also mismatches? CA
			last; # exit for loop
		    }
		    $num_nonbarcodes = 0; # num_nonbarcodes = 0 if lengths are the same CA
		    my $mask = $ref_match_substr ^ $sub_seq; # XOR operation that stores differences of seqs in mask CA
		    my $num_mismatches = 0;
		    while ($mask =~ /[^\0]/g) {
			$num_mismatches++;
		    }
		
		    $num_nonbarcodes = $num_mismatches; # num_nonbarcodes is determined by mismatches CA
		}

		last if ($num_nonbarcodes != -1); # exit for loop if num_nonbarcodes is found CA

		$num_nonbarcodes = -2;
		# repeat for rc CA
		if ($ref_seq_rc =~ /$start_seq(\S+)$end_seq/){

		    my $ref_match_substr = $start_seq.$1.$end_seq;
		    if (length($ref_match_substr) != length($sub_seq)){
			$num_nonbarcodes = -3 + -abs(length($ref_match_substr)-length($sub_seq));
			last;
		    }
		    $num_nonbarcodes = 0;
		    my $mask = $ref_match_substr ^ $sub_seq;
		    my $num_mismatches = 0;
		    while ($mask =~ /[^\0]/g) {
			$num_mismatches++;
		    }
		
		    $num_nonbarcodes = $num_mismatches;

		}

		(last) if ($num_nonbarcodes != -2); # exit for loop if num_nonbarcodes is found CA

#		print "NO! START: $start_seq END: $end_seq\n";
	    }


#	    print "Num nonbc: $num_nonbarcodes\n";


	    
            ############ DEAL WITH ERRORS ####################
	    #  EMPTY - no matched bar codes found
	    if ($num_barcodes == 0){ # if b1 codon b2 or r1 codon r2 not found CA
		$num_total_empty++;

		# Find out which barcode with largest number of instances
		my @max_partial_barcode_found = (sort {$which_barcodes{$a} <=> $which_barcodes{$b}} keys %which_barcodes);

		my $partial_barcodes_size = $partial_barcodes{1} +$partial_barcodes{2} +$partial_barcodes{3} +$partial_barcodes{4};
		print "$header [$seq_length] has ".$partial_barcodes_size." partial barcodes, [ ".$partial_barcodes{1}."]  [ ".$partial_barcodes{2}."]  [ ".$partial_barcodes{3}."]  [ ".$partial_barcodes{4}."] ".join(", ",@max_partial_barcode_found)." max(".$which_barcodes{$max_partial_barcode_found[0]}.")]\n";
#		print NOBAR ">$header\n$seq\n";

	    } else { 

		################### MATCHED BARCODES ###################################

#		print ">$header\n$seq\n";
#		printf "%30s %3d %2d %3s %3d\n", $header, $num_barcodes, $seq_length,$barcodes[0],$barcode_positions[0];
		printf OUT "%30s %30s %3d %4d %4d ", $fname, $header, $seq_length, $num_barcodes,$num_nonbarcodes; # print reads with barcodes to example reads data CA

		my $barcode_str     = "";
		my $barcode_pos_str = "";
		my $barcode_aa_str  = "";
		for ($i = 0; $i <= $#barcodes;$i++){
		    if ($barcode_str ne ""){
			$barcode_str     .= ",";
			$barcode_pos_str .= ",";
			$barcode_aa_str  .= ",";
		    }
		    $barcode_str     .= $barcodes[$i];
		    $barcode_pos_str .= $barcode_positions[$i];
		    $barcode_aa_str  .= $codons_all{$barcodes[$i]};

		}

		printf OUT "%15s %15s %15s\n", $barcode_str,$barcode_pos_str, $barcode_aa_str; # print codon, position, amino acid CA
		

	    }

	    
	}

	# stores h in header if header is empty CA
	$header = $h;
	$seq    = "";
	next;
    } 

    $seq .= $line; # stores seq in seq CA
}
close(TMP);
close(OUT);
#close(NOBAR);
print "Number of empty barcode: $num_total_empty\n"; # number of sequences without barcodes CA

# barcode = 
# nonbarcode = 