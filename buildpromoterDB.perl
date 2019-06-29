use strict;
use List::Util qw(min max);

if( @ARGV < 4) {
   print "makePromoterDB <ChromosomeFastaFile> <GeneGffFile> <Length> <UTR_Length>\n";
   exit();
}

open (DB, "gunzip -c $ARGV[0] |") || die "can`t open file $ARGV[0]";
open (my $gfffile, "$ARGV[1]") || die "can't open file $ARGV[1]";

my $seq;
my $prom_len = $ARGV[2];
my $UTR_len = $ARGV[3];


# read genome
my @chromosomes;
my @chromosomes_names;
my $cur_chr_num = 0;
my $cur_chr = <DB>;
my $cur_chrom_seq ="";
sub clean_chr {
	my $chr = shift(@_);
	chomp($chr);
	$chr = substr($chr,1);
	if(index($chr, " ")!= -1 || index($chr, "\t") != -1) {
		$chr = (split(' ', $chr))[0];
	}
	return $chr;
}

my $cur_line;
while($cur_line = <DB>) {
	if(substr($cur_line,0,1) eq ">") {
		$chromosomes[$cur_chr_num] = $cur_chrom_seq;
		$cur_chr = clean_chr($cur_chr);
		$chromosomes_names[$cur_chr_num] = $cur_chr;
		$cur_chr_num++;
		$cur_chr = $cur_line;
		$cur_chrom_seq = "";
	} else {
		chomp($cur_line);
		$cur_chrom_seq = $cur_chrom_seq . $cur_line;
	}
}
# dump last chromosome
$chromosomes[$cur_chr_num] = $cur_chrom_seq;
$cur_chr = clean_chr($cur_chr);
$chromosomes_names[$cur_chr_num] = $cur_chr;

my %index;
@index{@chromosomes_names} = (0..$#chromosomes_names);

my $mylastgeneend_pos=0;
my $mylastchromosome_num=0;

while(my $gffline = <$gfffile>)
{
  if(substr($gffline, 1,1) ne "#") {
  	chomp($gffline);
  	my @genefields = split(/\t/, $gffline);
  	my %geneinfo = split /[;=]/, $genefields[8];
  	my $genename = "";
				
	if(exists $geneinfo{"Name"}) {
		$genename = $geneinfo{"Name"};
		# get rid of the odd inconsistent appendix for gene/protein name
		
		$genename =~ s/\.g$//g;
		$genename =~ s/\.[0-9]$//g;

	} elsif(exists $geneinfo{"Parent"}) {
		$genename = $geneinfo{"Parent"};
		# get rid of weird embalishments
		$genename =~ s/mRNA\.|CDS//g;

	} elsif(exists $geneinfo{"ID"}) {
		$genename = $geneinfo{"ID"};
	}
	
				
	if($genefields[2] eq "gene") {
  		my $chromosome = $genefields[0];
		if( !(exists $index{$chromosome})) {
			print (STDERR "ERROR! " . $chromosome . " does not exists! aborting.\n");
			die;
		}
		my $chromosome_num = $index{$chromosome};
  		my $dir = $genefields[6];
  		
  		if($chromosome_num != $mylastchromosome_num) {
  			$mylastchromosome_num = $chromosome_num;
  			$mylastgeneend_pos=0;
  		}

  		print ">". $genename . "\n";
  		my $prom;


		if($dir eq "+") {
 	 	     	my $loc = $genefields[3];
 	 	     	my $pstart = $loc - $prom_len;
 	 	     	if ($pstart <0 ) { $pstart=0; }

 	 	     	# get sequence until the previous gene
 	 	     	if($pstart < $mylastgeneend_pos) {
 	 	     		# one exterme case - overlapping genes. In this case, promoter is NULL. Maybe there is a better way of dealing with this. The problem
 	 	     		# of usign gene coding sequence as promoter is that coding sequences have much higher conservation rates that promoters which will confiuse downstream analysis
 	 	     		if($mylastgeneend_pos >= $loc) {
 	 	     			$prom = "";
 	 	     		} else {
 						$prom = substr($chromosomes[$chromosome_num], $mylastgeneend_pos, $prom_len - $mylastgeneend_pos + $pstart + $UTR_len);
 					}
 				} else {
	     	    	$prom = substr($chromosomes[$chromosome_num], $pstart, $prom_len + $UTR_len);
	     	    }
  		} else {
     	     	my $loc = $genefields[4]; 
  				## peek at the next line to check when the next gene start. That's not nice programming.
  				my $curfilepos = tell $gfffile;
  				my $peekgene = <$gfffile>;
  				my @peekgenefields = split(/\t/, $peekgene);
  				my $nextgenestart = $peekgenefields[3];
  				seek($gfffile, $curfilepos,0);
  				
  				if($loc + $prom_len > $nextgenestart) {
  					# check for overlapping genes...
  					if($nextgenestart < $loc) {
  						$prom="";
  					} else {
   	     				$prom = substr($chromosomes[$chromosome_num], $loc - $UTR_len, $nextgenestart - $loc + $UTR_len);
   	     			}
  				} else {
   	     			$prom = substr($chromosomes[$chromosome_num], $loc - $UTR_len, $prom_len + $UTR_len);
   	     		}
      	     	# reverse complement sequence
             	$prom = reverse $prom;
     	     	$prom =~ tr/ACGT/TGCA/;
		}
   		print ('N' x ($prom_len - length($prom)) . $prom . "\n");
   		$mylastgeneend_pos = max($genefields[3], $genefields[4]);

   		}
    }
}
