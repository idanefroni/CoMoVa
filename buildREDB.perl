#use POSIX;
use strict;
$|=1;

if( @ARGV < 2) {
   print "buildREDB <PromoterFasta> <MotifsFasta>\n";
   exit();
}
open UP1KB, $ARGV[0];
open PLACE, $ARGV[1];

my @motifs;
my @motifname;
my @motifshort;

my $motname;
my $motifcount;
my $revmot;
my $seq;
my $mot;
my $name;
my $start;
my $explorewid=0;

my $prm_start= 0;

my $prm_end=3000;
my $flank_length=20;

my $prm_len= $prm_end - $prm_start;
my $i=0;

#read motifs
while(my $line = <PLACE>) { 
   $seq = <PLACE>;

   chomp($line);
   if(substr($line,0,1) eq '>') {

   	$line = substr($line,1);

   	$seq =~s/\n//;
   	push(@motifshort, $line);
   	push(@motifname, $seq);

   	$seq =~s/N/[A|C|G|T]/g;
   	$seq =~s/X/[A|C|G|T]/g;
   	$seq =~s/W/[A|T]/g;
   	$seq =~s/R/[A|G]/g;
   	$seq =~s/Y/[C|T]/g;
   	$seq =~s/M/[A|C]/g;
   	$seq =~s/S/[G|C]/g;
   	$seq =~s/s/[G|C]/g;
   	$seq =~s/K/[G|T]/g;
   	$seq =~s/H/[A|C|T]/g;
   	$seq =~s/V/[A|C|G]/g;
   	$seq =~s/h/[A|C|T]/g;
   	$seq =~s/D/[A|G|T]/g;
   	$seq =~s/B/[C|G|T]/g;
 
   	push(@motifs, $seq);
   	$i++;
   }
}

close(PLACE);
#read seqeunce

print "Locus,Motif,Variant,Context,Pos,Rev,Region\n";

seek UP1KB,0,0;
$name = <UP1KB>;
chomp($name);

while(my $line = <UP1KB>) {
   $seq = "";
   chomp($line);
   while((substr $line,0,1) ne ">" && $line ne "") {
	$seq = $seq . $line;
        $line = <UP1KB>;
	chomp($line);
   } 
   if(index($name,";")>0 ) {
   	$name = substr(substr ($name, 0, index($name,";")),0);
   }

   $name = substr($name,1);
 
   my $motindex=0;

   foreach(@motifs) {
     $motifcount=0;

     $mot = $_;
     $revmot = reverse $mot;
     $revmot =~ tr/ACGT/TGCA/;
     $revmot =~ tr/[]/][/;
     $revmot =~ tr/{}/}{/;

     $motname = $motifname[$motindex];
     my @xpos;

     while( $motname =~ /[X|R|Y|B|K|D|M|V|S|H|W]/g) { push @xpos, scalar (pos($motname)); }
	
     if($seq=~ m/($mot)/) {

        while($seq =~ m/($mot)/g) {
	      #collapse pairs variants
	      my $flank5;
	      my $flank3;
	      if(pos($seq)> $flank_length+length($1)) {
	          $flank5 = substr($seq, (pos($seq))- $flank_length - length($1), $flank_length);
	      } else {
		  $flank5 = substr($seq, 0, pos($seq)-length($1));
	          $flank5 = "-" x ($flank_length-length($flank5)) . $flank5;
              }

	      if(pos($seq) < length($seq)- $flank_length) {
	        $flank3 = substr($seq, (pos($seq)), $flank_length);
	      } else {
	        $flank3 = substr($seq, pos($seq));
	      	$flank3 = $flank3 . "-" x ($flank_length-length($flank3));
	      }

      	      # get X's
	      my @varX;
	      foreach(@xpos) { 
	      		push @varX, substr($seq, pos($seq) - length($1) + scalar $_ -1, 1);
	      }

	      print $name . "," . $motifshort[$motindex] . "," . join("",@varX) . "," . $1 . "," . ($prm_end - (pos($seq)- length($1))) . ",F," . $flank5 . $1 . $flank3 . "\n";
        }
      }

     if($seq=~ m/($revmot)/) {
        while($seq =~ m/($revmot)/g) {
	      my $tmot = $motname;
	      $tmot = reverse $tmot;
	      $tmot =~ tr/ACGT/TGCA/;
	      my $flank5;
	      my $flank3;
	      if(pos($seq)> $flank_length+length($1)) {
	          $flank5 = substr($seq, (pos($seq))- $flank_length - length($1), $flank_length);
	      } else {
		  $flank5 = substr($seq, 0, pos($seq)-length($1));
	          $flank5 = "-" x ($flank_length-length($flank5)) . $flank5;
              }

	      if(pos($seq) < length($seq)- $flank_length) {
		 $flank3 = substr($seq, (pos($seq)), $flank_length);
	      } else {
		$flank3 = substr($seq, pos($seq));
	      	$flank3 = $flank3 . "-" x ($flank_length-length($flank3));
	      }

      	      # get X's
	      my @varX;
	      foreach(@xpos) { 
		push @varX, substr($seq, pos($seq) - scalar $_, 1);
	      }
	      my $vars = join("",@varX);
	      $vars =~ tr/ACGT/TGCA/;

	      print $name . "," . $motifshort[$motindex] . "," . $vars . "," . $1 . "," . ($prm_end - pos($seq)- length($1)) .  ",R," .  $flank5 . $1 . $flank3 . "\n";
        }
      }
      $motindex++;
   }
   $name = $line;
}

