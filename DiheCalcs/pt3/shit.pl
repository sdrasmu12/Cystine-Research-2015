use HackaMol;
use Modern::Perl;
use YAML::XS qw(LoadFile);
use POSIX;

# setup

my $t1 = time;

my $file = shift;
#my $mol = HackaMol->new->read_file_mol($file);
say 'Done Reading';

my $root = $file;
$root =~ s/\.xyz//;
my $histo = LoadFile("$root\_dist_chi3.yaml");
my %histo = %{$histo};

my $tot= 0;
foreach my $dist (sort { $a <=> $b } keys %histo){
  foreach my $chi3 (sort{ $a <=> $b } keys %{$histo{$dist}}){
    $tot++;
    say "total configs: $tot ", "number of configs for $dist $chi3:", scalar(@{$histo{$dist}{$chi3}}); 
  }
}

=cut;

foreach my $t (0 .. $mol->tmax){
  $mol->gt($t);
  my $nnlength = ceil(($nnbond->bond_length)/$partit2);
  my $chi3 = ceil($dihe3->dihe_deg/$partit);
  push @{$fresults{$nnlength}{$chi3}}, $t;
}



open(my $fchi3, ">", "$root\_dist_chi3.txt") or die "couldn't open";
open(my $fchi12, ">", "$root\_chi1_chi2.txt") or die "couldn't open";


my $tot= 0;
foreach my $dist (sort { $a <=> $b } keys %fresults){
  foreach my $chi3 (sort{ $a <=> $b } keys %{$fresults{$dist}}){
    $tot++;
    say "total configs: $tot ", "number of configs for $dist $chi3:", $fresults{$dist}{$chi3}; 
  }
}

say "total number of configs ", $tot;

#use Data::Dumper;

#print Dumper \%results;
foreach (0 .. $mol->tmax){
  $mol->gt($_);
  my $chi1 = ceil($dihe1->dihe_deg/$partit);
  my $chi1p = ceil($dihe1p->dihe_deg/$partit);
  my $chi2 = ceil($dihe2->dihe_deg/$partit);
  my $chi2p = ceil($dihe2p->dihe_deg/$partit);
  my $nnlength = ceil(($nnbond->bond_length)/$partit2);
  my $chi3 = ceil($dihe3->dihe_deg/$partit);
  $fresults{$nnlength}{$chi3}++;


  # sets variables to the appropriate angle rounded up and in partitions set above

#  unless( exists( $results{$chi1}{$chi2}{$chi3}{$chi1p}{$chi2p} ) ) {
  unless( exists( $results{$chi3}{$nnlength} ) ) {
    $results{$chi1}{$chi2}++;
    $results{$chi1}{$chi2p}++;
    $results{$chi1p}{$chi2}++;
    $results{$chi1p}{$chi2p}++;



#    $results{$chi3}{$nnlength}++;
  
#  my @energies = $orca2->opt;
#  printf ("%10.3f %10.3f %10.3f %10.3f %10.3f %14.6f\n",$dihe1->dihe_deg, $dihe1p->dihe_deg, $dihe2->dihe_deg, $dihe2p->dihe_deg, $dihe3->dihe_deg, $energies[-1]*627.51);
#  printf $in ("%10.3f %10.3f %10.3f %10.3f %10.3f %14.6f\n",$dihe1->dihe_deg, $dihe1p->dihe_deg, $dihe2->dihe_deg, $dihe2p->dihe_deg, $dihe3->dihe_deg, $energies[-1]*627.51);  
#    printf $in ("%10.3f %10.3f %10.3f\n",$dihe1->dihe_deg, $dihe2->dihe_deg, $nnlength*$partit2); 
#    printf $in ("%10.3f %10.3f %10.3f\n",$dihe1p->dihe_deg, $dihe2p->dihe_deg, $nnlength*$partit2);
  #my $mol2 = HackaMol->new->read_file_mol("tmp/mol.xyz");
  #$mol2->print_xyz($fh);
  
#  else
#  );
#  }
}
#};

my $root = $file;
$root =~ s/\.xyz//;

open(my $fchi3, ">", "$root\_dist_chi3.txt") or die "couldn't open";
open(my $fchi12, ">", "$root\_chi1_chi2.txt") or die "couldn't open";


my $tot= 0;
foreach my $dist (sort { $a <=> $b } keys %fresults){
  foreach my $chi3 (sort{ $a <=> $b } keys %{$fresults{$dist}}){
    $tot++;
    say "total configs: $tot ", "number of configs for $dist $chi3:", $fresults{$dist}{$chi3}; 
  }
}

say "total number of configs ", $tot;

#use Data::Dumper;

#print Dumper \%results;

