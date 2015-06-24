use HackaMol;
use Modern::Perl;
use Math::Trig;


open(my $in1, ">", "chi1vchi3.txt") or die "couldn't open";
my $chi1 = 0;
my $chi2 = 0;
my $chi3 = 0;
my $ener = 0;
my $DSE = 0;

foreach my $chi1 (-180..180){
  #my $chi1 = $_ ;
  my $DSE1 = 2*(1+cos(3*deg2rad($chi1)));#+2*(1+cos(3*deg2rad($chi1)));
  foreach my $chi3 (-180..180){
                 #my $chi_3 = $_ ;
                 my $DSE3 =  3.5*(1+cos(2*deg2rad($chi3)))
                             +0.6*(1+cos(3*deg2rad($chi3)));
                 $DSE = $DSE3 + $DSE1;  
                 printf $in1 ("%10.2f %10.2f %10.2f\n", $chi1, $chi3, $DSE);
                 printf ("%10.2f\n", $DSE);
    };

};


open(my $in2, ">", "chi1vchi2.txt") or die "couldn't open";


foreach my $chi1 (-180..180){
  #my $chi1 = $_ ;
  my $DSE1 = 2*(1+cos(3*deg2rad($chi1)));#+2*(1+cos(3*deg2rad($chi1)));
  foreach my $chi2 (-180..180){
                 #my $chi_3 = $_ ;
                 my $DSE2 =  (1+cos(3*deg2rad($chi2)));
                 $DSE = $DSE2 + $DSE1;
                 printf $in2 ("%10.2f %10.2f %10.2f\n", $chi1, $chi2, $DSE);
                 printf ("%10.2f\n", $DSE);
    };

};

open(my $in3, ">", "chi3vchi2.txt") or die "couldn't open";


foreach my $chi3 (-180..180){
  #my $chi1 = $_ ;
  my $DSE3 = 3.5*(1+cos(2*deg2rad($chi3)))+0.6*(1+cos(3*deg2rad($chi3)));
  foreach my $chi2 (-180..180){
                 #my $chi_3 = $_ ;
                 my $DSE2 =  (1+cos(3*deg2rad($chi2)));
                 $DSE = $DSE2 + $DSE3;
                 printf $in3 ("%10.2f %10.2f %10.2f\n", $chi3, $chi2, $DSE);
                 printf ("%10.2f\n", $DSE);
    };

};

  
#  my $DSE1 = 2*(1+cos(3*$chi1))+2*(1+cos(3*$chi1p));
# my $DSE2 = (1+cos(3*$chi2))+(1+cos(3*$chi2p));
#  my $DSE3 = 3.5*(1+cos(2*$chi3))+0.6*(1+cos(3*$chi3));
#  my $DSE = $DSE1+$DSE2+$DSE3;
#  printf ("%10.2f\n", $DSE);
#  printf $in ("%10.2f\n", $DSE);
  
