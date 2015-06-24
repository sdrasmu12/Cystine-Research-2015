use HackaMol;
use HackaMol::X::Orca;
use Modern::Perl;
use POSIX;

# setup

my $t1 = time;

my $mol = HackaMol->new->read_file_mol(shift);

$mol->push_charges(0);
$mol->multiplicity(1);


#establish dihedrals

my $dihe1 = HackaMol::Dihedral->new(
                                is_constrained => 1,
                                    atoms=> [
                                $mol->get_atoms(1),
                                $mol->get_atoms(0),
                                $mol->get_atoms(3),
                                $mol->get_atoms(4),
                                             ],
);

my $dihe2 = HackaMol::Dihedral->new(
                               is_constrained => 1,
                                    atoms=> [
                                $mol->get_atoms(0),
                                $mol->get_atoms(3),
                                $mol->get_atoms(4),
                                $mol->get_atoms(10),
                                             ],
);

my $dihe3 = HackaMol::Dihedral->new(
                                is_constrained => 1,
                                    atoms=> [
                                $mol->get_atoms(3),
                                $mol->get_atoms(4),
                                $mol->get_atoms(10),
                                $mol->get_atoms(9),
                                             ],
);

my $dihe1p = HackaMol::Dihedral->new(
                                is_constrained => 1,
                                   atoms=> [
                                $mol->get_atoms(7),
                                $mol->get_atoms(6),
                                $mol->get_atoms(9),
                                $mol->get_atoms(10),
                                             ],
);


my $dihe2p = HackaMol::Dihedral->new(
                                is_constrained => 1,
                                    atoms=> [
                                $mol->get_atoms(7),
                                $mol->get_atoms(9),
                                $mol->get_atoms(10),
                               $mol->get_atoms(4),
                                             ],
);



# push dihedrals
$mol->push_dihedrals($dihe1);
$mol->push_dihedrals($dihe2);
$mol->push_dihedrals($dihe3);
$mol->push_dihedrals($dihe1p);
$mol->push_dihedrals($dihe2p);

#checks
say $dihe1->dihe_deg;
say $dihe2->dihe_deg;
say $dihe1p->dihe_deg;
say $dihe2p->dihe_deg;
say $dihe3->dihe_deg;


# get ready for the foreach loop
$mol->gt(0);
open(my $in, ">", "plotly.txt") or die "couldn't open";

foreach (0 .. $mol->tmax){
  $mol->gt($_);
  my $chi1 = $dihe1->dihe_deg;
  my $chi1p = $dihe1p->dihe_deg;
  my $chi2 = $dihe2->dihe_deg;
  my $chi2p = $dihe2p->dihe_deg;
  my $chi3 = $dihe3->dihe_deg;
  
  my $DSE1 = 2*(1+cos(3*$chi1))+2*(1+cos(3*$chi1p));
  my $DSE2 = (1+cos(3*$chi2))+(1+cos(3*$chi2p));
  my $DSE3 = 3.5*(1+cos(2*$chi3))+0.6*(1+cos(3*$chi3));
  my $DSE = $DSE1+$DSE2+$DSE3;
  printf ("%10.2f\n", $DSE);
  printf $in ("%10.2f\n", $DSE);
  

}

my $t2 = time;

printf ("%10.2f\n", ($t2-$t1)/60);