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

# set up multiplicities and force constants
$dihe1->dihe_fc(2) ;
$dihe1->dihe_mult(3);
$dihe1p->dihe_fc(2) ;
$dihe1p->dihe_mult(3);
$dihe2->dihe_fc(1) ;
$dihe2->dihe_mult(3);
$dihe2p->dihe_fc(1) ;
$dihe2p->dihe_mult(3);


#checks
say $dihe1->dihe_deg;
say $dihe2->dihe_deg;
say $dihe1p->dihe_deg;
say $dihe2p->dihe_deg;
say $dihe3->dihe_deg;


# get ready for the foreach loop
$mol->gt(0);
open(my $in, ">", "plotly.txt") or die "couldn't open";

$dihe3->torsion_efunc(\&chi3_torsion_efunc);
my @dihedrals = ($dihe1, $dihe1p, $dihe2, $dihe2p, $dihe3);
my $DSE = 0;
$DSE += $_->torsion_energy foreach @dihedrals;
say $DSE;


#foreach (0 .. $mol->tmax){
#  $mol->gt($_);
#  my $ener1 = $dihe1->torsion_energy;
#  my $ener1p = $dihe1p->torsion_energy;
#  my $ener2 = $dihe2->torsion_energy;
#  my $ener2p = $dihe2p->torsion_energy; 
#  my $ener3 = $dihe3->torsion_energy; 


#  my $chi3 = $dihe3->dihe_deg;
#  my $ener3 = 3.5*(1+cos(2*$chi3))+0.6*(1+cos(3*$chi3));



 # my $DSE = $ener1+$ener1p+$ener2+$ener2p+$ener3;
 # printf ("%10.2f\n", $DSE);
 # printf $in ("%10.2f\n", $DSE);
  
#}

sub chi3_torsion_efunc {
  my $dihedral = shift;
  my $ang = $dihedral->dihe_rad;
  my $energy = 3.5*(1+cos(2*$ang))+0.6*(1+cos(3*$ang));
  return ($energy);
}



