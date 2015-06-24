use Modern::Perl;
use HackaMol::X::Orca;
use Math::Vector::Real;
use HackaMol;
use Time::HiRes qw(time);
use Data::Dumper;

my $t1 = time;
# not sure

my $mol = HackaMol->new->read_file_mol(shift);
# read into $mol

$mol->push_charges(0);
$mol->multiplicity(1);
#stuff we need


my $dihe = HackaMol::Dihedral->new(
                                is_constrained => 1,
                                    atoms=> [
                                $mol->get_atoms(3),
                                $mol->get_atoms(2),
                                $mol->get_atoms(1), 
                                $mol->get_atoms(0),
                                             ],
);
#creates the dihedral object

$mol->push_dihedrals($dihe);
#pushes the dihedral object into the molecule

my $orca2 = HackaMol::X::Orca->new(
      mol             => $mol,
#      has_constraints => 1,      
      theory          => 'HF-3c',
      exe             => '/Users/chem_student/perl5/apps/orca_3_0_3_macosx_openmpi165/orca',
      scratch         => 'tmp',
);
#sets up orca

#$mol->gt(0);
# set time 0 for the molecule I don't think we will use t

my $fh = $mol->print_xyz('shit.xyz');
#prints the molecule to shit.xyz
open(my $in, ">", "plotly.txt") or die "couldn't open";
# opens plotly.txt for printing
my $min = 9999;
#sets my min 

foreach (0 .. $mol->tmax){
  $mol->gt($_);
  my @energies = $orca2->opt;
  my $e_rel = ($energies[-1]-$min)*627.51;
  printf ("%10.3f %10.2f %14.6f\n", $dihe->dihe_deg, $e_rel, $energies[-1]*627.51);
  printf $in ("%10.3f %10.2f %14.6f\n", $dihe->dihe_deg, $e_rel, $energies[-1]*627.51);
  my $mol2 = HackaMol->new->read_file_mol("tmp/mol.xyz");
  $mol2->print_xyz($fh);
}

#$orca->map_input;
#$orca2->load_engrad;

my $t2 = time;





printf ("%10.2f\n", $t2-$t1):
