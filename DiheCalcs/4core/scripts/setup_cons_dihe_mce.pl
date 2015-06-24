use Modern::Perl;
use HackaMol::X::Orca;
use Math::Vector::Real;
use HackaMol;
use Time::HiRes qw(time);
use Data::Dumper;

my $t1 = time;

my $mol = HackaMol->new->read_file_mol(shift);

foreach my $t (0 .. $mol->tmax){
  $mol->gt($t); 
  $mol->print_xyz("mol_$t.xyz");
}

exit;
foreach my $t (0 .. $mol->tmax){
  # because we want to run in parallel, we must be sure we pull the correct state of the molecule
  # so we will write out the state and then read into a new molecule
  my $mol2 = HackaMol->new->read_file_mol("mol_$t.xyz");

  my $dihe2 = HackaMol::Dihedral->new(
                                is_constrained => 1,
                                    atoms=> [
                                $mol2->get_atoms(3),
                                $mol2->get_atoms(2),
                                $mol2->get_atoms(1),
                                $mol2->get_atoms(0),
                                             ],
  ); 

  $mol2->push_dihedrals($dihe2);
  $mol2->push_charges(0);
  $mol2->multiplicity(1);
  
  my $orca = HackaMol::X::Orca->new(
      mol             => $mol2,
      has_constraints => 1,      
      theory          => 'HF-3c',
      exe             => '/scratch/riccade/bin/orca_3_0_3_linux_x86-64/orca',
      scratch         => "tmp_$t",
  );

  my $ang = $dihe2->dihe_deg;
  my @energies = $orca->opt;

  mce_printf ("%10.3f %14.6f\n", $dihe2->dihe_deg, $energies[-1]*627.51);
}

my $t2 = time;

printf ("time %10.2f\n", $t2-$t1);
