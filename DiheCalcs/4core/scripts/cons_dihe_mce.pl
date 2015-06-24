use Modern::Perl;
use HackaMol::X::Orca;
use Math::Vector::Real;
use HackaMol;
use MCE::Loop max_workers => 4, chunk_size => 1;
use MCE::Subs qw( :worker );
use Time::HiRes qw(time);
use Data::Dumper;

my $t1 = time;

my @xyzs = glob("mol_*xyz");

mce_loop_s {
  my $t = $_;
  my $mol = HackaMol->new->read_file_mol("mol_$t.xyz");

  my $dihe = HackaMol::Dihedral->new(
                                is_constrained => 1,
                                    atoms=> [
                                $mol->get_atoms(3),
                                $mol->get_atoms(2),
                                $mol->get_atoms(1),
                                $mol->get_atoms(0),
                                             ],
  ); 

  $mol->push_dihedrals($dihe);
  $mol->push_charges(0);
  $mol->multiplicity(1);
  
  my $orca = HackaMol::X::Orca->new(
      mol             => $mol,
      has_constraints => 1,      
      theory          => 'HF-3c',
      exe             => '/Users/chem_student/perl5/apps/orca_3_0_3_macosx_openmpi165/orca',
      scratch         => "tmp_$t",
  );

  my $ang = $dihe->dihe_deg;
  my @energies = $orca->opt;
  mce_printf ("%10.3f %14.6f\n", $dihe->dihe_deg, $energies[-1]*627.51);
} 0 , $#xyzs;

my $t2 = time;

printf ("time %10.2f\n", $t2-$t1);
