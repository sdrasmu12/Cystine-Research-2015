use HackaMol;
use HackaMol::X::Orca;
use  HackaMol::X::NERF;
use Modern::Perl;
use YAML::XS qw(LoadFile DumpFile);
use POSIX;

# setup

my $t1 = time;

my $file = shift;
my $mol = HackaMol->new->read_file_mol($file);

say 'Done Reading';

#my $root = $file;
#$root =~ s/\.xyz//;
#my $histo = LoadFile("best\cystine_nmrs_best.txt");
#my %histo = %{$histo};


#establish dihedrals
my $dihe1 = HackaMol::Dihedral->new(
                                is_constrained => 1,
                                    atoms=> [
                                $mol->get_atoms(2),
                                $mol->get_atoms(1),
                                $mol->get_atoms(4),
                                $mol->get_atoms(5),
                                             ],
);


my $dihe2 = HackaMol::Dihedral->new(
                               is_constrained => 1,
                                    atoms=> [
                                $mol->get_atoms(1),
                                $mol->get_atoms(4),
                                $mol->get_atoms(5),
                                $mol->get_atoms(11),
                                             ],
);

my $dihe3 = HackaMol::Dihedral->new(
                                is_constrained => 1,
                                    atoms=> [
                                $mol->get_atoms(4),
                                $mol->get_atoms(5),
                                $mol->get_atoms(11),
                                $mol->get_atoms(10),
                                             ],
);

my $dihe1p = HackaMol::Dihedral->new(
                                is_constrained => 1,
                                   atoms=> [
                                $mol->get_atoms(11),
                                $mol->get_atoms(10),
                                $mol->get_atoms(7),
                                $mol->get_atoms(8),
                                             ],
);


my $dihe2p = HackaMol::Dihedral->new(
                                is_constrained => 1,
                                    atoms=> [
                                $mol->get_atoms(5),
                                $mol->get_atoms(11),
                                $mol->get_atoms(10),
                               $mol->get_atoms(7),
                                             ],
);

say 'established dihes';

# push dihedrals
$mol->push_dihedrals($dihe1);
$mol->push_dihedrals($dihe2);
$mol->push_dihedrals($dihe3);
$mol->push_dihedrals($dihe1p);
$mol->push_dihedrals($dihe2p);
say 'pushed dihes';
my $nnbond = HackaMol::Bond->new(
    atoms=> [
    $mol->get_atoms(6),
   $mol->get_atoms(0),
  ]
);
say 'setup nnbond';

#set up orca
my $orca2 = HackaMol::X::Orca->new(
      mol             => $mol,
      has_constraints => 1,      
      theory          => 'HF-3c',
      exe             => '/Users/chem_student/perl5/apps/orca_3_0_3_macosx_openmpi165/orca',
      scratch         => 'tmp',
);
say 'setup orca';

open(my $in1, ">", "chi3nnlength.txt") or die "couldn't open";
#open(my $in2 ">", "chi3nnlength.txt") or die "couldn't open";
#open(my $in1, ">", "chi3nnlength.txt") or die "couldn't open";

$mol->t(0);

foreach my $t (0..$mol->tmax) {
    $mol->t($t);
    my $dist = $nnbond->bond_length;
    my $chi3 = $dihe3->dihe_deg;
    my @energies = (0);
    @energies = $orca2->opt;
    say foreach @energies;
    my $mol2 = $orca2->load_trj;
    $mol2->print_xyz_ts([0 .. $mol2->tmax]);
    printf ("%10.3f %10.3f %10.2f\n", $dist, $chi3,$mol2->get_energy($mol2->tmax)*627.51);
      printf $in1 ("%10.3f %10.3f %10.2f\n",$dist, $chi3,$energies[-1]*627.51); 
  } 


