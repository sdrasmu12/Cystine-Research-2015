use HackaMol;
use Modern::Perl;
use YAML::XS qw(DumpFile);
use POSIX;

# setup

my $t1 = time;

my $file = shift;
my $mol = HackaMol->new->read_file_mol($file);

say 'Done Reading';

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
                                $mol->get_atoms(8),
                                $mol->get_atoms(7),
                                $mol->get_atoms(10),
                                $mol->get_atoms(11),
                                             ],
);


my $dihe2p = HackaMol::Dihedral->new(
                                is_constrained => 1,
                                    atoms=> [
                                $mol->get_atoms(7),
                                $mol->get_atoms(10),
                                $mol->get_atoms(11),
                               $mol->get_atoms(5),
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
my %fresults;
my $partit = 5;
my $partit2 = 1;

my $nnbond = HackaMol::Bond->new(
    atoms=> [
    $mol->get_atoms(6),
   $mol->get_atoms(0),
  ]
);


say $nnbond->bond_length;
my $nnlength = 0;

# step 1 bin by distance and chi3

foreach my $t (0 .. $mol->tmax){
  $mol->gt($t);
  my $nnlength = ceil(($nnbond->bond_length)/$partit2);
  my $chi3 = ceil($dihe3->dihe_deg/$partit);
  push @{$fresults{$nnlength}{$chi3}}, $t;
}


my $root = $file;
$root =~ s/\.xyz//;

DumpFile("$root\_dist_chi3.yaml", \%fresults);

