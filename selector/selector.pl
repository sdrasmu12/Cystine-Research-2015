use HackaMol;
use Modern::Perl;
use YAML::XS qw(LoadFile DumpFile);
use POSIX;

# setup

my $t1 = time;

my $file = shift;
my $mol = HackaMol->new->read_file_mol($file);
say 'Done Reading';

my $root = $file;
$root =~ s/\.xyz//;
say $root;
my $histo = LoadFile("stacks/cystine_nmrs_dist_chi3.yaml");
my %histo = %{$histo};


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

#set up mult and fc
$dihe1->dihe_fc(2) ;
$dihe1->dihe_mult(3);
$dihe1p->dihe_fc(2) ;
$dihe1p->dihe_mult(3);
$dihe2->dihe_fc(1) ;
$dihe2->dihe_mult(3);
$dihe2p->dihe_fc(1) ;
$dihe2p->dihe_mult(3);

my %best;
my $tot= 0;
foreach my $dist (sort { $a <=> $b } keys %histo){
  foreach my $chi3 (sort{ $a <=> $b } keys %{$histo{$dist}}){
    my $dse  = 1800;
    my $ldse = 1800;
    $best{$dist}{$chi3} = 9999 unless (exists($best{$dist}{$chi3}));
    foreach my $t (@{$histo{$dist}{$chi3}}){
      #my $t = $histo{$dist}{$chi3}[0];
      $mol->t($t);
      my $ldse = $dihe1->torsion_energy + $dihe1p->torsion_energy + $dihe2->torsion_energy +$dihe2p->torsion_energy; 
      say "shitt $t shitdse $dse  shitldse $ldse shitdist $dist shitchi3 $chi3"."shitbest: ".$best{$dist}{$chi3};
      if ($ldse < $dse){
        $best{$dist}{$chi3}= $t;
        $dse = $ldse;
      } 
  }
}
}
DumpFile("best/nmrs_best.txt", \%best);


