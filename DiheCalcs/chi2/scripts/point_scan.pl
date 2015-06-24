use HackaMol;
use HackaMol::X::Orca;
use Modern::Perl;

my $bldr  = HackaMol->new;
my $mol   = $bldr->read_file_mol(shift);

my $dihe = HackaMol::Dihedral->new(
                               is_constrained => 1,
                                    atoms=> [
                                $mol->get_atoms(4),
                                $mol->get_atoms(0),
                                $mol->get_atoms(1),
                                $mol->get_atoms(2),
                                             ],
);


say $dihe->dihe_deg;


$mol->push_dihedrals($dihe);


