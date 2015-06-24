use HackaMol;
use HackaMol::X::Orca;
use Modern::Perl;

my $bldr  = HackaMol->new;
my $mol   = $bldr->read_file_mol(shift);
#  create bldr class, read file into mol


my @atoms = $mol->all_atoms;
# puts atoms from molecule into the @atoms array



#my @atoms = HackaMol->new->file_read_atoms(shift);




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


my @atoms_rotate = map {$mol->get_atoms($_)}(4,5,6);


my $atoms_rotate =
  HackaMol::AtomGroup->new( atoms => [ @atoms_rotate ] );


##############################################################################
#    Loop over rotations to generate inputs at several dihedral angles       #
#    this can and should be done in steps for larger molecules               #
##############################################################################

my $orca = HackaMol::X::Orca->new(
      mol             => $mol,
      has_constraints => 1,
     # theory          => 'AM1',
      theory          => 'HF-3c',
      exe             => '/Users/chem_student/perl5/apps/orca_3_0_3_macosx_openmpi165/orca',
      scratch         => "tmp",
);

my $dang = 15;
my $ceil = int( 360 / $dang );
#dang is increments, int counts to create the max for the foreachloop 

my $t = 0;
open(my $in, ">", "plotly.txt") or die "couldn't open";
#setup

foreach ( 1 .. $ceil ) {

    $mol->dihedral_rotate_groups( $dihe, $dang ,
        $atoms_rotate );

#rotates by increments established above

    my @energies = $orca->opt;
# optimizes using orca

    $mol->print_xyz('tmp/mol.xyz');
# prints the optimized molecule to the file

    $bldr->read_file_push_coords_mol("tmp/mol.xyz",$mol); 
# reads the optimized molecule back into $mol

    $t++;
    $mol->t($t); 
# increments t and changes t
    
   printf $in ("%10.3f %14.6f\n", $dihe->dihe_deg, $energies[-1]*627.51);
   #printf ("%10.3f %14.6f\n", $dihe->dihe_deg, $energies[-1]*627.51);
   $mol->print_xyz;
# prints energies to plotly and STDOUT    

}


sub qrotatable {
    my $atoms   = shift;
    my $iroot   = shift;
    my $visited = shift;

# sets these variables to the three arguments written 

    $visited->{$iroot}++;

# visited points to iroot, not sure how the incrementing works

    my @cands;
    foreach my $at (@$atoms) {
        push @cands, $at unless ( grep { $at->iatom == $_ } keys %{$visited} );
    }
# each atom is pushed to cands unless..


    #find S-C bonds
    my @bonds = $bldr->find_bonds_brute(
        bond_atoms => [ $atoms->[$iroot] ],
        candidates => [@cands],
        fudge      => 0.45,
    );

    foreach my $cand ( map { $_->get_atoms(1) } @bonds ) {
        next if $visited->{ $cand->iatom };
        my $visited = qrotatable( $atoms, $cand->iatom, $visited );
    }
    return ($visited);
}

