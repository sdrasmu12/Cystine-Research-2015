use HackaMol;
use HackaMol::X::Orca;
use Modern::Perl;

my $bldr  = HackaMol->new;
my $mol   = $bldr->read_file_mol(shift);
#  create bldr class, read file into mol


my @atoms = $mol->all_atoms;
# puts atoms from molecule into the @atoms array



#my @atoms = HackaMol->new->file_read_atoms(shift);

my @Ss = grep { $_->symbol eq 'S' } @atoms;
my @Cs = grep { $_->symbol eq 'C' } @atoms;
# sorts atoms into carbons and sulfurs



my @SCs = $bldr->find_bonds_brute(
    bond_atoms => [@Ss],
    candidates => [@Cs],
    fudge      => 0.45,
);
# finds S-C bonds 


#bond_atoms (S) are first in the group! wanted: C-S -- S-C
my ($dihe) =
  $bldr->build_dihedrals( reverse( $SCs[0]->all_atoms ), $SCs[1]->all_atoms );
# builds dihedrals out of reverse order of all atoms in the first and second SC bond
# note: build_dihedrals can make dihedrals out of more than just 4 inputs, if 5 it wild make two of 1234 and 2345.
 
say $dihe->dihe_deg; 

$dihe->is_constrained(1);
$mol->push_dihedrals($dihe);
#constrains the dihedral


my $init = {
    $dihe->get_atoms(1)->iatom => 1,
    $dihe->get_atoms(2)->iatom => 1,
};
# I'm not sure what iatom is

my $atoms_rotate = qrotatable( $mol->atoms, $dihe->get_atoms(2)->iatom, $init );
delete( $atoms_rotate->{ $Ss[0]->iatom } );
delete( $atoms_rotate->{ $Ss[1]->iatom } );

my $group_rotate =
  HackaMol::AtomGroup->new( atoms => [ @atoms[ keys %{$atoms_rotate} ] ] );


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
#set up orca

#my $dang = 10;
#my $ceil = int( 180 / $dang );
#foreach my $ang ( map { $dang * $_ } 0 .. $ceil ) {
#    $mol->dihedral_rotate_groups( $dihe, $dihe->dihe_deg - $ang,
#        $group_rotate );
#    $mol->print_xyz;
#}

#exit;

my $dang = 10;
my $ceil = int( 180 / $dang );
#dang is increments, int counts to create the max for the foreachloop 

my $t = 0;
open(my $in, ">", "plotly.txt") or die "couldn't open";
#setup

foreach ( 1 .. $ceil ) {

    $mol->dihedral_rotate_groups( $dihe, $dang ,
        $group_rotate );

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

