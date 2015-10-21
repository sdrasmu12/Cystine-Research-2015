use HackaMol;
use HackaMol::X::Orca;
use HackaMol::X::NERF;
use Modern::Perl;
use YAML::XS qw(Dump LoadFile DumpFile);
use POSIX;

my $file = shift;
my $mol  = HackaMol->new->read_file_mol($file);

my $root = $file;
$root =~s/\/stacks//;
$root =~ s/\.xyz//;
my $histo  = LoadFile("best/nmrs_best.txt");
my %histo  = %{$histo};
my $qm_mol = 0;
my $t      = 0;

open( my $fh,  ">", "sorted/nmrs_bnd.xyz" );
open( my $fh2, ">", "sorted/nmrs_unbnd.xyz" );

foreach my $dist ( sort { $a <=> $b } keys %histo ) {
    foreach my $chi3 ( sort { $a <=> $b } keys %{ $histo{$dist} } ) {
        my $t = $histo{$dist}{$chi3};
        $mol->t($t);

        next if ( $mol->get_atoms(5)->distance( $mol->get_atoms(11)) > 2.5  );

        if ( $mol->get_atoms(2)->distance( $mol->get_atoms(6) ) < 1.8 ) {
            $mol->print_xyz($fh);
        }
        else {
            $mol->print_xyz($fh2);
        }

    }
}

