use HackaMol;
use HackaMol::X::Orca;
use HackaMol::X::NERF;
use Modern::Perl;
use YAML::XS qw(LoadFile DumpFile);
use POSIX;

my $file = shift;
my $mol  = HackaMol->new->read_file_mol($file);
my %distmat;
my $dist = 0;
my @mistakes;

$mol->t(0);
open( my $in, ">", "mistakes.txt" ) or die "couldn't open";

foreach my $t ( 0 .. $mol->tmax ) {
    my $bad =0;
    $mol->t($t);
    my $ssdist = $mol->get_atoms(5)->distance($mol->get_atoms(11));
    if ($ssdist > 2.5){
      print "REALLY BAD time $t\n"; next; 
    }
       foreach my $i ( 0 .. 22 ) {
           foreach my $j ( $i + 1 .. 22 ) {
               $dist = $mol->get_atoms($i)->distance($mol->get_atoms($j));
               $distmat{$i}{$j} = $dist;
               $bad++ if ($dist < 0.9);
            }
       }
    print "BAD time $t\n" if $bad;
}

#DumpFile( 'distmat.txt', \%distmat );

