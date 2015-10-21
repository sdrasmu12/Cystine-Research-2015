use HackaMol;
use HackaMol::X::Orca;
use HackaMol::X::NERF;
use Modern::Perl;
use YAML::XS qw(Dump LoadFile DumpFile);
use POSIX;

# setup


my $file = shift;
my $mol = HackaMol->new->read_file_mol($file);
$mol->t(0);

my $qm_mol = 0;

foreach my $t (0..$mol->tmax){
      $mol->t($t);
         if ($mol->get_atoms(2)->distance ($mol->get_atoms(6)) < 1.8 ){
                 $qm_mol = bndcystine_protonate($mol, $t);
                  }else{
                  $qm_mol = cystine_protonate($mol, $t);
                  };
   $qm_mol->print_xyz;
  }

sub bndcystine_protonate {
  $mol = shift;
  my $t   = shift;

  $mol->t($t);


  my $nerf = HackaMol::X::NERF->new;
  my @xyzs = map{$_->xyz} $mol->all_atoms;
  my $hmol = HackaMol::Molecule->new(atoms=>[$mol->all_atoms]);

  @xyzs = map{$_->xyz} $hmol->all_atoms;
  my $h1  = $nerf->extend_abc( $xyzs[2],  $xyzs[1], $xyzs[0],  1.0,  120, 0 );
  my $h2  = $nerf->extend_abc( $xyzs[2],  $xyzs[1], $xyzs[0],  1.0,  120, 180 );
  my $h3  = $nerf->extend_abc( $xyzs[7],  $xyzs[2], $xyzs[6],  1.0,  109,  180);
 # my $h4  = $nerf->extend_abc( $xyzs[7],  $xyzs[2], $xyzs[6],  1.0,  120, -120);
  my $h5  = $nerf->extend_abc( $xyzs[2],  $xyzs[0], $xyzs[1],  1.0,  109, -120);
 my $h6  = $nerf->extend_abc( $xyzs[6],  $xyzs[8], $xyzs[7],  1.0,  109, -120);
  my $h7  = $nerf->extend_abc( $xyzs[1],  $xyzs[5], $xyzs[4],  1.0,  109, -120);
  my $h8  = $nerf->extend_abc( $xyzs[1],  $xyzs[5], $xyzs[4],  1.0,  109, 120);
 my $h9  = $nerf->extend_abc( $xyzs[9],  $xyzs[7], $xyzs[8],  1.0,  109, 180);
 my $h10  = $nerf->extend_abc( $xyzs[7], $xyzs[11], $xyzs[10], 1.0,  109, 120);
 my $h11  = $nerf->extend_abc( $xyzs[7], $xyzs[11], $xyzs[10], 1.0,  109, -120); 

 if ($h5->dist($hmol->get_atoms(4)->xyz) < 0.9 ){ 
     $h5  = $nerf->extend_abc( $xyzs[2],  $xyzs[0], $xyzs[1],  1.0,  109, 120);
 }

 if ($h6->dist ($hmol->get_atoms(10)->xyz) < 0.9){
   $h6  = $nerf->extend_abc( $xyzs[6],  $xyzs[8], $xyzs[7],  1.0,  109, 120)
 }
  $hmol->push_atoms(
    map{ HackaMol::Atom->new(Z => 1, coords => [$_]) } ($h1,$h2,$h3,$h5,$h6,$h7,$h8,$h9,$h10,$h11)
 );

#$hmol->print_xyz;

return ($hmol);

}


sub cystine_protonate {
  my $mol = shift;
  my $t   = shift;

  $mol->t($t);

  my $nerf = HackaMol::X::NERF->new;
  #build in two oxygen
  my @xyzs = map{$_->xyz} $mol->all_atoms;
  use Data::Dumper;

  my $o1  = $nerf->extend_abc( $xyzs[3],  $xyzs[1], $xyzs[2],  1.3,  120, 180 );
  my $o2  = $nerf->extend_abc( $xyzs[9],  $xyzs[7], $xyzs[8],  1.3,  120, 180 );
  my $hmol = HackaMol::Molecule->new(atoms=>[$mol->all_atoms,
        map{ HackaMol::Atom->new(Z => 8, coords => [$_]) } ($o1,$o2)
      ]
  );

  @xyzs = map{$_->xyz} $hmol->all_atoms;
  my $h1  = $nerf->extend_abc( $xyzs[2],  $xyzs[1], $xyzs[0],  1.0,  120, 0 );
  my $h2  = $nerf->extend_abc( $xyzs[2],  $xyzs[1], $xyzs[0],  1.0,  120, 180 );
  my $h3  = $nerf->extend_abc( $xyzs[8],  $xyzs[7], $xyzs[6],  1.0,  120, 1 );
  my $h4  = $nerf->extend_abc( $xyzs[8],  $xyzs[7], $xyzs[6],  1.0,  120, 180 );
  my $h5  = $nerf->extend_abc( $xyzs[2],  $xyzs[0], $xyzs[1],  1.0,  109, -120);
  my $h6  = $nerf->extend_abc( $xyzs[6],  $xyzs[8], $xyzs[7],  1.0,  109, -120);
  my $h7  = $nerf->extend_abc( $xyzs[1],  $xyzs[5], $xyzs[4],  1.0,  109, -120);
  my $h8  = $nerf->extend_abc( $xyzs[1],  $xyzs[5], $xyzs[4],  1.0,  109, 120);
  my $h9  = $nerf->extend_abc( $xyzs[7], $xyzs[11], $xyzs[10], 1.0,  109, 120);
  my $h10 = $nerf->extend_abc( $xyzs[7], $xyzs[11], $xyzs[10], 1.0,  109, -120);
  my $h11 = $nerf->extend_abc( $xyzs[3],  $xyzs[2], $xyzs[12], 1.0,  109, 180);
  my $h12 = $nerf->extend_abc( $xyzs[9],  $xyzs[8], $xyzs[13], 1.0,  109, 54);



 if ($h5->dist($hmol->get_atoms(4)->xyz) < 0.9 ){
     $h5  = $nerf->extend_abc( $xyzs[2],  $xyzs[0], $xyzs[1],  1.0,  109, 120);
 }

 if ($h6->dist ($hmol->get_atoms(10)->xyz) <0.9){
   $h6  = $nerf->extend_abc( $xyzs[6],  $xyzs[8], $xyzs[7],  1.0,  109, 120)
 }

  $hmol->push_atoms(
    map{ HackaMol::Atom->new(Z => 1, coords => [$_]) } ($h1,$h2,$h3,$h4,$h5,$h6,$h7,$h8,$h9,$h10,$h11,$h12)
 );

#$hmol->print_xyz;
  
return ($hmol);
}


