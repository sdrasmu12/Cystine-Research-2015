use Modern::Perl;
use HackaMol;
use POSIX;
my $bldr = HackaMol->new; 1;


my $cystine_xtal = $bldr->read_file_mol('cystine_nmrs.xyz');


my $chi3  = HackaMol::Dihedral->new(atoms => [ map { $cystine_xtal->get_atoms($_)} 4,5,11,10 ]);
my $chi1  = HackaMol::Dihedral->new(atoms => [ map { $cystine_xtal->get_atoms($_)} 1,4,5,11]);
my $chi1p = HackaMol::Dihedral->new(atoms => [ map { $cystine_xtal->get_atoms($_)} 7,10,11,5]);
my $chi2 = HackaMol::Dihedral->new(atoms => [ map { $cystine_xtal->get_atoms($_)} 5,11,10,7]);
my $chi2p = HackaMol::Dihedral->new(atoms => [ map { $cystine_xtal->get_atoms($_)} 11,5,4,1]);


my %chi13_bin = ();
#my $gran = 1;
my $total;

foreach my $t (0 .. $cystine_xtal->tmax ){ #$cystine_xtal->tmax){
    $cystine_xtal->t($t);
#    my $dchi3  = floor($chi3->dihe_deg/$gran);
#    my $dchi1 = floor($chi1->dihe_deg/$gran);
#    my $dchi1p = floor($chi1p->dihe_deg/$gran);
    #say "$dchi3 $dchi1 $dchi1p";
    my $dchi3  = $chi3->dihe_deg;
    my $dchi1  = $chi1->dihe_deg;
    my $dchi1p = $chi1p->dihe_deg;
    $chi13_bin{$dchi3}{$dchi1}++;
    $chi13_bin{$dchi3}{$dchi1p}++;
    $total += 2;
}

say "1";

open (my $fh , ">", 'nmrs_output.txt');
printf $fh ("%10s %10s %10s\n", 'x','y','z');
foreach my $chi3 (sort {$a<=>$b} keys %chi13_bin){
    foreach my $chi1 (sort {$a<=>$b} keys %{$chi13_bin{$chi3}}){
        printf $fh ("%10.2f %10.2f %10.2f\n",$chi3,$chi1,$chi13_bin{$chi3}{$chi1} ); 
       # printf $fh ("%10.2f %10.2f %10.2f\n",$chi3*$gran,$chi1*$gran,100*$chi13_bin{$chi3}{$chi1}/$total); 


    }
}
close($fh);
