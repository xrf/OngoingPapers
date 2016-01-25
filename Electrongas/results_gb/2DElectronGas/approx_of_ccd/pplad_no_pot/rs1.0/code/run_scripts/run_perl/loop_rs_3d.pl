
#
#     This script loops over different values for the
#     number of particles. The script writes an input
#     file and runs a calculation for every given
#     particle number and number of shells
#
#     Run in background:      perl loop_nparticles.pl &
#
use strict;
use warnings;

#
#     Input parameters for the code:
#

# Harmonic oscillator strength
#
my $omega = 1.0;
# Mean distance between particles rs
#
my $rs;
my @rsArray = (0.5, 1.0, 2.0);
# Number of particles and total number of shells	  
# (nOccupied, nShells)
# 
my $nOccupied = 14;
my $nShells = 5;
# Single-particle basis
#
#    ho2d: Harmonic oscillator in two dimensions
#    pw2d: Plane wave basis for two-dimensional 
#          periodic systems
#    pw3d: Plane wave basis for three-dimensional
#          periodic systems
#
my $spBasis = 'pw3d';
# Spin state
#
#    normal: Spin fluid with total spin zero
#    polarized: Spin polarized fluid 
#
my $spinState = 'normal';
# Two-body interaction
# 
#    yukawa2d: Yukawa interaction for two-dimensional
#              periodic systems with plane wave basis
#    yukawa3d: Yukawa interaction for three-dimensional
#              periodic systems with plane wave basis
#    coulombHO2d: Coulomb interaction in two-dimensional
#                 harmonic oscillator basis 
#
my $interaction = 'yukawa3d';
# Screening factor mu in the Yukawa interaction
#
my $mu = 0.0;
# Approximation
#         ref: Only the reference energy
#         pt2: Second-order perturbation theory
#         ccd: Coupled-cluster CCD approximation
my $approx = 'ccd';
# tolerance:   Tolerance in the CC selfconsistency
#              loop (MeV)
#
my $tolerance = 0.00000001;
# maxIter:     Maximum number of iterations
#              in the CC selfconsistency loop
my $maxIter = 100;

my $rsEne = 'rs_ene';
my $opened = open my $fileHandle1, ">", $rsEne;
if(!$opened) {
    die "Could not open the input file";
}

my @lines = (
    "#\n",
    "# Energy per particle of three-dimensional electron gas\n",
    "#\n",
    "# * 14 electrons and three energy shells, mu = 0.0\n",
    "# * Normal spin fluid (not polarized)\n",
    "#\n",
    "# Energies in units of rydberg (e^2/(2*a_{0}))\n",
    "#\n",
    "# ene_ref = Reference energy per particle \n",
    "# ene_corr = CCD correlation energy per particle\n",
    "# ene_tot = ene_ref + ene_corr\n",
    "#\n",
    "# rs	ene_ref	      ene_corr     ene_tot\n",
    );

#     Write the lines to file
my $line = 0; 
while($line < scalar @lines) {
    print $fileHandle1 $lines[$line];
    $line++;
}
close $fileHandle1;

my $status = system("cp ".$rsEne." rs_ene_old");

#     Loop over elements in @rsArray.
#
my $count = 0;
while($count < scalar @rsArray) {
    
    $rs = $rsArray[$count];
    
    #     Open input file
    #
    my $inFile = 'input_rs'.$rs.'.dat';
    my $opened = open my $fileHandle, ">", $inFile; 
    if(!$opened) {
	die "Could not open the input file";
    }
    
    #     Array for the lines of the new input file
    #
    my @lines = (
	"# Harmonic oscillator strength\n",
	"#\n",
	$omega, "\n", 
	"# Mean distance between particles rs\n",
	"#\n",
	$rs, "\n",
	"# Number of particles and total number of shells\n",
	"# (nOccupied, nShells)\n",
	"#\n",
	$nOccupied, ", ", $nShells, "\n",
	"# Single-particle basis\n",
	"#\n",
	"#    ho2d: Harmonic oscillator in two dimensions\n",
	"#    pw2d: Plane wave basis for two-dimensional\n",
	"#          periodic systems\n",
	"#    pw3d: Plane wave basis for three-dimensional\n",
	"#          periodic systems\n",
	"#\n",
	$spBasis, "\n",
	"# Spin state\n",
	"#\n",
	"#    normal: Spin fluid with total spin zero\n",
	"#    polarized: Spin polarized fluid \n",
	"#\n",
	$spinState, "\n",
	"# Two-body interaction\n",
	"#\n",
	"#    yukawa2d: Yukawa interaction for two-dimensional\n",
	"#              periodic systems with plane wave basis\n",
	"#    yukawa3d: Yukawa interaction for three-dimensional\n",
	"#              periodic systems with plane wave basis\n",
	"#    coulombHO2d: Coulomb interaction in two-dimensional\n",
	"#                 harmonic oscillator basis\n",
	"#\n",
	$interaction, "\n",
	"# Screening factor mu in the Yukawa interaction\n",
	"#\n",	
	$mu, "\n",
	"# Approximation\n",
	"#         ref: Only the reference energy\n",
	"#         pt2: Second-order perturbation theory\n",
	"#         ccd: Coupled-cluster CCD approximation\n",
	$approx, "\n",
	"# tolerance:   Tolerance in the CC selfconsistency\n",
	"#              loop (MeV)\n",
	"#\n",
	$tolerance, "\n",
	"# maxIter:     Maximum number of iterations\n",
	"#              in the CC selfconsistency loop\n",
	$maxIter, "\n",
	);
    
    #     Write the lines to file
    my $line = 0; 
    while($line < scalar @lines) {
	print $fileHandle $lines[$line];
	$line++;
    }
    
    #     Run a simulation with the produced input file
    my $statusCp = system("cp ".$inFile." input.dat");
    my $outFile = "outlist_rs".$rs;
    my $statusCcd = system("./ccd.exe input.dat > ".$outFile." < /dev/null ");

    $opened = open my $fileHandleRs, "> rs";
    $opened = open my $fileHandleRsNew, "> rs_ene_new";
    print $fileHandleRs $rs;
    my $pasteCmd = "paste rs ene_ref_corr_tot | column -t > rs_ene_new";
    my $status = system($pasteCmd);
    
    $status = system("cat rs_ene_old rs_ene_new > ".$rsEne);
    $status = system("cp ".$rsEne." rs_ene_old");
    
    close $fileHandle;
    $count++;
    
}

$status = system("rm rs_ene_old rs_ene_new rs");
