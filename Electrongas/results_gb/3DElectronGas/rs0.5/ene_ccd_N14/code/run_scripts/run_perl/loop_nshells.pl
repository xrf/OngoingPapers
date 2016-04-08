
#
#     This script loops over different values for the
#     number of single-particle energy shells. The script 
#     writes an input file and runs a calculation for 
#     every given number of shells.
#
#     Run in background:      perl loop_nshells.pl &
#
use strict;
use warnings;

#
#     Input parameters for the code:
#

# Mean distance between particles rs
#
my $rs = 0.5;
# Number of particles and total number of shells	  
# (nOccupied, nShells)
# 
my $nOccupied = 14;
my $nShells;
# Single-particle basis
#
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
my $maxIter = 140;
# alpha:       t = alpha*t + (1-alpha)*t_old
#              alpha \in [0.0, 1.0]
my $alpha = 1.0;

# my @nShellsArray = (2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 
#  		    12, 13, 14, 15, 16, 17, 18, 19, 20, 
#  		    21, 22, 23, 24, 25, 26, 27, 28, 29,
# 		    30, 31, 32, 33, 34, 35, 36, 37, 38, 39,
# 		    40, 53, 58, 60, 70, 80, 100, 120, 200);
#my @nShellsArray = (4, 6, 8, 10, 20, 30, 40, 60, 80, 100, 120, 
#		    140, 160, 180, 200, 240, 280, 340, 400, 500);
my @nShellsArray = (36);
#my @nShellsArray = (120, 140, 160, 180, 200, 220, 240, 260, 280);
# my @nShellsArray = (6, 8, 10, 12, 14, 16, 18, 20, 
# 		    22, 24, 26, 28, 30, 32, 34, 36, 38, 40, 42,
# 		    50, 60, 80, 100, 120, 140, 160, 180, 200,
# 		    220, 240, 260, 280, 300);
#my @nShellsArray = (4, 6, 8, 10, 12, 14, 16, 18, 19, 20, 22, 24, 
#		    32, 37);


my $NEneTot = 'nShellsEneTot';
my $opened = open my $fileHandle1, ">", $NEneTot;
if(!$opened) {
    die "Could not open the input file";
}


my @lines = (
    "#\n",
    "# CC energies of three-dimensional electron gas\n",
    "#\n",
    "# rs = 0.5, N = 14, alpha = 1.0\n",
    "# nShells = Number of sp energy shells\n",
    "# ene_ref = Reference energy per particle in units\n",
    "#           of rydberg (e^2/(2*a_{0})) \n",
    "# ene_corr = CCD correlation energy\n",
    "#\n",
    "# nShells nOrbits ene_ref ene_corr\n",
    );

#     Write the lines to file
my $line = 0; 
while($line < scalar @lines) {
    print $fileHandle1 $lines[$line];
    $line++;
}
close $fileHandle1;

my $status = system("cp ".$NEneTot." nshells_ene_old");


#     Loop over elements in @nShellsArray.
#
my $count = 0;
while($count < scalar @nShellsArray) {

    $nShells = $nShellsArray[$count];
    
    #     Open input file
    #
    my $inFile = 'input_nShells'.$nShells.'.dat';
    my $opened = open my $fileHandle, ">", $inFile; 
    if(!$opened) {
	die "Could not open the input file";
    }

        
    #     Array for the lines of the new input file
    #
    my @lines = ( 
	"# Mean distance between particles rs\n",
	"#\n",
	$rs, "\n",
	"# Number of particles and total number of shells\n",
	"# (nOccupied, nShells)\n",
	"#\n",
	$nOccupied, ", ", $nShells, "\n",
	"# Single-particle basis\n",
	"#\n",
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
	"# alpha:       t = alpha*t + (1-alpha)*t_old\n",
	"#              alpha in [0.0, 1.0] \n",
	$alpha, "\n",
	);
    
    #     Write the lines to file
    my $line = 0; 
    while($line < scalar @lines) {
	print $fileHandle $lines[$line];
	$line++;
    }
    
    #     Run a simulation with the produced input file
    my $statusCp = system("cp ".$inFile." input.dat");
    my $outFile = "outlist_nOcc".$nOccupied."_nShells".$nShells;
    my $statusCcd = system("mpirun -np 1 ./ccd.exe input.dat > ".$outFile." < /dev/null ");
    
    my $status = system("cat nshells_ene_old nshells_ene > ".$NEneTot);
    $status = system("cp ".$NEneTot." nshells_ene_old");
    
    close $fileHandle;
    $count++; 
    
}