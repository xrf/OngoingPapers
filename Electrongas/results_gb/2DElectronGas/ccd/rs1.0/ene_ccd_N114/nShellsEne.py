from numpy import *
from matplotlib.pyplot import *
from matplotlib import rc
from matplotlib.ticker import MultipleLocator, FormatStrFormatter

#     Use Latex font as default
rc('font',**{'family':'serif','serif':['Helvetica']})
rc('text', usetex=True)
matplotlib.rcParams.update({'font.size': 14})
rc('legend',**{'fontsize':14})   

#     Initialize a figure
fig = figure(figsize=(6,4), dpi=60)

#     Read data file
fileName = "nShellsEne";
#fileNameSRG = "N26_srg";
#fileNameFCIQ = "N26_fciqmc_data";
fileNameCCD = "N114_ccd";

#dataSRG = loadtxt(fileNameSRG)
#dataFCIQ = loadtxt(fileNameFCIQ)
dataCCD = loadtxt(fileNameCCD)
#nShellsSRG = dataSRG[:,0]
#nShellsFCIQ = dataFCIQ[:,0]
nShellsCCD = dataCCD[:,0]
#eneSRG = dataSRG[:,2]
#eneFCIQ = dataFCIQ[:,1]
eneCCD = dataCCD[:,3] 

dashes1 = [5,2,2,2];
dashes2 = [5,2,1,2,1,2];

spring_green = (0, 1.0, 0)

#     Plot
#p1 = plot(nShellsSRG, eneSRG, color="red", linewidth=1.2, 
#          linestyle="-", label="SRG")
#p2 = plot(nShellsFCIQ, eneFCIQ, color=spring_green, 
#          marker=".", linestyle="noline", label="FCIQMC")
p3 = plot(nShellsCCD, eneCCD, color="blue", linewidth=1.2,
          linestyle="-.", label="CCD")

# plot([0, N_ele.max()], [0.5525, 0.5525], 'k', lw=1.2, linestyle=":")
# plot([0, N_ele.max()], [0.242, 0.242], 'k', lw=1.2, linestyle="-.")
# plot([0, N_ele.max()], [0.0945, 0.0945], 'k', lw=1.2, linestyle="-")


#     Legend
leg = legend(loc='upper right', labelspacing=0.1)
fr_leg = leg.get_frame()
fr_leg.set_lw(0.6)

#     Set minor and major ticks
majorLocatorX = MultipleLocator(50)
minorLocatorX = MultipleLocator(10)
majorLocatorY = MultipleLocator(0.05)
minorLocatorY = MultipleLocator(0.01)
axis = fig.gca()
axis.xaxis.set_major_locator(majorLocatorX)
axis.xaxis.set_minor_locator(minorLocatorX)
axis.yaxis.set_major_locator(majorLocatorY)
axis.yaxis.set_minor_locator(minorLocatorY)

#     Limits on axes
xlim(nShellsCCD.min(), nShellsCCD.max())
#ylim(ene_hf.min(), ene_hf.max())
#ylim(0.0, 0.9)

#     Set axis labels
xlabel('Number of sp energy shells', fontsize=16)
ylabel('Correlation energy per particle (Ry)', fontsize=16)
title("2D electron gas, $N = 114$, $r_{s} = 1.0$")

fig.tight_layout()
frame = axis.get_frame()
frame.set_linewidth(0.1)

#show()

fig.savefig(fileName+'.pdf', format='pdf')
fig.savefig(fileName+'.eps', format='eps')
