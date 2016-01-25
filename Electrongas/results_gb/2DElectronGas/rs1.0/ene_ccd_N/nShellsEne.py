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
fileNameSRG26 = "N26_srg";
fileNameSRG42 = "N42_srg";
fileNameSRG58 = "N58_srg";
fileNameSRG74 = "N74_srg";
fileNameFCIQ26 = "N26_fciqmc";
fileNameCCD26 = "N26_ccd";
fileNameCCD42 = "N42_ccd";
fileNameCCD74 = "N74_ccd";
fileNameCCD114 = "N114_ccd";
fileNameCCD122 = "N122_ccd";
fileNameCCD138 = "N138_ccd";
dataSRG26 = loadtxt(fileNameSRG26)
dataSRG42 = loadtxt(fileNameSRG42)
dataSRG58 = loadtxt(fileNameSRG58)
dataSRG74 = loadtxt(fileNameSRG74)
dataFCIQ26 = loadtxt(fileNameFCIQ26)
dataCCD26 = loadtxt(fileNameCCD26)
dataCCD42 = loadtxt(fileNameCCD42)
dataCCD74 = loadtxt(fileNameCCD74)
dataCCD114 = loadtxt(fileNameCCD114)
dataCCD122 = loadtxt(fileNameCCD122)
dataCCD138 = loadtxt(fileNameCCD138)
nShellsSRG26 = dataSRG26[:,0]
nShellsSRG42 = dataSRG42[:,0]
nShellsSRG58 = dataSRG58[:,0]
nShellsSRG74 = dataSRG74[:,0]
nShellsFCIQ26 = dataFCIQ26[:,0]
nShellsCCD26 = dataCCD26[:,0]
nShellsCCD42 = dataCCD42[:,0]
nShellsCCD74 = dataCCD74[:,0]
nShellsCCD114 = dataCCD114[:,0]
nShellsCCD122 = dataCCD122[:,0]
nShellsCCD138 = dataCCD138[:,0]
eneSRG26 = dataSRG26[:,2]
eneSRG42 = dataSRG42[:,2]
eneSRG58 = dataSRG58[:,2]
eneSRG74 = dataSRG74[:,2]
eneFCIQ26 = dataFCIQ26[:,1]
eneCCD26 = dataCCD26[:,3] 
eneCCD42 = dataCCD42[:,3]
eneCCD74 = dataCCD74[:,3]
eneCCD114 = dataCCD114[:,3]
eneCCD122 = dataCCD122[:,3]
eneCCD138 = dataCCD138[:,3]

dashes1 = [5,2,2,2];
dashes2 = [5,2,1,2,1,2];

spring_green = (0, 1.0, 0)

#     Plot
p1 = plot(nShellsSRG26, eneSRG26, color="red", linewidth=1.2, 
          linestyle=":", label="SRG, N = 26")
p3 = plot(nShellsFCIQ26, eneFCIQ26, color=spring_green, 
          linestyle="noline", marker=".", label="FCIQMC, N = 26")
p4 = plot(nShellsCCD26, eneCCD26, color="blue", linewidth=1.2,
          linestyle=":", label="CCD, N = 26")
p2 = plot(nShellsSRG42, eneSRG42, color="red", linewidth=1.2,
          linestyle="--", label="SRG, N = 42")
p5 = plot(nShellsCCD42, eneCCD42, color="blue", linewidth=1.2,
          linestyle="--", label="CCD, N = 42")
p10 = plot(nShellsSRG58, eneSRG58, color="red", linewidth=1.2,
           linestyle="-.", label="SRG, N = 58")
p6 = plot(nShellsCCD74, eneCCD74, color="blue", linewidth=1.2,
          linestyle="-", label="CCD, N = 74")
p11 = plot(nShellsSRG74, eneSRG74, color="red", linewidth=1.2,
           linestyle="-", label="SRG, N=74")
p7 = plot(nShellsCCD114, eneCCD114, color="k", linewidth=1.2,
          linestyle=":", label="CCD, N = 114")
p8 = plot(nShellsCCD122, eneCCD122, color="k", linewidth=1.2,
          linestyle="-.", label="CCD, N = 122")
p9 = plot(nShellsCCD138, eneCCD138, color="k", linewidth=1.2,
          linestyle="-", label="CCD, N = 138")
# p1 = plot(nShellsSRG26, eneSRG26, color="red", linewidth=1.2, 
#           linestyle=":", marker=".", label="SRG, N = 26")
# p2 = plot(nShellsSRG42, eneSRG42, color="red", linewidth=1.2,
#           linestyle="-.", marker="o", markerfacecolor="None",
#           label="SRG, N = 42")
# p3 = plot(nShellsFCIQ26, eneFCIQ26, color=spring_green, 
#           linestyle="noline", marker=".", label="FCIQMC, N = 26")
# p4 = plot(nShellsCCD26, eneCCD26, color="blue", linewidth=1.2,
#           linestyle=":", marker=".", label="CCD, N = 26")
# p5 = plot(nShellsCCD42, eneCCD42, color="blue", linewidth=1.2,
#           linestyle="-.", marker="o", markerfacecolor="None",
#           label="CCD, N = 42")
# p6 = plot(nShellsCCD74, eneCCD74, color="blue", linewidth=1.2,
#           linestyle="-", marker="s", markerfacecolor="None",
#           label="CCD, N = 74")


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
#xlim(nShellsCCD.min(), nShellsCCD.max())
xlim(0, nShellsCCD26.max())
#ylim(ene_hf.min(), ene_hf.max())
ylim(-0.25, 0.3)

#     Set axis labels
xlabel('Number of sp energy shells', fontsize=16)
ylabel('Correlation energy per particle (Ry)', fontsize=16)
title("2D electron gas, $r_{s} = 1.0$")

fig.tight_layout()
frame = axis.get_frame()
frame.set_linewidth(0.1)

#show()

fig.savefig(fileName+'.pdf', format='pdf')
fig.savefig(fileName+'.eps', format='eps')
