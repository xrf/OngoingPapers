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
fileName = "nOrbitsEneRs";
fileNameSRG1 = "srg_rs0.5";
fileNameSRG2 = "srg_rs1.0";
fileNameSRG3 = "srg_rs2.0";
fileNameSRG4 = "srg_rs5.0";
fileNameFCIQ1 = "fciqmc_rs0.5";
fileNameFCIQ2 = "fciqmc_rs1.0";
fileNameFCIQ3 = "fciqmc_rs2.0";
fileNameFCIQ4 = "fciqmc_rs5.0";
fileNameCCD1 = "ccd_rs0.5";
fileNameCCD2 = "ccd_rs1.0";
fileNameCCD3 = "ccd_rs2.0";
fileNameCCD4 = "ccd_rs5.0";

dataSRG1 = loadtxt(fileNameSRG1)
dataSRG2 = loadtxt(fileNameSRG2)
dataSRG3 = loadtxt(fileNameSRG3)
dataSRG4 = loadtxt(fileNameSRG4)
dataFCIQ1 = loadtxt(fileNameFCIQ1)
dataFCIQ2 = loadtxt(fileNameFCIQ2)
dataFCIQ3 = loadtxt(fileNameFCIQ3)
dataFCIQ4 = loadtxt(fileNameFCIQ4)
dataCCD1 = loadtxt(fileNameCCD1)
dataCCD2 = loadtxt(fileNameCCD2)
dataCCD3 = loadtxt(fileNameCCD3)
dataCCD4 = loadtxt(fileNameCCD4)

orbitsSRG1 = 1/dataSRG1[:,1]
orbitsSRG2 = 1/dataSRG2[:,1]
orbitsSRG3 = 1/dataSRG3[:,1]
orbitsSRG4 = 1/dataSRG4[:,1]
orbitsFCIQ1 = 1/dataFCIQ1[:,0]
orbitsFCIQ2 = 1/dataFCIQ2[:,0]
orbitsFCIQ3 = 1/dataFCIQ3[:,0]
orbitsFCIQ4 = 1/dataFCIQ4[:,0]
orbitsCCD1 = 1/dataCCD1[:,1]
orbitsCCD2 = 1/dataCCD2[:,1]
orbitsCCD3 = 1/dataCCD3[:,1]
orbitsCCD4 = 1/dataCCD4[:,1]

eneSRG1 = dataSRG1[:,3]*7.0
eneSRG2 = dataSRG2[:,3]*7.0
eneSRG3 = dataSRG3[:,3]*7.0
eneSRG4 = dataSRG4[:,3]*7.0
eneFCIQ1 = dataFCIQ1[:,1]
eneFCIQ2 = dataFCIQ2[:,1]
eneFCIQ3 = dataFCIQ3[:,1]
eneFCIQ4 = dataFCIQ4[:,1]
eneCCD1 = dataCCD1[:,3]*7.0 
eneCCD2 = dataCCD2[:,3]*7.0
eneCCD3 = dataCCD3[:,3]*7.0
eneCCD4 = dataCCD4[:,3]*7.0

dashes1 = [5,2,2,2];
dashes2 = [5,2,1,2,1,2];

spring_green = (0, 1.0, 0)

#     Plot
#p1 = plot(orbitsSRG, eneSRG, color="red", linewidth=1.2, 
#          linestyle=":", marker=".", label="SRG")
p9 = plot(orbitsSRG4, eneSRG4, color="red", linewidth=1.2,
          linestyle="-", marker=".", label="SRG, rs = 5.0")
p8 = plot(orbitsCCD4, eneCCD4, color="blue", linewidth=1.2,
          linestyle="-", marker=".", label="CCD, rs = 5.0")
p7 = plot(orbitsFCIQ4, eneFCIQ4, color=spring_green, linewidth=1.2,
          linestyle="-", marker=".", label="FCIQMC, rs = 5.0")

p10 = plot(orbitsSRG3, eneSRG3, color="red", linewidth=1.2,
           linestyle="--", marker=".", label="SRG, rs = 2.0")
p6 = plot(orbitsCCD3, eneCCD3, color="blue", linewidth=1.2,
          linestyle="--", marker=".", label="CCD, rs = 2.0")
p5 = plot(orbitsFCIQ3, eneFCIQ3, color=spring_green, linewidth=1.2,
          linestyle="--", marker=".", label="FCIQMC, rs = 2.0")

p11 = plot(orbitsSRG2, eneSRG2, color="red", linewidth=1.2,
           linestyle="-.", marker=".", label="SRG, rs = 1.0")
p4 = plot(orbitsCCD2, eneCCD2, color="blue", linewidth=1.2,
          linestyle="-.", marker=".", label="CCD, rs = 1.0")
p3 = plot(orbitsFCIQ2, eneFCIQ2, color=spring_green, linewidth=1.2,
          linestyle="-.", marker=".", label="FCIQMC, rs = 1.0")

p12 = plot(orbitsSRG1, eneSRG1, color="red", linewidth=1.2,
           linestyle=":", marker=".", label="SRG, rs = 0.5")
p2 = plot(orbitsCCD1, eneCCD1, color="blue", linewidth=1.2,
          linestyle=":", marker=".", label="CCD, rs = 0.5")
p1 = plot(orbitsFCIQ1, eneFCIQ1, color=spring_green, linewidth=1.2,
          linestyle=":", marker=".", label="FCIQMC, rs = 0.5")


#plot([0, 0.021], [-0.217, -0.217], 'k', lw=1.2, linestyle="--",
#    label="VMC, TDL")
# plot([0, N_ele.max()], [0.242, 0.242], 'k', lw=1.2, linestyle="-.")
# plot([0, N_ele.max()], [0.0945, 0.0945], 'k', lw=1.2, linestyle="-")


#     Legend
#leg = legend(loc='upper left', labelspacing=0.1)
#fr_leg = leg.get_frame()
#fr_leg.set_lw(0.6)

#     Set minor and major ticks
#majorLocatorX = MultipleLocator(50)
#minorLocatorX = MultipleLocator(10)
majorLocatorY = MultipleLocator(0.10)
minorLocatorY = MultipleLocator(0.02)
axis = fig.gca()
#axis.xaxis.set_major_locator(majorLocatorX)
#axis.xaxis.set_minor_locator(minorLocatorX)
axis.yaxis.set_major_locator(majorLocatorY)
axis.yaxis.set_minor_locator(minorLocatorY)

xticks([0.0, 1./778, 1./358, 1./114],
       [r'$\infty^{-1}$', r'$778^{-1}$', r'$358^{-1}$', r'$114^{-1}$'])
#yticks([-0.217], [-0.217])

#     Limits on axes
#xlim(nShellsCCD.min(), nShellsCCD.max())
#xlim(0.0, orbitsSRG.max())
xlim(0.0, 0.009)
#ylim(ene_hf.min(), ene_hf.max())
#ylim(-0.62, 0.6)

#     Set axis labels
xlabel('Orbits$^{-1}$', fontsize=16)
ylabel('Correlation energy per particle (Ry)', fontsize=16)
title("3D electron gas, $N = 14$")

fig.tight_layout()
frame = axis.get_frame()
frame.set_linewidth(0.1)

#show()

fig.savefig(fileName+'.pdf', format='pdf')
fig.savefig(fileName+'.eps', format='eps')
