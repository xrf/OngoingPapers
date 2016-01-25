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
fileName = "NEneRef3d_rs2.0";
data_ene = loadtxt(fileName)
N_ele = data_ene[:,0]
ene_kin = data_ene[:,1]
ene_pot = data_ene[:,2] + 0.7
ene_ref = data_ene[:,3] 

dashes1 = [5,2,2,2];
dashes2 = [5,2,1,2,1,2];

spring_green = (0, 1.0, 0)

#     Plot
p1= plot(N_ele, ene_kin, color=spring_green, linewidth=1.2, linestyle=":",
         label="Kinetic energy")
p2 = plot(N_ele, ene_pot, color="blue", linewidth=1.2, linestyle="-.",
          label="Potential energy + 0.7")
p3 = plot(N_ele, ene_ref, color="red", linewidth=1.2, linestyle="-",
          label="Reference energy")

plot([0, N_ele.max()], [0.5525, 0.5525], 'k', lw=1.2, linestyle=":")
plot([0, N_ele.max()], [0.242, 0.242], 'k', lw=1.2, linestyle="-.")
plot([0, N_ele.max()], [0.0945, 0.0945], 'k', lw=1.2, linestyle="-")


#     Legend
leg = legend(loc='upper right', labelspacing=0.1)
fr_leg = leg.get_frame()
fr_leg.set_lw(0.6)

#     Set minor and major ticks
majorLocatorX = MultipleLocator(1000)
minorLocatorX = MultipleLocator(100)
majorLocatorY = MultipleLocator(0.2)
minorLocatorY = MultipleLocator(0.05)
axis = fig.gca()
axis.xaxis.set_major_locator(majorLocatorX)
axis.xaxis.set_minor_locator(minorLocatorX)
axis.yaxis.set_major_locator(majorLocatorY)
axis.yaxis.set_minor_locator(minorLocatorY)

#     Limits on axes
xlim(0, N_ele.max())
#ylim(ene_hf.min(), ene_hf.max())
ylim(0.0, 0.9)

#     Set axis labels
xlabel('Number of particles', fontsize=16)
ylabel('Energy per particle (Ry)', fontsize=16)
title("The three-dimensional electron gas, $r_{s} = 2.0$")

fig.tight_layout()
frame = axis.get_frame()
frame.set_linewidth(0.1)

#show()

fig.savefig(fileName+'.pdf', format='pdf')
fig.savefig(fileName+'.eps', format='eps')
