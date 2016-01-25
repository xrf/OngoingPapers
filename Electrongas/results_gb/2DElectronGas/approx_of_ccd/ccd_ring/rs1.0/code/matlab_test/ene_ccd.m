function corr_energy = ene_ccd(t2matrix)

%
% Calculate the CCD correlation energy
%

v0125 = 0.313328534432575;
v0134 = -0.313328534432575;

corr_energy = (v0125*t2matrix(1,1) + v0134*t2matrix(1,2));