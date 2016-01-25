function ene_ref = ene_ref_num_3d(rs)

nParticles = 14;
volume = nParticles*4.0*pi*rs^3/3.0; 
area = volume^(2.0/3.0);
length = volume^(1.0/3.0);

ene_kin = 2.0*24.0*pi^2/(area*nParticles)

ene_pot = -2.0*25.5/(pi*length*nParticles)

%-1/(pi*length)

ene_ref = ene_kin + ene_pot;
