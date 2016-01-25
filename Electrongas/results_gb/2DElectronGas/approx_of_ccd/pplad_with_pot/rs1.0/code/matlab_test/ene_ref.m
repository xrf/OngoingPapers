function ene = ene_ref_2d(rs)

% Energies in units of Rydbergs (e^2/(2*a_{0}))
ene_kin = 1.0/(rs^2)
ene_pot = -8.0*sqrt(2)/(3.0*pi*rs)
ene = ene_kin + ene_pot;
