function ene = ene_ref_3d(rs)

% Energies in units of Rydbergs (e^2/(2*a_{0}))
ene_kin = 2.21/(rs^2)
ene_pot = - 0.916/rs
ene = ene_kin + ene_pot;



