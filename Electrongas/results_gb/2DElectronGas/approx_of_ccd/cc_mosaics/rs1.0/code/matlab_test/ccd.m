clear all;
format long;

max_iter = 100;
tolerance = 0.000000000001;

%
% Do the two first iteration steps of the CCD self-consistency
% loop.
% 

% Initialize CCD T2 matrix
t2mat_old = zeros(2,1);

iter = 1;
difference = 100.0;
energyOld = 100000.0;
while ((iter < max_iter) & (difference > tolerance))

    % Calculate T2 matrix
    t2mat = t2_matrix(t2mat_old);
    % Calculate CCD correlation energy
    energy = ene_ccd(t2mat)
    t2mat_old = t2mat;
    
    difference = abs(energy - energyOld);
    energyOld = energy;
    
    iter = iter + 1
end


