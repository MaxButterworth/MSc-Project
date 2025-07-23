%--------------------------------------------------------------------------
%
% Propagate quantum state vectors attached to trajectories
% The time evolution is given by the Schroedinger equation
% 
%    d           i   ^         
%   -- c(t) = - ---- H(t) c(t) 
%   dt -        hbar =    -    
% 
% for a time-step TAU using the exponential series
% 
%                      i       ^            1    ^2
%   c(t+tau) = c(t) - ---- tau H c(t) - -------- H  c(t) - ...
%                     hbar     = -      2 hbar^2 =  -
% 
%  
%
%--------------------------------------------------------------------------

% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2017-20.. Burkhard Schmidt's group
%
% see the README file for license details.

function psi_new = tdse_adi_trapezoidal (obj, step_size, mom, pot_mats, frc_mats, D, U, mom_old, pot_mats_old, frc_mats_old, D_old, U_old, psi_old, psi_oldold)
global hamilt space
psi_new = cell ( size(psi_old) );

n_traj = length(psi_old{1});

nac_all     = ham.nac_all (pot_mats,frc_mats,U,D);
nac_all_old = ham.nac_all (pot_mats_old,frc_mats_old,U_old,D_old);

H       = zeros ( n_traj , hamilt.coupling.n_eqs , hamilt.coupling.n_eqs );
H_trans = zeros ( size(H) ); % H^T
for m=1:hamilt.coupling.n_eqs
    for n=1:hamilt.coupling.n_eqs
        if m==n
            H_mn = 0.5 * ( D(m,:) + D_old(m,:))';
        else
            coupling = zeros(n_traj,1);
            for d = 1:space.n_dim
                coupling = coupling - 1i * 0.5 * (nac_all{m,n,d} .* mom{d} + nac_all_old{m,n,d} .* mom_old{d}) / space.dof{d}.mass;
            end
            H_mn = coupling;
        end
        H(:,m,n)        = -1i * step_size * H_mn;
        H_trans(:,n,m)  = H(:,m,n);
    end
end

% Can be odd or even
maxdegree = 3;%11;

H_power             = cell(maxdegree,1);    % the j-th cell entry will be H^j
H_power_trans       = cell(size(H_power));  % the j-th cell entry will be (H^j)^T
H_power{1}          = H;
H_power_trans{1}    = H_trans;


% Computation of H^j for j>1
for j = 2:maxdegree
    H_power{j}          = zeros ( size(H) );
    H_power_trans{j}    = zeros ( size(H) );
    for m=1:hamilt.coupling.n_eqs
        for n=1:hamilt.coupling.n_eqs
            H_power{j}(:,m,n)       = sum( H_power_trans{j-1}(:,:,m) .* H_power{1}(:,:,n) , 2 );
            H_power_trans{j}(:,n,m) = H_power{j}(:,m,n);
        end
    end
end

% Computation of the exponential series
for m=1:hamilt.coupling.n_eqs
    psi_new{m} = psi_old{m};
    for n=1:hamilt.coupling.n_eqs
        for j = 1:maxdegree
            psi_new{m} = psi_new{m} + 1 / factorial(j) * H_power{j}(:,m,n) .* psi_old{n};
        end
    end
end

end
