%--------------------------------------------------------------------------
%
% Propagate quantum state vectors attached to trajectories
% The time evolution is given by the Schroedinger equation
% 
%    d           i   ^         
%   -- c(t) = - ---- H(t) c(t) 
%   dt -        hbar =    -    
% 
% for a time-step TAU using second (or higher) order differencing
% 
%                                                      +
%   c(t+tau) = U(t+tau) * exp(-i E(t+tau) tau) U(t+tau)  c(t)
% 
% This version for diabatic representation 
%
%--------------------------------------------------------------------------

% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2017-20.. Burkhard Schmidt's group
%
% see the README file for license details.
function psi_new = tdse_dia_2 (obj, step_size, pos, psi_old, psi_oldold)
global hamilt
psi_new = cell ( size(psi_old) );

n_traj = length(pos{1});

pot_mats = ham.pot_dia(pos);

[U_new , D] = ham.pot_eig_adi(pot_mats);

H = zeros ( n_traj , hamilt.coupling.n_eqs , hamilt.coupling.n_eqs );
U = zeros ( n_traj , hamilt.coupling.n_eqs , hamilt.coupling.n_eqs );
for m=1:hamilt.coupling.n_eqs
    for n=m:hamilt.coupling.n_eqs
        U(:,m,n) = U_new(m,n,:); % different order of the dimensions compared to obj.U_new
        U(:,n,m) = U_new(n,m,:);
    end
    H(:,m,m) = exp(-1i * step_size * D(m,:)');
end

% Computation of H*U^T*U_old
U_H_U_trans = zeros ( size(H) );
for m=1:hamilt.coupling.n_eqs
    for n=1:hamilt.coupling.n_eqs
        for k=1:hamilt.coupling.n_eqs
            U_H_U_trans(:,m,n) = U_H_U_trans(:,m,n) + U(:,m,k) .* H(:,k,k) .* U(:,n,k);
        end
    end
end

% Computation of the exponential matrix
for m=1:hamilt.coupling.n_eqs
    psi_new{m} = zeros ( size(psi_old{m}) );
    for n=1:hamilt.coupling.n_eqs
        psi_new{m} = psi_new{m} + U_H_U_trans(:,m,n) .* psi_old{n};
    end
end

end
