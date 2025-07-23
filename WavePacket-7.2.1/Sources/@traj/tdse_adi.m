%--------------------------------------------------------------------------
%
% Propagate quantum state vectors attached to trajectories
% The time evolution is given by the Schroedinger equation
%
%    d           i   ^
%   -- c(t) = - ---- H(t) c(t)
%   dt -        hbar =    -
%
% for a time-step TAU by direct diagonalization of the Hamiltonian
%
%                     (    i          )  +
%   c(t+tau) = S  exp ( - ---- E tau  ) S  c(t)
%   -          =      (   hbar -      ) =  -
%
%            +
% where S H S = diag(E) defines the adiabatic representation
%
%--------------------------------------------------------------------------

% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2017-20.. Burkhard Schmidt's group
%
% see the README file for license details.

function psi_new = tdse_adi (obj, step_size, D, obj_U_new, obj_U_old, psi_old)
global hamilt

n_traj = length(D(1,:));

% Note that obj.U_old is always available, even at the very first time
% This is because @pot\eval_V_F.m has already been called from inside
% method "hamilton" (calculating mean energies) which is called during  
% the initialization of each of the propagators (method "traj_init").
U_new       = zeros(n_traj, hamilt.coupling.n_eqs, hamilt.coupling.n_eqs);
U_old       = zeros(n_traj, hamilt.coupling.n_eqs, hamilt.coupling.n_eqs);
H           = zeros(n_traj, hamilt.coupling.n_eqs, hamilt.coupling.n_eqs);
psi_new     = cell ( size(psi_old) );

for m=1:hamilt.coupling.n_eqs
    for n=1:hamilt.coupling.n_eqs
        U_new(:,m,n)        = obj_U_new(m,n,:); % different order of the dimensions compared to obj.U_new
        U_old(:,m,n)        = obj_U_old(m,n,:); % different order of the dimensions compared to obj.U_old
    end
    H(:,m,m) = exp(-1i * step_size * D(m,:)');
end

% Computation of H*U^T*U_old
H_U_trans_U_old = zeros ( size(H) );
for m=1:hamilt.coupling.n_eqs
    for n=1:hamilt.coupling.n_eqs
        H_U_trans_U_old(:,m,n) = H(:,m,m) .* sum( U_new(:,:,m) .* U_old(:,:,n) , 2 );
    end
end

% Computation of the exponential matrix
for m=1:hamilt.coupling.n_eqs
    psi_new{m} = zeros ( size(psi_old{m}) );
    for n=1:hamilt.coupling.n_eqs
        psi_new{m} = psi_new{m} + H_U_trans_U_old(:,m,n) .* psi_old{n};
    end
end

end
