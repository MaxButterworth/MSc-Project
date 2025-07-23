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
%                            i       ^            i    ^3
%   c(t+tau) = c(t-tau) - 2 ---- tau H c(t) + -------- H  c(t) - ...
%                           hbar     = -      3 hbar^3 =  -
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
function psi_new = tdse_dia (obj, step_size, pos, psi_old, psi_oldold)
global hamilt
psi_new = cell ( size(psi_old) );

n_traj = length(pos{1});

pot_mats = ham.pot_dia(pos);

H       = zeros ( n_traj , hamilt.coupling.n_eqs , hamilt.coupling.n_eqs );
H_trans = zeros ( size(H) ); % H^T
for m=1:hamilt.coupling.n_eqs
    for n=1:hamilt.coupling.n_eqs
        H(:,m,n)        = -1i * step_size * pot_mats(m,n,:);
        H_trans(:,n,m)  = H(:,m,n);
    end
end

% Has to be odd
maxdegree = 11;

if(mod(maxdegree,2)==0)
    prt.error ( 'maxdegree has to be odd' )
end

H_power             = cell(maxdegree,1);    % the j-th cell entry will be H^j
H_power_trans       = cell(size(H_power));  % the j-th cell entry will be (H^j)^T
H_power{1}          = H;
H_power_trans{1}    = H_trans;

% Computation of H^2
if(maxdegree > 2)
    H_power{2}          = zeros ( size(H) );
    H_power_trans{2}    = zeros ( size(H) );
    for m=1:hamilt.coupling.n_eqs
        for n=1:hamilt.coupling.n_eqs
            % Computation of H^2 and (H^2)^T
            % A * B = element-wise product of A^T and B with
            % the summmation dimension over the right dimension
            H_power{2}(:,m,n)       = sum( H_trans(:,:,m) .* H(:,:,n) , 2 );
            H_power_trans{2}(:,n,m) = H_power{2}(:,m,n);
        end
    end
end

% Computation of H^j for odd j
for j = 3:2:maxdegree
    H_power{j}          = zeros ( size(H) );
    H_power_trans{j}    = zeros ( size(H) );
    for m=1:hamilt.coupling.n_eqs
        for n=1:hamilt.coupling.n_eqs
            H_power{j}(:,m,n)       = sum( H_power_trans{j-2}(:,:,m) .* H_power{2}(:,:,n) , 2 );
            H_power_trans{j}(:,n,m) = H_power{j}(:,m,n);
        end
    end
end

% Computation of the exponential series
for m=1:hamilt.coupling.n_eqs
    psi_new{m} = psi_oldold{m};
    for n=1:hamilt.coupling.n_eqs
        for j = 1:2:maxdegree
            psi_new{m} = psi_new{m} + 2 / factorial(j) * H_power{j}(:,m,n) .* psi_old{n};
        end
    end
end

end
