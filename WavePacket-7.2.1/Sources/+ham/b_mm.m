%--------------------------------------------------------------------------
%
% Compute second order non-adiabatic coupling element at adiabatic state m 
% from derivative of diabatic force matrices, 
% eigenvectors, eigenvalues, force matrices, nac-vectors.
% 
%--------------------------------------------------------------------------

% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2018-.... Burkhard Schmidt and Leonardo Cancissu Araujo
%
% see the README file for license details.

function b_mm = b_mm (nac_all,m)
global space hamilt

n_traj = length(nac_all{m,m,1});

if n_traj>0
    
    b_mm = zeros(n_traj,1);
    
    if (hamilt.coupling.n_eqs >= 2)
        for l = 1:hamilt.coupling.n_eqs
            if( l ~= m )
                for d = 1:space.n_dim
                    b_mm = b_mm - nac_all{l,m,d}.^2;
                end
            end
        end
    end
end
