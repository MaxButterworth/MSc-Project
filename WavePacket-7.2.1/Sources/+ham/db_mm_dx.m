%--------------------------------------------------------------------------
%
% Compute derivative of second order non-adiabatic coupling element at 
% adiabatic state m from derivative of diabatic force matrices, 
% eigenvectors, eigenvalues, force matrices, nac-vectors.
% 
%--------------------------------------------------------------------------

% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2018-.... Burkhard Schmidt and Leonardo Cancissu Araujo
%
% see the README file for license details.

function db_mm_dx = db_mm_dx (nac_all,dnac_dx_all,m)
global space hamilt

n_traj = length(nac_all{m,m,1});

if n_traj>0
    
    db_mm_dx = cell(space.n_dim,1);
    
    if (hamilt.coupling.n_eqs >= 2)
        
        for k=1:space.n_dim
            db_mm_dxk = zeros (n_traj,1);
            
            for l = 1:hamilt.coupling.n_eqs
                if( l ~= m )
                    for i = 1:space.n_dim
                        db_mm_dxk = db_mm_dxk - 2.* dnac_dx_all{l,m,k,i} .* nac_all{l,m,i};
                    end
                end
            end
            
            db_mm_dx{k} = db_mm_dxk / ( 2 * space.dof{1}.mass );
        end
    end
end
