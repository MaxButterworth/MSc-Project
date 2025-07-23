%--------------------------------------------------------------------------
%
% Compute force along m-th potential energy surface 
% from diabatic force matrices and eigenvector matrices
%
%--------------------------------------------------------------------------

% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2018-.... Burkhard Schmidt and Leonardo Cancissu Araujo
%
% see the README file for license details.

function frc_adi_all = frc_adi_all (pot_mats,frc_mats,U)
global space hamilt

n_traj = length(U(1,1,:));

if n_traj>0
    
    frc_adi_all = cell(hamilt.coupling.n_eqs,space.n_dim);
    
    switch hamilt.coupling.n_eqs
        case 1
            for d=1:space.n_dim
                frc_adi_all{1,d} = frc_mats{d};
            end
        
        otherwise
            
            for n = 1:hamilt.coupling.n_eqs
                frc_n = ham.frc_adi (pot_mats,frc_mats,U,n);
                for d=1:space.n_dim
                    frc_adi_all{n,d} = frc_n{d};
                end
            end
    end
end
        
