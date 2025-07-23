%--------------------------------------------------------------------------
%
% Compute all first order non-adiabatic coupling vector between 
% adiatic states m and n from diabatic potential matrices, 
% eigenvectors and eigenvalues thereof, and force matrices.
% 
%--------------------------------------------------------------------------

% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2018-.... Burkhard Schmidt and Leonardo Cancissu Araujo
%
% see the README file for license details.

function nac_all = nac_all (pot_mats,frc_mats,U,D)
global space hamilt

n_traj = length(U(1,1,:));

if n_traj>0
    
    nac_all = cell(hamilt.coupling.n_eqs,hamilt.coupling.n_eqs,space.n_dim);
    
    switch hamilt.coupling.n_eqs
        case 1
            for m=1:hamilt.coupling.n_eqs
                for n=1:hamilt.coupling.n_eqs
                    for d=1:space.n_dim
                        nac_all{m,n,d} = zeros(n_traj,1);
                    end
                end
            end
        
        otherwise
            
            % Calculate non-adiabatic couplings, using also the diabatic(!) forces
            for m=1:hamilt.coupling.n_eqs
                for d=1:space.n_dim
                    nac_all{m,m,d} = zeros(n_traj,1);
                end
                
                for n=(m+1):hamilt.coupling.n_eqs
                    
                    nac_mn = ham.nac_mn (pot_mats,frc_mats,U,D,m,n);
                    
                    for d=1:space.n_dim
                        nac_all{m,n,d} =   nac_mn{d};
                        nac_all{n,m,d} = - nac_mn{d};
                    end
                end
            end
    end
end
