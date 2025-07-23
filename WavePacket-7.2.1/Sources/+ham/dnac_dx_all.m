%--------------------------------------------------------------------------
%
% Compute derivative of first order non-adiabatic coupling vector between 
% all adiatic states from derivative of diabatic force matrices, 
% eigenvectors, eigenvalues, force matrices, nac-vectors.
% 
%--------------------------------------------------------------------------

% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2018-.... Burkhard Schmidt and Leonardo Cancissu Araujo
%
% see the README file for license details.

function dnac_dx_all = dnac_dx_all (dfrc_dx_mats,U,D,frc_adi,nac_all)
global space hamilt

n_traj = length(U(1,1,:));

if n_traj>0
    
    dnac_dx_all = cell(hamilt.coupling.n_eqs,hamilt.coupling.n_eqs,space.n_dim,space.n_dim);
    
    switch hamilt.coupling.n_eqs
        case 1
            for m=1:hamilt.coupling.n_eqs
                for n=1:hamilt.coupling.n_eqs
                    for k=1:space.n_dim
                        for i=1:space.n_dim
                            dnac_dx_all{m,n,k,i} = zeros(n_traj,1);
                        end
                    end
                end
            end
        
        otherwise
            
            % Calculate non-adiabatic couplings, using also the diabatic(!) forces
            for m=1:hamilt.coupling.n_eqs
                for k=1:space.n_dim
                    for i=1:space.n_dim
                        dnac_dx_all{m,m,k,i} = zeros(n_traj,1);
                    end
                end
                
                for n=(m+1):hamilt.coupling.n_eqs
                    
                    dnac_dx_mn = ham.dnac_dx_mn (dfrc_dx_mats,U,D,frc_adi,nac_all,m,n);
                    
                    for k=1:space.n_dim
                        for i=1:space.n_dim
                            dnac_dx_all{m,n,k,i} =   dnac_dx_mn{k,i};
                            dnac_dx_all{n,m,k,i} = - dnac_dx_mn{k,i};
                        end
                    end
                end
            end
    end
end
