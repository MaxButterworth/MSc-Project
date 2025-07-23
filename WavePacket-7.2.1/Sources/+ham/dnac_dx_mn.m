%--------------------------------------------------------------------------
%
% Compute derivative of first order non-adiabatic coupling vector between 
% adiatic states l and m from derivative of diabatic force matrices, 
% eigenvectors, eigenvalues, force matrices, nac-vectors.
% 
%--------------------------------------------------------------------------

% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2018-.... Burkhard Schmidt and Leonardo Cancissu Araujo
%
% see the README file for license details.

function dnac_dx_lm = dnac_dx_mn (dfrc_dx_mats,U,D,frc_all,nac_all,l,m)
global space hamilt

n_traj = length(U(1,1,:));

if n_traj>0
    
    dnac_dx_lm = cell(space.n_dim);
    
    switch hamilt.coupling.n_eqs
        case 1
            for k=1:space.n_dim
                for i=1:space.n_dim
                    dnac_dx_lm{k,i} = zeros(n_traj,1);
                end
            end
            
        otherwise
            
            % Calculate non-adiabatic couplings, using also the diabatic(!) forces
            for k=1:space.n_dim
                for i=1:space.n_dim
                    
                    sum_nac_F_diff =   nac_all{l,m,k} .* ( frc_all{l,i} - frc_all{m,i} ) ...
                                     + nac_all{l,m,i} .* ( frc_all{l,k} - frc_all{m,k} );
                    
                    % Computation of u_l^T * G_ki * u_m
                    G_ki = dfrc_dx_mats{k,i}(:,:,:);
                    u_l_G_ki = reshape( sum( U(:,l,:) .* G_ki , 1 ) , [hamilt.coupling.n_eqs,n_traj] );
                    u_m = squeeze(U(:,m,:));
                    u_l_G_ki_u_n = sum( u_l_G_ki .* u_m ,1)';
                    
                    sum_nac_E_diff = zeros(n_traj,1);
                    for q = 1:hamilt.coupling.n_eqs 
                        if(q ~= l && q~= m)
                            sum_nac_E_diff =   nac_all{q,l,k} .* nac_all{q,m,i} .* ( D(q,:)' - D(m,:)' ) ...
                                             + nac_all{q,m,k} .* nac_all{q,l,i} .* ( D(q,:)' - D(l,:)' ) ...
                                             + sum_nac_E_diff;
                        end
                    end
                    
                    dnac_dx_lm{k,i} = ( sum_nac_F_diff + u_l_G_ki_u_n + sum_nac_E_diff ) ./ (D(l,:)'-D(m,:)');
                end
            end
    end
end
