%--------------------------------------------------------------------------
% 
% Compute adiabatic potentials, i.e. eigenvalues and 
% eigenvector matrix of diabatic potential matrix
% 
%--------------------------------------------------------------------------

% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2018-.... Burkhard Schmidt and Leonardo Cancissu Araujo
%
% see the README file for license details.

function u_sa = pot_eig_sa(mom, pot_mats, frc_mats, D, U)
global hamilt space

n_traj = length(pot_mats(1,1,:));
u_sa = zeros ( hamilt.coupling.n_eqs , hamilt.coupling.n_eqs , n_traj );

nac_norm_all = ham.nac_norm_all (pot_mats,frc_mats,U,D);

H = zeros ( hamilt.coupling.n_eqs , hamilt.coupling.n_eqs , n_traj);
for m=1:hamilt.coupling.n_eqs
    for n=1:hamilt.coupling.n_eqs
        if m==n
            H_mn = D(m,:)' .* ( D(2,:) - D(1,:) )';
        else
            coupling = zeros(n_traj,1);
            for d = 1:space.n_dim
                coupling = coupling - 1i * nac_norm_all{m,n,d} .* mom{d} / space.dof{d}.mass;
            end
            H_mn = coupling;
        end
        H(m,n,:)        = H_mn;
    end
end

switch hamilt.coupling.n_eqs
        
    case 2
        V_11(:,1) = H(1,1,:);
        V_22(:,1) = H(2,2,:);
        V_12(:,1) = 1i * H(1,2,:);
        
        % Adiabatic potential energy curves|surfaces
        % Analytic solutions for a 2-state problem
        % See e.g. Eq. (4.8) in doi:10.1063/1.1522712
        dlt = (V_11 - V_22)/2;
        eta = (V_11 + V_22)/2;
        rho = sqrt ( dlt.^2 + V_12.^2 );
        
        % Compute eigenvectors with explicit formulas (not for SSSH, variant 2)
        ind_n = V_12~=0;
        n21   = (V_11 > V_22)*1.0;
        n22   = (V_11 < V_22)*1.0;
        
        b21          = V_11(ind_n)-V_22(ind_n)+2*rho(ind_n);
        b22          = 2*V_12(ind_n);
        n22(ind_n)   = b22 ./ sqrt(b21.^2 + b22.^2);
        n21(ind_n)   = b21 ./ sqrt(b21.^2 + b22.^2);
        
        % Save previous and new eigenvectors
        u_sa(1,1,:) =  - n22;
        u_sa(2,1,:) =   1i * n21;
        u_sa(1,2,:) =   1i * n21;
        u_sa(2,2,:) =  - n22;
        
    otherwise
        
        % Diagonalize diabatic potential matrices for each trajectory
        for t = 1:n_traj
            % Entries of vector D are eigenvalues
            % Columns of matrix U are right eigenvectors
            [u_sa(:,:,t),D] = eig(H(:,:,t),'vector');
        end
end

% time.counter.diag.total = time.counter.diag.total + n_traj;

end

