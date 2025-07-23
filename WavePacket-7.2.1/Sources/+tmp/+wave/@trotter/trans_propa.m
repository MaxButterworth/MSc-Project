%----------------------------------------------------------------
% Perform propagation associated with transition dipoles (matrix)
%----------------------------------------------------------------
function trans_propa ( ~, psi, efield)
global hamilt

if isfield(hamilt, 'dip')
    for p = 1:length(hamilt.dip)
        if ~isempty (hamilt.dip{p})
            if hamilt.dip_trans(p)
                if abs(efield(p))>0

                    % Transform to adiabatic (transformation matrix D+)
                    for m = 1:hamilt.coupling.n_eqs
                        psi.new{m} = zeros(size(psi.dvr{1}));
                        for n=1:hamilt.coupling.n_eqs
                            psi.new{m} = psi.new{m} ...
                                + conj(hamilt.dip_eig_vecs{p}{n,m}) .* psi.dvr{n};
                        end
                        
                        % Propagate in adiabatic representation
                        psi.new{m} = psi.new{m} .* exp( - efield(p) * hamilt.dip_eig_vals{p}{m});
                    end
    
                    % Transform back to diabatic (transformation matrix D)
                    for m = 1:hamilt.coupling.n_eqs
                        psi.dvr{m} = zeros(size(psi.new{1}));
                        for n=1:hamilt.coupling.n_eqs
                            psi.dvr{m} = psi.dvr{m} ...
                                + hamilt.dip_eig_vecs{p}{m,n} .* psi.new{n};
                        end
                    end
    
                end
            end
        end
    end
end

end
            
   