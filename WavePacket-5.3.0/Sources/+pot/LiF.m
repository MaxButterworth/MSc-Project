% This file is written to be used with the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2012-2016 Adam Kirrander
%
% see the README file for license details.
%--------------------------------------------------------------------------------------------
%
% LiF: Lithium Flouride
% Diabatic 2-state potential matrix for LiF 1Sigma ground and first
% excited states in atomic units. Ab initio data from Zeiri,
% Balint-Kurti JMS 99, 1-24 (1983).

% VERSION HISTORY:
% 2016 Dec: Current version adapted from pot_LiF_1S_Sadeghpour2.m for WavePacket program package.
% 2012 Apr: pot_LiF_1S_Sadeghpour2.m based on pot_rb2_1Sg_2Coulomb.m .

% INPUT:
% R = space.dvr.grid_ND{1}          Grid of internuclear distances R
% L = hamilt.pot.params.L           Rotational angular momentum L
% Nc = hamilt.coupling.n_eqs        Number of chnls
% Np = length(space.dvr.grid_ND{1}) Number radial grid points
% mred = space.dof{1}.mass;         Reduced mass LiF
% DEBUG = hamilt.pot.params.debug   DEBUG-flag (=1 for debugging output)  

% OUTPUT:
% DIABATIC POTENTIAL MATRIX WRITTEN TO GLOBAL PARAMETER hamilt.pot.grid_ND{i,j}  (1<ij<2)

% OUTPUT (internal):
%   Vd(Nc,Nc,Np)  diabatic potential matrix
%   VV(Nc,Np)     adiabatic energies (~diag of W(Nc,Nc) matrix)
%   B(Nc,Nc,Np)   non-adiabatic coupling (~off-diag of W(Nc,Nc) matrix), B-matrix
%   AA(Nc,Nc,Np)  non-adiabatic coupling (~off-diag of W(Nc,Nc) matrix), A-matrix
%   Elim(Nc)      asymptotic dissociation energy in each chnl
%   Vang(Np)      angular momentum component in adiabatic potential
%--------------------------------------------------------------------------------------------




function LiF
global hamilt space

util.disp (' ')
util.disp ('******************************************************************************************')
util.disp ('LiF 1Sigma PECs (ground and 1st excited): Zeiri, Balint-Kurti JMS 99, 1-24 (1983)')
util.disp ('')
util.disp ('Potentials given in atomic units and rotational angular momentum term added automatically.')
util.disp ('******************************************************************************************')
util.disp ( [ 'Rotational Angular Momentum L : ' num2str(hamilt.pot.params.L) ] )
util.disp ( [ 'Reduced mass (au)             : ' num2str(space.dof{1}.mass) ] )
util.disp ( [ 'DEBUG-flag                    : ' num2str(hamilt.pot.params.debug) ] )

% Check validity
if space.size.n_dim ~= 1
    util.error ('This potential is only for 1 dimension')
end

if hamilt.coupling.n_eqs ~= 2
    util.error ('This potential is only for 2 states')
end

if space.dvr.grid_ND{1}(1) < 0
    util.error ('Potential only defined for R>0')
end

if space.dof{1}.mass ~= 9265.d0 % not necessarily an error, but worth checking
    util.error ('Reduced mass of LiF is 9265 a.u. according to Sadeghpour Phys. Rev. Lett. (1999)') 
end


% Rename parameters for convenience
Nc      = hamilt.coupling.n_eqs;        % number of chnls
Np      = length(space.dvr.grid_ND{1}); % number radial grid points
R       = space.dvr.grid_ND{1};         % R grid
L       = hamilt.pot.params.L;          % rotational ang mom L
mred    = space.dof{1}.mass;            % reduced mass LiF
DEBUG   = hamilt.pot.params.debug;      % DEBUG-flag (=1 for debugging output)  


% GET POTENTIAL
% NB. rotational angular momentum term is added to adiabatic and diabatic potentials
    
% LOAD DATA  (READ DIABATIC PECs FROM SADEGHPOUR)
% xp         = load('LIF_pot_Sadeghpour.dat');
xp         =  load ('LIF_PECS.dat');
rin        = xp(:,1);   % 1st column has distances
[nrr, ncc] = size(xp);
xp         = xp(1:nrr,2:ncc);
ncc        = ncc - 1;  
if ncc ~= Nc+1 
    disp('unexpected number of columns in input')
    keyboard; 
end
  
% energy asymptotic limit
Elim(1) = xp(nrr,1);
Elim(2) = xp(nrr,2);   
Elim(2) = 0.d0; % adjust ion-pair limit, Sadeghpour et al set it to zero

  
% CALC ADIABATIC PECS
Ci  = zeros(Nc,Nc,nrr); % rotation matrix
Vad = zeros(Nc,nrr);    % adiabatic energies
Vdi = zeros(Nc,Nc,nrr); % diabatic energies
for i = 1:nrr
    Vdi(:,:,i)      = [xp(i,1) xp(i,3); xp(i,3) xp(i,2)];
    [eve,eva]       = eig(Vdi(:,:,i));
    Vad(1:Nc,i)     = diag(eva);
    Ci(1:Nc,1:Nc,i) = eve;
end



% align eigenvectors for maximum overlap (nb: arbitrary phase overall)
if Ci(1,1,nrr) < 0 Ci(1:Nc,1,nrr) = -Ci(1:Nc,1,nrr); end
if Ci(2,2,nrr) < 0 Ci(1:Nc,2,nrr) = -Ci(1:Nc,2,nrr); end  
for i = nrr-1:-1:1
    O = Ci(:,:,i).' * Ci(:,:,i+1);
    if O(1,1) < 0  Ci(1:Nc,1,i) = -Ci(1:Nc,1,i); end
    if O(2,2) < 0  Ci(1:Nc,2,i) = -Ci(1:Nc,2,i); end    
end
% want rotation matrix to be in the form
% [cos(th) -sin(th); sin(th) cos(th)]
% (not [cos(th) sin(th); -sin(th) cos(th)])
for i = 1:nrr
    Ci(:,:,i) = Ci(:,:,i).';
end


% OPTIONAL DEBUG/CHECK OF ADIABATIC TO DIABATIC TRANSFORMATION
if (DEBUG)
    error(1:nrr) = 0.d0;
    Vdt = zeros(Nc,Nc,nrr); % test
    for i = 1:nrr
        Vdt(:,:,i) = Ci(:,:,i).' * diag(Vad(:,i)) * Ci(:,:,i); % adiab to diab transform
        error(i) = abs(sum(sum(Vdt(:,:,i)-Vdi(:,:,i))));
    end
    clf
    plot(rin,squeeze(Vdi(1,1,:)),'k',rin,squeeze(Vdi(2,2,:)),'r',rin,squeeze(Vdi(1,2,:)),'b',rin,error,'g');
    legend('V(1,1)','V(2,2)','V(1,2)','Error')
    hold on
    plot(rin,squeeze(Vdt(1,1,:)),'--r',rin,squeeze(Vdt(2,2,:)),'--k',rin,squeeze(Vdt(1,2,:)),'--r');  
    xlabel('R (au)')
    ylabel('Diabatic matrix')
    title(['Adiab-diab transformation, Cumulative error=' num2str(sum(error))])
    hgsave('LiF_debug_adiab-diab-transf.fig');
end




%------------------------------------------------  
% MORE ABOUT COUPLINGS:
% calculate derivative of rotation matrix dC/dR
dCi   = zeros(Nc,Nc,nrr);
drin = diff(rin);
for i = 1:nrr-1
    dCi(:,:,i) = (Ci(:,:,i+1)-Ci(:,:,i))./drin(i);
end

% calculate first derivative coupling matrix, B
Bi = zeros(Nc,Nc,nrr);
for i = 1:nrr
    Bi(:,:,i)  = -dCi(:,:,i)*Ci(:,:,i).';
    Bi2(:,:,i) = -dCi(:,:,i)*Ci(:,:,i);
end

% we know that Bi should be antisymmetric
% i.e. that diagonal should be zero
Bi(1,1,:) = 0.d0;
Bi(2,2,:) = 0.d0;

% calculate derivative of Bi
dBi = zeros(Nc,Nc,nrr);
for i = 1:nrr-1
    dBi(:,:,i) = (Bi(:,:,i+1)-Bi(:,:,i))./drin(i);
end
  
% calculate A-matrix, <i|d2/dR2|j>
% A = B*B + dB/dR
Ai = zeros(Nc,Nc,nrr);
for i = 1:nrr
    Ai(1:Nc,1:Nc,i) = Bi(1:Nc,1:Nc,i) * Bi(1:Nc,1:Nc,i); % B^2
end
Ai = Ai + dBi;
% end couplings
%------------------------------------------------    

  
%------------------------------------------------  
% SPLINE OUTPUTS ONTO REQUIRED GRID
VV = zeros(Nc,Np);      % adiabatic energies
Vd = zeros(Nc,Nc,Np);   % diabatic matrix  
AA = zeros(Nc,Nc,Np);   % symmetric, second derivative couplings
B  = zeros(Nc,Nc,Np);   % antisymmetric B matrix, 1st derivative cplgs
C  = zeros(Nc,Nc,Np);   % transformation adiab-diab, Vd=C.'*Vad*C

% Check that requested distances R fall within defined region rin():
if R(1)<rin(1) || R(Np)>rin(nrr)
    util.error ([num2str(R(1)) '-' num2str(R(Np))  ' outside range ' num2str(rin(1)) '-' num2str(rin(nrr))]) 
end

% spline couplings
for i = 1:Nc    
    for j = 1:Nc
        AA(i,j,1:Np) = spline(rin, Ai(i,j,1:nrr), R(1:Np));
        B(i,j,1:Np)  = spline(rin, Bi(i,j,1:nrr), R(1:Np));
        C(i,j,1:Np)  = spline(rin, Ci(i,j,1:nrr), R(1:Np));
    end
end

% spline adiabatic potential curves
for i = 1:Nc    VV(i,1:Np) = spline(rin, Vad(i,:), R(1:Np));  end

% add rotational angular momentum term to adiabatic potential
Vang(1:Np) = (0.5d0*L.*(L+1.d0)/mred) ./ (R.^2);
for i = 1:Nc    
    VV(i,1:Np) = VV(i,1:Np) + Vang;  
end  

% diabatic potential by transformation from adiabatic
for i = 1:Np
    Vd(:,:,i) = C(:,:,i).' * diag(VV(:,i)) * C(:,:,i); % adiab to diab transform
end
%------------------------------------------------  


% OPTIONAL DEBUG/CHECK
if (DEBUG) % check spline
    clf
    plot(R,squeeze(Vd(1,1,:)),'k',R,squeeze(Vd(2,2,:)),'r',R,squeeze(Vd(1,2,:)),'b');
    hold on
    plot(rin,squeeze(Vdi(1,1,:)),'xr',rin,squeeze(Vdi(2,2,:)),'xk',rin,squeeze(Vdi(1,2,:)),'xb');
    legend('V(1,1)','V(2,2)','V(1,2)')
    xlabel('R (au)')
    ylabel('Diabatic matrix')
    title('Spline')
    hgsave('LiF_debug_spline.fig');  
end





% PLACE OUTPUT IN GLOBAL VARIABLES
hamilt.pot.grid_ND{1,1} = squeeze(Vd(1,1,1:Np));
hamilt.pot.grid_ND{2,2} = squeeze(Vd(2,2,1:Np));
hamilt.pot.grid_ND{1,2} = squeeze(Vd(1,2,1:Np));

%keyboard

% OPTIONAL DEBUG/CHECK
if (DEBUG)
    % SAVE DATA
    save LiF_debug_data.mat R Np Nc mred L VV Vang B AA Elim Vd C Vad Vdi rin
    % PLOT DATA
    cc = {'k','r','b','g','y','m','c','k--','r--','b--','g--','y--','m--','c--'};
    %[ion,Ediss] = calc_ionic(R); % calc ionic reference  
    clf
    for i = 1:Nc
        subplot(2,1,1), plot(R,VV(i,:),cc{i}); hold on
    end
    %subplot(2,1,1),plot(R,ion,'--r'); 
    for i = 1:Nc-1
        for j = i+1:Nc
            subplot(2,1,2), plot(R,squeeze(B(i,j,:)),'k'); hold on
        end
    end
    subplot(2,1,1), axis([min(R) max(R) min(min(VV)) max(max(VV))]); grid
    %subplot(2,1,2), axis([min(R) max(R) min(min(min(B))) max(max(max(B)))]); grid; 
    subplot(2,1,2), axis([min(R) max(R) min(min(min(B))) 0]); grid; 
    subplot(2,1,1), ylabel('PECs LiF2 1Sg')
    subplot(2,1,2), xlabel('R (au)'); ylabel('couplings');
    hgsave('LiF_debug_adiab+coupl.fig')
end

% End Function LiF
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



