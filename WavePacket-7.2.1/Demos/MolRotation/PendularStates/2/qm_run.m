% Copyright (C) 2017-.... Burkhard Schmidt
% Copyright (C) 2008-2011 Ulf Lorenz
%
% This file is under public domain. You may modify or incorporate it into other
% works without any restrictions.


%---------------------------------------------------------------
% This script intends to reproduce the second graph in Figure 1
% of [J.Chem.Phys 110:3870].
%---------------------------------------------------------------

%% Do the propagation for all three values of Delta omega and
%% save the resulting expectation values.
global time expect

for ii=1:3
    qm_setup();
    qm_init(ii^2*100);
    qm_propa();
    ali{ii} = expect.amo{1}.cha{1};
    kin{ii} = expect.kin.cha{1};
    tim = time.steps.m_grid;
    qm_cleanup();
end

%% Plot it in the same way as in the paper.
figure(1);
clf;
plot(tim, ali{1}, 'k-', ...
     tim, ali{2}, 'k--', ...
     tim, ali{3}, 'k:', ...
     'LineWidth', 2);
xlabel('t');
ylabel('<cos^2 \theta>');
saveas(gcf, 'alignment.jpg');

figure(2);
clf;
plot(tim, kin{1}, 'k-', ...
     tim, kin{2}, 'k--', ...
     tim, kin{3}, 'k:', ...
     'LineWidth', 2);
xlabel('t');
ylabel('<J^2>');
saveas(gcf, 'j2.jpg');