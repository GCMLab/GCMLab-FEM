function plotLoadVsDispl(F, d, Control)
% PLOTLOADVSDISPL Plots displacement vs load curve
%   --------------------------------------------------------------------
%   Input
%   --------------------------------------------------------------------
%   F:  Matrix of externally applied loads (size nDOF x nsteps)
%   d: Matrix of displacements (size nDOF x nsteps)
%   Control: Structure array with the following field
%           .plotAt: selected DOF to plot curve (defined in Config File) 

% DOF to plot
plotAt = Control.plotAt;
d_plot = d(plotAt, :);
F_plot = F(plotAt, :);

% plot
figure;
plot(d_plot, F_plot, 'k-o', 'LineWidth', 1.5);
hold on
title('Load vs. Displacement curve');
xlabel('Displacement (m)');
ylabel('Load (N)');
hold off

end