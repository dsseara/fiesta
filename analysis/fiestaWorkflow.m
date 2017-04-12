% Workflow to analyze FIESTA data

%%
% First find what fiesta data you want, and call the .mat file name
% fname. For example:
% fname = '560_lessthan_0.6tc(20170404T173034436)'

run plotFilaments.m

%%
% Now rotate the filaments and plot them along with their angles to
% get an idea of what's going on

run rotateFilaments.m

%%
% Now get the elastohydrodynamic mode coefficients.
% This depends on the function 'elastohydroModes.m'

run getModeCoeffs.m

%%
% Now plot and bin the data in a coarse-grained phase space

run plotModeCoeffs.m

%%
% Now find the fluxes from one coarse-grained phase space to another

run fluxAnalysis.m