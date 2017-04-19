% Workflow to analyze FIESTA data

%%
% First find and load the fiesta data you want

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