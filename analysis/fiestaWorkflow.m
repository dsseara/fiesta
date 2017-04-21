% Workflow to analyze FIESTA data

%%
% First find what fiesta data you want, and call the .mat file name
% fname. For example:
% fname = '560_lessthan_0.6tc(20170404T173034436)'

run plotFilaments.m

%%
% Now rotate the filaments and plot them along with their angles to
% get an idea of what's going on
% Remember to save the new Filaments structure! Adds a new field called
% "theta" which has all the tangent information

run rotateFilaments.m

%%
% Now clean up the filament data and then perform PCA on the aggregated
% data. First have to set a cutoff for number of segments in the filaments
% to consider

cutoff = 50;


%%
% Now plot and bin the data in a coarse-grained phase space

run plotModeCoeffs.m

%%
% Now find the fluxes from one coarse-grained phase space to another

run fluxAnalysis.m