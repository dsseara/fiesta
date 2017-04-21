% Workflow to analyze FIESTA data

%%
% First find and load the fiesta data you want

run plotFilaments.m

%%
% Now rotate the filaments and plot them along with their angles to
% get an idea of what's going on
% Remember to save the new Filaments structure! Adds a new field called
% "theta" which has all the tangent information

run rotateFilaments.m

%%
% Now clean up the filament data and aggregate them
% First have to set a cutoff for number of segments in the filaments
% to consider

cutoff = 50;
run cleanFilamentAngleData.m

%%
% Now, flux analysis performs PCA on the aggregated data
% and finds the fluxes from one coarse-grained phase space to another

run fluxAnalysis.m