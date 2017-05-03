% Script for getting the flux field for bead displacments with average
% displacement subtracted off
%clear, close all;

macpath = '/Volumes/Storage/'; % Path to server on mac
ubupath = '/mnt/llmStorage/'; % path to server on LLM ubuntu desktop
datapath = 'Danny/122311_5ulRED_2pt0uMact_phall_beads_nta_nochol_mc_1to300xlink_2/adding6pt2ulunspunskmuscmyo_5s/msd/';

if ismac
	path = [macpath datapath];
elseif isunix
	path = [ubupath datapath];
else
	error('Where is the data?')
end

%nbins = 50;%[25,15];
dt = 5;
rg_cutoff = [0,Inf]; % radius of gyration cutoffs (essentially none with [0,Inf])
tc = 60; % only look at time < 2/3 tc for now


[probMap, fluxField, xEdges,yEdges] = beadMotionFluxLoop(path,nbins,dt,rg_cutoff,tc);

dx.pre = xEdges.pre(2) - xEdges.pre(1);
dy.pre = yEdges.pre(2) - yEdges.pre(1);
binCentersX.pre = xEdges.pre(1:end-1)+dx.pre/2;
binCentersY.pre = yEdges.pre(1:end-1)+dy.pre/2;

netCurl.pre = sum(sum(curl(binCentersX.pre, binCentersY.pre,...
	fluxField.pre(:,:,1), fluxField.pre(:,:,2))))*dx.pre*dy.pre;

figure,
pcolor(xEdges.pre(1:end-1), yEdges.pre(1:end-1), probMap.pre)
colorbar;
hold on;

quiver(binCentersX.pre, binCentersY.pre,...
	fluxField.pre(:,:,1), fluxField.pre(:,:,2),'w')%, 'LineWidth', 1.5);
xlabel('\Deltax (\mum)')
ylabel('\Deltay (\mum)')
title(['Active systems, net curl=', num2str(netCurl.pre)])

% Do the control on a system without myosin
datapath ='/Danny/122311_5ulRED_2pt0uMact_phall_beads_nta_nochol_mc_1to300xlink_2/tc1_5s/msd/';

if ismac
	path = [macpath datapath];
elseif isunix
	path = [ubupath datapath];
else
	error('Where is the data on this platform?')
end

%nbins = 25;
dt = 5;
rg_cutoff = [0,Inf];
tc = []; 


[probMap.ctrl, fluxField.ctrl, xEdges.ctrl,yEdges.ctrl] = beadMotionFluxLoop(path,nbins,dt,rg_cutoff,tc);

dx.ctrl = xEdges.ctrl(2) - xEdges.ctrl(1);
dy.ctrl = yEdges.ctrl(2) - yEdges.ctrl(1);
binCentersX.ctrl = xEdges.ctrl(1:end-1)+dx.ctrl/2;
binCentersY.ctrl = yEdges.ctrl(1:end-1)+dy.ctrl/2;

netCurl.ctrl = sum(sum(curl(binCentersX.ctrl, binCentersY.ctrl, fluxField.ctrl(:,:,1), fluxField.ctrl(:,:,2))))*dx.ctrl*dy.ctrl;

figure,
pcolor(xEdges.ctrl(1:end-1), yEdges.ctrl(1:end-1), probMap.ctrl)
colorbar;
hold on;

quiver(binCentersX.ctrl, binCentersY.ctrl, fluxField.ctrl(:,:,1), fluxField.ctrl(:,:,2),'w', 'LineWidth', 1.5);
xlabel('\Deltax (\mum)')
ylabel('\Deltay (\mum)')
title(['Control, net curl=', num2str(netCurl.ctrl)])