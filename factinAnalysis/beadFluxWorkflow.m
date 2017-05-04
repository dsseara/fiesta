% Script for getting the flux field for bead displacments with average
% displacement subtracted off
%clear, close all;

macpath = '/Volumes/Storage/'; % Path to server on mac
ubupath = '/mnt/llmStorage/'; % path to server on LLM ubuntu desktop

activeDataPath = 'Danny/122311_5ulRED_2pt0uMact_phall_beads_nta_nochol_mc_1to300xlink_2/adding6pt2ulunspunskmuscmyo_5s/msd/';

if ismac
	path = [macpath activeDataPath];
elseif isunix
	path = [ubupath activeDataPath];
else
	error('Where is the data on this platform?')
end

dt = 5;
rg_cutoff = [0,Inf]; % radius of gyration cutoffs (essentially none with [0,Inf])
tc = 50; % only look at time < 2/3 tc for now
dbin = 0.250; % Size of bins in um
cutoff = []; % number of times std to cut flux field off by, [] for no cutoff

[xActive, yActive] = beadAggregator(path,rg_cutoff,tc);

[probMap.active, fluxField.active, xEdges.active, yEdges.active] =...
                probabilityFlux([xActive, yActive], dt, dbin, cutoff);

binCentersX.active = xEdges.active(1:end-1)+dbin/2;
binCentersY.active = yEdges.active(1:end-1)+dbin/2;

netCurl.active = sum(sum(curl(binCentersX.active, binCentersY.active,...
	fluxField.active(:,:,1), fluxField.active(:,:,2))))*dbin^2;

figure,
pcolor(xEdges.active(1:end-1), yEdges.active(1:end-1), probMap.active)
colorbar;
hold on;

quiver(binCentersX.active, binCentersY.active,...
	fluxField.active(:,:,1), fluxField.active(:,:,2),'w')%, 'LineWidth', 1.5);
xlabel('$x- \langle x \rangle$ ($\mu$m)', 'Interpreter','latex')
ylabel('$y- \langle y \rangle$ ($\mu$m)', 'Interpreter','latex')
title(['Active system, net curl=', num2str(netCurl.active)])
saveas(gcf,['phaseSpace_active_tc_',num2str(tc)],'fig')
saveas(gcf,['phaseSpace_active_tc_',num2str(tc)],'tif')
saveas(gcf,['phaseSpace_active_tc_',num2str(tc)],'epsc')

% Do the control on a system without myosin
ctrlDataPath ='/Danny/122311_5ulRED_2pt0uMact_phall_beads_nta_nochol_mc_1to300xlink_2/tc1_5s/msd/';

if ismac
	path = [macpath ctrlDataPath];
elseif isunix
	path = [ubupath ctrlDataPath];
else
	error('Where is the data on this platform?')
end

dt = 5;
rg_cutoff = [0,Inf];
tc = [];

[xCtrl, yCtrl] = beadAggregator(path,rg_cutoff,tc);

[probMap.ctrl, fluxField.ctrl, xEdges.ctrl, yEdges.ctrl] =...
                probabilityFlux([xCtrl, yCtrl], dt, dbin, cutoff);

binCentersX.ctrl = xEdges.ctrl(1:end-1)+dbin/2;
binCentersY.ctrl = yEdges.ctrl(1:end-1)+dbin/2;

netCurl.ctrl = sum(sum(curl(binCentersX.ctrl, binCentersY.ctrl, fluxField.ctrl(:,:,1), fluxField.ctrl(:,:,2))))*dbin^2;

figure,
pcolor(xEdges.ctrl(1:end-1), yEdges.ctrl(1:end-1), probMap.ctrl)
colorbar;
hold on;

quiver(binCentersX.ctrl, binCentersY.ctrl, fluxField.ctrl(:,:,1), fluxField.ctrl(:,:,2),'w', 'LineWidth', 1.5);
xlabel('$x- \langle x \rangle$ ($\mu$m)', 'Interpreter','latex')
ylabel('$y- \langle y \rangle$ ($\mu$m)', 'Interpreter','latex')
title(['Control, net curl=', num2str(netCurl.ctrl)])
saveas(gcf,'phaseSpace_ctrl','fig')
saveas(gcf,'phaseSpace_ctrl','tif')
saveas(gcf,'phaseSpace_ctrl','epsc')