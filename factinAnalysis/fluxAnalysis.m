% Flux analysis now calculates the PCA components, and performs the
% flux analysis all in one go on the aggregated data that is in a matrix
% called aggregateData
%
% All data on server
%
% Data set 1
% /Danny/122311_5ulRED_2pt0uMact_phall_beads_nta_nochol_mc_1to300xlink_2/adding6pt2ulunspunskmuscmyo_5s/560/gaussSmoothed_2px/560_lessthan_0.6tc(20170411T131818268)
%
% Data set 2
% /Danny/122311_5ulRED_2pt0uMact_phall_beads_nta_nochol_mc_1to300xlink/adding6pt2ulunspunskmuscmyo_5s/560/gaussSmoothed_2px/560_lessThan0.6tc(20170419T170302783)

[coeff,score,latent,tsquared,explained,mu]=...
    pca(aggregateData, 'NumComponents', 3);

nbins = 100;
nbinx = nbins;
nbiny = nbins;

% Experiment Specific parameters
tframe = 5;
nFrames = 50;

% Collect data on how long the filaments spend in each state
time1 = zeros(nbins, nbins);
time2 = zeros(nbins, nbins);

% place to store all the flux vectors
fluxField = zeros(nbins, nbins, 2);

% Split the data according to the set it comes from to
% accurately measure fluxes and time spent in each discrete state
data1 = aggregateData(aggregateData(:,end-1)==1,:);
filaID1 = unique(data1(:,end));
score1 = score(aggregateData(:,end-1)==1,:);

data2 = aggregateData(aggregateData(:,end-1)==2,:);
filaID2 = unique(data2(:,end));
score2 = score(aggregateData(:,end-1)==2,:);

% Cut histogram off by mean and std in each direction
means = [mean(score(:,1)),mean(score(:,2))];
stds  = [ std(score(:,1)), std(score(:,2))];

minEdge1 = means(1)-2*stds(1);
maxEdge1 = means(1)+2*stds(1);
minEdge2 = means(2)-2*stds(2);
maxEdge2 = means(2)+2*stds(2);

% Use the +-Inf to catch outliers that lay outside mu +- 3*sigma
pca1Edges = [-Inf, minEdge1:(maxEdge1 - minEdge1)/(nbins-3):maxEdge1,Inf];
pca2Edges = [-Inf, minEdge2:(maxEdge2 - minEdge2)/(nbins-3):maxEdge2,Inf];

%pca1Edges = linspace(min(score(:,1)), max(score(:,1)), nbins+1);
%pca2Edges = linspace(min(score(:,2)), max(score(:,2)), nbins+1);

% Find all the transitions that happen
% Transitions matrix will be nbins^2 long, so each state in the 2d state
% is described by its linear index. The transition matrix has elements
% t(i,j), which is a transition FROM STATE I TO STATE J
transitions = zeros(nbins*nbins);

% go through each filament separately for data 1 (location given above)
for ii = 2%1:numel(filaID1)
    filament = filaID1(ii);
    % Get time series for first and second modes
    pca1t = score1(data1(:,end)==filament,1);
    pca2t = score1(data1(:,end)==filament,2);

    for jj = 2:numel(pca1t)
        % Get prior state
        pca1P = pca1t(jj-1);
        pca2P = pca2t(jj-1);
        priorState   = histcounts2(pca1P, pca2P, pca1Edges, pca2Edges);
        [pRow, pCol] = find(priorState);
        
        % Get current state
        pca1C = pca1t(jj);
        pca2C = pca2t(jj);
        currentState = histcounts2(pca1C, pca2C, pca1Edges, pca2Edges);
        [cRow, cCol] = find(currentState);

        % Get all points traversed in phase space using Bresenham algorithm
        [pca1Path, pca2Path] = bresenham(pRow,pCol, cRow, cCol);
        % keep track of time spent at each state, divided equally between all the states
        % on the path. Time between frames is 5 seconds
        time1(sub2ind(size(time1), pca1Path, pca2Path)) =...
            time1(sub2ind(size(time1), pca1Path, pca2Path)) + tframe/numel(pca1Path);

        % If no transition, move on
        if numel(pca1Path)==1 || numel(pca2Path) ==1
            continue
        else
            % turn path coordinates into linear indices to put into transition matrix
            linearInd = sub2ind([nbins nbins], pca1Path, pca2Path);
            transitions(sub2ind(size(transitions), linearInd(1:end-1), linearInd(2:end)))...
                = transitions(sub2ind(size(transitions), linearInd(1:end-1), linearInd(2:end))) + 1;
        end
    end
end

% % go through each filament separately for data 2 (location given above)
% for ii = 1:numel(filaID2)
%     filament = filaID2(ii);
%     % Get time series for first and second modes
%     pca1t = score2(data2(:,end)==filament,1);
%     pca2t = score2(data2(:,end)==filament,2);

%     for jj = 2:numel(pca1t)
%         % Get prior state
%         pca1P = pca1t(jj-1);
%         pca2P = pca2t(jj-1);
%         priorState   = histcounts2(pca1P, pca2P, pca1Edges, pca2Edges);
%         [pRow, pCol] = find(priorState);
        
%         % Get current state
%         pca1C = pca1t(jj);
%         pca2C = pca2t(jj);
%         currentState = histcounts2(pca1C, pca2C, pca1Edges, pca2Edges);
%         [cRow, cCol] = find(currentState);

%         % Get all points traversed in phase space using Bresenham algorithm
%         [pca1Path, pca2Path] = bresenham(pRow,pCol, cRow, cCol);
%         % keep track of time spent at each state, divided equally between all the states
%         % on the path. Time between frames is 5 seconds
%         time2(sub2ind(size(time2), pca1Path, pca2Path)) =...
%             time2(sub2ind(size(time2), pca1Path, pca2Path)) + tframe/numel(pca1Path);

%         % If no transition, move on
%         if numel(pca1Path)==1 || numel(pca2Path) ==1
%             continue
%         else
%             % turn path coordinates into linear indices to put into transition matrix
%             linearInd = sub2ind([nbins nbins], pca1Path, pca2Path);
%             transitions(sub2ind(size(transitions), linearInd(1:end-1), linearInd(2:end)))...
%                 = transitions(sub2ind(size(transitions), linearInd(1:end-1), linearInd(2:end))) + 1;
%         end
%     end
% end

time = (time1 + time2);

% This matrix can be multiplied against the localFlow matrix, then take
% the trace, and you have the x-component of the flux
xFluxMat = [-1/sqrt(2), -1, -1/sqrt(2); 0,0,0;1/sqrt(2), 1, 1/sqrt(2)];
% Same for y component, call upwards the positive direction
yFluxMat = -xFluxMat';

for state = 1:nbins^2
    % Recall, transitions is built so that element (i,j) is
    % the number of transitions from i to j
    [stateRow,stateCol] = ind2sub([nbins, nbins], state);
    outFlow = reshape(transitions(state,:), [nbins, nbins]); % from state to all
    inFlow  = reshape(transitions(:,state), [nbins, nbins]); % from all to state
    netFlow = outFlow - inFlow;

    if stateRow ~= 1 && stateCol ~= 1 && stateRow ~= nbins && stateCol ~= nbins
        % center region
        localFlow = netFlow(stateRow-1:stateRow+1,stateCol-1:stateCol+1);
        flowVec = 1/2.*[trace(localFlow*xFluxMat), trace(localFlow*yFluxMat)];
    elseif stateRow==1 && stateCol ~=1 && stateCol ~= nbins
        % top row
        localFlow = netFlow(stateRow:stateRow+1, stateCol-1:stateCol+1);
        flowVec = 1/2.*[trace(localFlow*xFluxMat(:,2:3)), trace(localFlow*yFluxMat(:,2:3))];
    elseif stateRow==nbins && stateCol ~=1 && stateCol ~= nbins
        % bottom row
        localFlow = netFlow(stateRow-1:stateRow, stateCol-1:stateCol+1);
        flowVec = 1/2.*[trace(localFlow*xFluxMat(:,1:2)), trace(localFlow*yFluxMat(:,1:2))];
    elseif stateCol==1 && stateRow ~=1 && stateRow ~= nbins
        % left column
        localFlow = netFlow(stateRow-1:stateRow+1,stateCol:stateCol+1);
        flowVec = 1/2.*[trace(localFlow*xFluxMat(2:3,:)), trace(localFlow*yFluxMat(2:3,:))];
    elseif stateCol==nbins && stateRow ~=1 && stateRow ~= nbins
        % right column
        localFlow = netFlow(stateRow-1:stateRow+1,stateCol-1:stateCol);
        flowVec = 1/2.*[trace(localFlow*xFluxMat(1:2,:)), trace(localFlow*yFluxMat(1:2,:))];
    elseif stateRow==1 && stateCol==1
        % northwest corner
        localFlow = netFlow(stateRow:stateRow+1,stateCol:stateCol+1);
        flowVec = 1/2.*[trace(localFlow*xFluxMat(2:3,2:3)), trace(localFlow*yFluxMat(2:3,2:3))];
    elseif stateRow==1 && stateCol==nbins
        % northeast corner
        localFlow = netFlow(stateRow:stateRow+1,stateCol-1:stateCol);
        flowVec = 1/2.*[trace(localFlow*xFluxMat(2:3,1:2)), trace(localFlow*yFluxMat(2:3,1:2))];
    elseif stateRow==nbins && stateCol == 1
        % southwest corner
        localFlow = netFlow(stateRow-1:stateRow,stateCol:stateCol+1);
        flowVec = 1/2.*[trace(localFlow*xFluxMat(1:2,2:3)), trace(localFlow*yFluxMat(1:2,2:3))];
    elseif stateRow==nbins && stateCol == nbins
        % southeast corner
        localFlow = netFlow(stateRow-1:stateRow,stateCol-1:stateCol);
        flowVec = 1/2.*[trace(localFlow*xFluxMat(1:2,1:2)), trace(localFlow*yFluxMat(1:2,1:2))];
    end
    % if isempty(find(localFlow,1))
    %     flux = [0, 0];
    % else
    %     [row,col,val] = find(localFlow);
    %     flux = [(col-stateCol),(row-stateRow)].*val;
    % end
    % The flux vector is given by mean of the value of netflow
    % (positive, negative, or zero) times the distance 
    % distX = (1:nbins) - stateCol;
    % distY = (1:nbins) - stateRow;

    % % want to multiply every row by each element in distX
    % vx = mean(mean(netFlow .* repmat(distX',[1,nbins])));

    % % want to multiply every col by each element in distY
    % vy = mean(mean(repmat(distY, [nbins,1]) .* netFlow));
    % if ~isempty(find(netFlow,1))
    %     keyboard
    % end
    % vx = mean(sum(netFlow,2));
    % vy = mean(sum(netFlow,1));
    if ~isempty(find(flowVec,1))
    	%keyboard
	    fluxField(stateRow, stateCol, :) = flowVec./time(stateRow,stateCol);
	else
		fluxField(stateRow, stateCol, :) = flowVec;
	end
end

% Plot the vectors over the phase space
figure, hold on
dx = pca1Edges(3)-pca1Edges(2);
dy = pca2Edges(3)-pca2Edges(2);
xlim([pca1Edges(1), pca1Edges(end)])
ylim([pca2Edges(1), pca2Edges(end)])
% shift the edges a bit so that they are in the center of the bins
binCenters = [pca1Edges(1:end-1) + diff(pca1Edges)/2; pca2Edges(1:end-1)+diff(pca2Edges)/2];
% full time of experiment is tFrame*nFrames
% pcolor(binCenters(1,:), binCenters(2,:), time(:,:)./(tframe*size(score,1))); 
pcolor(pca1Edges,pca2Edges, time./(tframe*size(score,1)));
% shading interp
colorbar
%quiver(pca1Edges(2:end), pca2Edges(2:end), fluxField(:,:,1), fluxField(:,:,2), 'w', 'LineWidth', 1.5 )
quiver(pca1Edges+(dx./2), pca2Edges+(dy./2), fluxField(:,:,1), fluxField(:,:,2), 'w', 'LineWidth', 1.5 )
title('Filament phase space');
xlabel('PCA component 1')
ylabel('PCA component 2')

% saveas(gcf,['phaseSpacePlots' filesep 'pca' filesep 'fig' filesep 'fluxMap_fila_',num2str(filaID)],'fig')
% saveas(gcf,['phaseSpacePlots' filesep 'pca' filesep 'tif' filesep 'fluxMap_fila_',num2str(filaID)],'tif')
% saveas(gcf,['phaseSpacePlots' filesep 'pca' filesep 'eps' filesep 'fluxMap_fila_',num2str(filaID)],'epsc')