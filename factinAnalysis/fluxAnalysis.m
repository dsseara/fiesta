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

nbins = 25;

% Experiment Specific parameters
tframe = 5;
nFrames = 50;

% Collect data on how long the filaments spend in each state
time = zeros(nbins, nbins);

% place to store all the flux vectors
meanVector = zeros(nbins, nbins, 2);

% Split the data according to the set it comes from to
% accurately measure fluxes and time spent in each discrete state
data1 = aggregateData(aggregateData(:,end-1)==1,:);
filaID1 = unique(data1(:,end));
score1 = score(aggregateData(:,end-1)==1,:);

data2 = aggregateData(aggregateData(:,end-1)==2,:);
filaID2 = unique(data2(:,end));
score2 = score(aggregateData(:,end-1)==2,:);

pca1Edges = linspace(min(score(:,1)), max(score(:,1)), nbins+1);
pca2Edges = linspace(min(score(:,2)), max(score(:,2)), nbins+1);

% go through each filament separately for data 1
for ii = 1:numel(filaID1)
    filament = filaID1(ii);
    % Get time series for first and second modes
    pca1t = score1(data1(:,end)==filament,1);
    pca2t = score1(data1(:,end)==filament,2);

    % Find all the transitions that happen
    % Transitions matrix will be nbins^2 long, so each state in the 2d state
    % is described by its linear index. The transition matrix has elements
    % t(i,j), which is a transition FROM STATE I TO STATE J
    transitions = zeros(nbins^2);

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

        % Get all points traversed in phase space
        [pca1Path, pca2Path] = bresenham(pRow,pCol, cRow, cCol);
        % keep track of time spent at each state, divided equally between all the states
        % on the path
        % time between frames is 5 seconds
        time(sub2ind(size(time), pca1Path, pca2Path)) =...
            time(sub2ind(size(time), pca1Path, pca2Path)) + tframe/numel(pca1Path);

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

    for state = 1:nbins^2
        % Recall, transitions is built so that element (i,j) is
        % the number of transitions from i to j
        [stateRow,stateCol] = ind2sub([nbins, nbins], state);
        outFlow = reshape(transitions(state,:), [nbins, nbins]); % from state to all
        inFlow  = reshape(transitions(:,state), [nbins, nbins]); % from all to state
        netFlow = inFlow - outFlow;
        % The flux vector is given by mean of the value of netflow
        % (positive, negative, or zero) times the distance 
        distX = (1:nbins) - stateCol;
        distY = (1:nbins) - stateRow;

        % want to multiply every row by each element in distX
        vx = mean(mean(netFlow .* repmat(distX',[1,nbins])));

        % want to multiply every col by each element in distY
        vy = mean(mean(repmat(distY, [nbins,1]) .* netFlow));
        
        meanVector(stateRow, stateCol, :,ii) = [vx, vy];
        
    end
end

% go through each filament separately for data 2
for ii = 1:numel(filaID2)
    filament = filaID2(ii);
    % Get time series for first and second modes
    pca1t = score2(data2(:,end)==filament,1);
    pca2t = score2(data2(:,end)==filament,2);

    % Find all the transitions that happen
    % Transitions matrix will be nbins^2 long, so each state in the 2d state
    % is described by its linear index. The transition matrix has elements
    % t(i,j), which is a transition FROM STATE I TO STATE J
    transitions = zeros(nbins^2);

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

        % Get all points traversed in phase space
        [pca1Path, pca2Path] = bresenham(pRow,pCol, cRow, cCol);
        % keep track of time spent at each state, divided equally between all the states
        % on the path
        % time between frames is 5 seconds
        time(sub2ind(size(time), pca1Path, pca2Path)) =...
            time(sub2ind(size(time), pca1Path, pca2Path)) + tframe/numel(pca1Path);

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

    for state = 1:nbins^2
        % Recall, transitions is built so that element (i,j) is
        % the number of transitions from i to j
        [stateRow,stateCol] = ind2sub([nbins, nbins], state);
        outFlow = reshape(transitions(state,:), [nbins, nbins]); % from state to all
        inFlow  = reshape(transitions(:,state), [nbins, nbins]); % from all to state
        netFlow = inFlow - outFlow;
        % The flux vector is given by mean of the value of netflow
        % (positive, negative, or zero) times the distance 
        distX = (1:nbins) - stateCol;
        distY = (1:nbins) - stateRow;

        % want to multiply every row by each element in distX
        vx = mean(mean(netFlow .* repmat(distX',[1,nbins])));

        % want to multiply every col by each element in distY
        vy = mean(mean(repmat(distY, [nbins,1]) .* netFlow));
        
        meanVector(stateRow, stateCol, :,ii) = [vx, vy];
        
    end
end

% Plot the vectors over the phase space
figure, hold on
xlim([pca1Edges(1),pca1Edges(end)])
ylim([pca2Edges(1), pca2Edges(end)])
% shift the edges a bit so that they are in the center of the bins
binCenters = [pca1Edges(1:end-1) + diff(pca1Edges)/2; pca2Edges(1:end-1)+diff(pca2Edges)/2];
% full time of experiment is tFrame*nFrames
pcolor(binCenters(1,:), binCenters(2,:), time(:,:,ii)./250) 
% shading interp
colorbar
quiver(pca1Edges(2:end), pca2Edges(2:end), meanVector(:,:,1,ii), meanVector(:,:,2,ii), 1.5, 'w', 'LineWidth', 1.5 )

title(['Filament ', num2str(filaID)]);
xlabel('PCA component 1')
ylabel('PCA component 2')

% saveas(gcf,['phaseSpacePlots' filesep 'pca' filesep 'fig' filesep 'fluxMap_fila_',num2str(filaID)],'fig')
% saveas(gcf,['phaseSpacePlots' filesep 'pca' filesep 'tif' filesep 'fluxMap_fila_',num2str(filaID)],'tif')
% saveas(gcf,['phaseSpacePlots' filesep 'pca' filesep 'eps' filesep 'fluxMap_fila_',num2str(filaID)],'epsc')