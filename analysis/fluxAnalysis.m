% Calculate the flux from one point in phase space to another over time
clearvars -except fname
close all;

load(fname)

% make chosenOnes an array before running, depending on FIESTA file you want

% 560_lessthan_0.6tc(20170411T131818268)
chosenOnes = [3,4,5,6,7,11,13,15,20,22,25,34,35,42,45,48,49];

nbins = 25;

% Experiment Specific parameters
tframe = 5;
nFrames = 50;


% Collect data on how long the filament spends in each state
time = zeros(nbins, nbins, numel(chosenOnes));

% place to store all the flux vectors
meanVector = zeros(nbins, nbins, 2, numel(chosenOnes));

% Start loop over all the filaments
for ii = 1:numel(chosenOnes)
    filaID = chosenOnes(ii);

    % Get time series for first and second modes
    pca1t = pcaData(ii).score(:,1);
    pca2t = pcaData(ii).score(:,2);
    
    nbins = 25;

    pca1Edges = linspace(min(pca1t), max(pca1t), nbins+1);
    pca2Edges = linspace(min(pca2t), max(pca2t), nbins+1);
    
    
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

        % Try just building an algorithm based entirely around the Bresenham Line algorithm.
        % Give the state prior and state after as row and column, and then interpolate all
        % the states that go between them, and count them all

        % Get all points traversed in phase space
        [pca1Path, pca2Path] = bresenham(pRow,pCol, cRow, cCol);
        % keep track of time spent at each state, divided equally between all the states
        % on the path
        % time between frames is 5 seconds
        time(sub2ind(size(time), pca1Path, pca2Path, repmat(ii,size(pca1Path)))) =...
            time(sub2ind(size(time), pca1Path, pca2Path, repmat(ii,size(pca1Path)))) + tframe/numel(pca1Path);

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

    % Get the flux vectors for every state in the phase space
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

        % want to multiple every row by each element in distX
        vx = mean(mean(netFlow .* repmat(distX',[1,nbins])));

        % want to multiple every col by each element in distY
        vy = mean(mean(repmat(distY, [nbins,1]) .* netFlow));
        
        meanVector(stateRow, stateCol, :,ii) = [vx, vy];
        
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
end



