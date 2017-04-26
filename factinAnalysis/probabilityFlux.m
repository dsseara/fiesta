% function [probMap, fluxField] = probabilityFluxLoop(x,y)
%
% This function takes a 2D time series and calculates the probability
% flux of the time series moving through the 2D space. Inspired by
% Battle et al. (Science 2016)
%
% INPUTS:     tseries - Nx2 matrix, each column is a time series in 
%                       one variable
%               nbins - Either an integer umber of bins in both directions,
%                       or a 1x2 array with nbins in x and y direction
%                  dt - time interval between each row of data in tseries
%              cutoff - (optional), if passed, an integer number of standard
%                       deviations out to which the bins go. Pass [] to just use
%                       max and min of data given.
%
% OUTPUTS:    probMap - (nbiny)x(nbinx) histogram of probability distribution of 
%                       phase space, estimated as time spent in each coarse
%                       grained bin in phase space 
%           fluxField - (nbiny)x(nbinx)x2 matrix, containing x and y components 
%                       of the time averaged flux vector field
%
% Created by Daniel Seara at 2017/04/25 10:29
function [probMap, fluxField,xEdges,yEdges] = probabilityFlux(tseries, dt, nbins, cutoff)
    
    switch numel(nbins)
    case 1
        nbinx = nbins;
        nbiny = nbins;
    case 2
        nbinx = nbins(1);
        nbiny = nbins(2);
    end

    % Find all the transitions that happens. Transitions matrix
    % will be (nbiny)(nbinx) long, so each state in the 2d state
    % is described by its linear index. The transition matrix has elements
    % t(i,j), which is a transition FROM STATE I TO STATE J
    transitions = zeros(nbinx*nbiny);

    if size(tseries,1)<size(tseries,2)
        % make sure data entered in right orientation
        tseries = tseries';
    end

    n = size(tseries,1); % total number of time points
    totalTime = dt*n;
    probMap = zeros(nbiny,nbinx);
    fluxField = zeros(nbiny, nbinx, 2);

    if ~isempty(cutoff)
        means = [mean(tseries(:,1)),mean(tseries(:,2))];
        stds  = [ std(tseries(:,1)), std(tseries(:,2))];

        minX = means(1)-cutoff*stds(1);
        maxX = means(1)+cutoff*stds(1);
        minY = means(2)-cutoff*stds(2);
        maxY = means(2)+cutoff*stds(2);

        % Use the +-Inf to catch outliers that lay outside mu +- cutoff*sigma
        xEdges = [-Inf, minX:(maxX - minX)/(nbinx-2):maxX,Inf];
        yEdges = [-Inf, minY:(maxY - minY)/(nbiny-2):maxY,Inf];
    else
        minX = min(tseries(:,1));
        maxX = max(tseries(:,1));
        minY = min(tseries(:,2));
        maxY = max(tseries(:,2));

        xEdges = (minX:(maxX-minX)/(nbinx):maxX);
        yEdges = (minY:(maxY-minY)/(nbiny):maxY);
    end

    for jj = 2:n
        % Get prior state
        px = tseries(jj-1,1);
        py = tseries(jj-1,2);
        pState = histcounts2(px,py, xEdges, yEdges);
        [pRow, pCol] = find(pState);
        
        % Get current state
        cx = tseries(jj,1);
        cy = tseries(jj,2);
        cState = histcounts2(cx,cy, xEdges, yEdges);
        [cRow, cCol] = find(cState);

        % Get all points traversed in phase space using Bresenham algorithm
        [xPath, yPath] = bresenham(pRow,pCol, cRow, cCol);
        % keep track of time spent at each state, divided equally between all 
        % the states on the path. Time between frames is 5 seconds
        
        % this part was tricky, switch order of path cuz x direction changes
        % in columns not rows, similar for y
        inds = sub2ind([nbiny, nbinx],yPath, xPath); 

        probMap(inds) = probMap(inds) + dt/numel(xPath);

        % If no transition, move on
        if numel(xPath)==1 || numel(yPath) ==1
            continue
        else
            % turn path coordinates into linear indices to put into transition matrix
            transitions(sub2ind(size(transitions), inds(1:end-1), inds(2:end)))...
                = transitions(sub2ind(size(transitions), inds(1:end-1), inds(2:end))) + 1;
        end
    end

    % This matrix can be multiplied against a local 3x3 matrix , then take
    % the trace, and you have the x-component of the flux
    xFluxMat = [-1/sqrt(2), -1, -1/sqrt(2); 0,0,0;1/sqrt(2), 1, 1/sqrt(2)];
    % Same for y component, call upwards the positive direction
    yFluxMat = -xFluxMat';

    % Now actually use th transitions to get the averaged flux
    for state = 1:nbinx*nbiny
        % Recall, transitions is built so that element (i,j) is
        % the number of transitions from i to j
        [stateRow,stateCol] = ind2sub([nbiny, nbinx], state);
        outFlow = reshape(transitions(state,:), [nbiny, nbinx]); % from state to all
        inFlow  = reshape(transitions(:,state), [nbiny, nbinx]); % from all to state
        netFlow = outFlow - inFlow;

        if stateRow ~= 1 && stateCol ~= 1 && stateRow ~= nbiny && stateCol ~= nbinx
            % center region %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            localFlow = netFlow(stateRow-1:stateRow+1,stateCol-1:stateCol+1);
            flowVec = 1/2.*[trace(localFlow*xFluxMat), trace(localFlow*yFluxMat)];

        elseif stateRow == 1 && stateCol ~= 1 && stateCol ~= nbinx
            % top row %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            localFlow = netFlow(stateRow:stateRow+1, stateCol-1:stateCol+1);
            flowVec = 1/2.*[trace(localFlow*xFluxMat(:,2:3)), trace(localFlow*yFluxMat(:,2:3))];

        elseif stateRow == nbiny && stateCol ~= 1 && stateCol ~= nbinx
            % bottom row %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            localFlow = netFlow(stateRow-1:stateRow, stateCol-1:stateCol+1);
            flowVec = 1/2.*[trace(localFlow*xFluxMat(:,1:2)), trace(localFlow*yFluxMat(:,1:2))];

        elseif stateCol == 1 && stateRow ~= 1 && stateRow ~= nbiny
            % left column %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            localFlow = netFlow(stateRow-1:stateRow+1,stateCol:stateCol+1);
            flowVec = 1/2.*[trace(localFlow*xFluxMat(2:3,:)), trace(localFlow*yFluxMat(2:3,:))];

        elseif stateCol == nbinx && stateRow ~= 1 && stateRow ~= nbiny
            % right columnm %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            localFlow = netFlow(stateRow-1:stateRow+1,stateCol-1:stateCol);
            flowVec = 1/2.*[trace(localFlow*xFluxMat(1:2,:)), trace(localFlow*yFluxMat(1:2,:))];

        elseif stateRow == 1 && stateCol == 1
            % northwest corner %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            localFlow = netFlow(stateRow:stateRow+1,stateCol:stateCol+1);
            flowVec = 1/2.*[trace(localFlow*xFluxMat(2:3,2:3)), trace(localFlow*yFluxMat(2:3,2:3))];

        elseif stateRow == 1 && stateCol == nbinx
            % northeast corner %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            localFlow = netFlow(stateRow:stateRow+1,stateCol-1:stateCol);
            flowVec = 1/2.*[trace(localFlow*xFluxMat(2:3,1:2)), trace(localFlow*yFluxMat(2:3,1:2))];

        elseif stateRow == nbiny && stateCol == 1
            % southwest corner %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            localFlow = netFlow(stateRow-1:stateRow,stateCol:stateCol+1);
            flowVec = 1/2.*[trace(localFlow*xFluxMat(1:2,2:3)), trace(localFlow*yFluxMat(1:2,2:3))];

        elseif stateRow == nbiny && stateCol == nbinx
            % southeast corner %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            localFlow = netFlow(stateRow-1:stateRow,stateCol-1:stateCol);
            flowVec = 1/2.*[trace(localFlow*xFluxMat(1:2,1:2)), trace(localFlow*yFluxMat(1:2,1:2))];

        end

        if ~isempty(find(flowVec,1))
            fluxField(stateRow, stateCol, :) = flowVec./probMap(stateRow,stateCol);
        else
            fluxField(stateRow, stateCol, :) = flowVec;
        end
    end

    probMap = probMap./totalTime;

