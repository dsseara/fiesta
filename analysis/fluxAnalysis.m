% Calculate the flux from one point in phase space to another over time

% Chosen ones for 560_lessthan_0.6tc(20170411T131818268)
chosenOnes = [49, 48, 45, 42, 34, 25, 20, 13, 11, 7, 6, 4, 35, 22, 15, 3]; %5,
% Turns out that filament 5 has a couple of wildly changing phase space points, so we omit it here

for ii = chosenOnes
    maxa1(ii) = max(coeffs(ii).modeCoeffs(1,:));
    maxa2(ii) = max(coeffs(ii).modeCoeffs(2,:));
    maxa3(ii) = max(coeffs(ii).modeCoeffs(3,:));
    
    mina1(ii) = min(coeffs(ii).modeCoeffs(1,:));
    mina2(ii) = min(coeffs(ii).modeCoeffs(2,:));
    mina3(ii) = min(coeffs(ii).modeCoeffs(3,:));
end

maxa1 = max(maxa1);
mina1 = min(mina1);
maxa2 = max(maxa2);
mina2 = min(mina2);
maxa3 = max(maxa3);
mina3 = min(mina3);

nbins = 25;
a1Edges = linspace(mina1,maxa1,nbins+1);
a2Edges = linspace(mina2,maxa2,nbins+1);
a3Edges = linspace(mina3,maxa3,nbins+1);

% Transitions matrix will be nbins^2 long, so each state in the 2d state
% is described by its linear index. The transition matrix has elements
% t(i,j), which is a transition FROM STATE I TO STATE J
transitions = zeros(nbins^2);

for ii = chosenOnes
    a1t = coeffs(ii).modeCoeffs(1,:);
    a2t = coeffs(ii).modeCoeffs(2,:);

    for jj = 2:numel(a1t)
        a1Prior   = a1t(jj-1);
        a1Current = a1t(jj);
        a2Prior   = a2t(jj-1);
        a2Current = a2t(jj);

        priorState   = histcounts2(  a1Prior,   a2Prior, a1Edges, a2Edges);
        [pRow, pCol] = find(priorState);
        pLinearInd = find(priorState);

        currentState = histcounts2(a1Current, a2Current, a1Edges, a2Edges);
        [currentRow, currentCol] = find(currentState);
        currentLinearInd = find(currentState);

        % Try just building an algorithm based entirely around the Bresenham Line algorithm.
        % Give the state prior and state after as row and column, and then interpolate all
        % the states that go between them, and count the flux going from all the states there.

        

    end
end



