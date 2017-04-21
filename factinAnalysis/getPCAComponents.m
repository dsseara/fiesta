% Script that uses PCA to track the state of the filaments over time
% Create a matrix of tangent angles, where each row is a single time point,
% and each column is a point along the arc length.

% Chosen ones for 560_lessthan_0.6tc(20170411T131818268)
%chosenOnes = [3,4,5,6,7,8,11,12,13,15,20,22,25,30,34,36,42,45,48,49,55,62];

% Chosen ones for 560_lessThan0.6tc(20170419T170302783)
% chosenOnes = [3,8,9,12,13,15,16,17,18,21,22,24,27,28,36,39,40,57,69,73,75,80,115,116,125];

pcaData = struct([]);

for jj = 1:numel(chosenOnes)
    filaID = chosenOnes(jj);
    figure
    pcaData(jj).filament = num2str(filaID);

    allTheta = Filament(filaID).theta;
    
    nFrames = size(allTheta,2);
    ntheta = cellfun(@numel, allTheta);
    mean_size = round(mean(ntheta));
    matrix = zeros(nFrames, mean_size);
    for ii = 1:nFrames
        theta = allTheta{ii};
        % Interpolate the theta vector to be the same size as the mean theta vector
        interpedTheta = interp1(1:ntheta(ii), theta, 1: (ntheta(ii)-1)/(mean_size-1) :ntheta(ii));
        matrix(ii,:) = interpedTheta; 
    end
    pcaData(jj).tangents = matrix;

    [pcaData(jj).coeff,pcaData(jj).score,pcaData(jj).latent,...
        pcaData(jj).tsquared,pcaData(jj).explained,pcaData(jj).mu]...
        = pca(pcaData(jj).tangents);

    colorline(pcaData(jj).score(:,1), pcaData(jj).score(:,2), 1:nFrames);
    title(['Filament ', num2str(filaID)]);
    xlabel('PCA component 1')
    ylabel('PCA component 2')

    saveas(gcf,['phaseSpacePlots' filesep 'pca' filesep 'fig' filesep 'pcaTimeSeries_fila_',num2str(filaID)],'fig')
    saveas(gcf,['phaseSpacePlots' filesep 'pca' filesep 'tif' filesep 'pcaTimeSeries_fila_',num2str(filaID)],'tif')
    saveas(gcf,['phaseSpacePlots' filesep 'pca' filesep 'eps' filesep 'pcaTimeSeries_fila_',num2str(filaID)],'epsc')
end