% Script that uses PCA to track the state of the filaments over time
% Create a matrix of tangent angles, where each row is a single time point,
% and each column is a point along the arc length.

% Chosen ones for 560_lessthan_0.6tc(20170411T131818268)
chosenOnes = [3,4,5,6,7,11,13,15,20,22,25,34,35,42,45,48,49];

pcaData = struct([]);

for jj = chosenOnes
    figure
    [~,ind] = ismember(jj,chosenOnes);
    pcaData(ind).filament = num2str(jj);

    allTheta = Filament(jj).theta;
    
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
    pcaData(ind).tangents = matrix;
    [pcaData(ind).coeff,pcaData(ind).score,pcaData(ind).latent,...
        pcaData(ind).tsquared,pcaData(ind).explained,pcaData(ind).mu]...
        = pca(pcaData(ind).tangents);

    colorline(pcaData(ind).score(:,1), pcaData(ind).score(:,2), 1:nFrames);
    title(['Filament ', num2str(jj)]);
    xlabel('PCA component 1')
    ylabel('PCA component 2')

    saveas(gcf,['pcaTimeSeries_fil_',num2str(jj)],'fig')
    saveas(gcf,['pcaTimeSeries_fil_',num2str(jj)],'tif')
    saveas(gcf,['pcaTimeSeries_fil_',num2str(jj)],'epsc')
end