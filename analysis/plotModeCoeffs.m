% Plot mode coefficients of beams

% First get an idea of what places in phase space are visited by all the
% filaments

% Chosen ones for 560_lessthan_0.6tc(20170411T131818268)
chosenOnes = [49, 48, 45, 42, 34, 25, 20, 13, 11, 7, 6, 4, 35, 22, 15, 3]; %5
% Turns out that filament 5 has a couple of wildly changing phase space points, so we omit it here

figure

colors = colormap(hsv(numel(Filament)));

for ii = chosenOnes
    nframes = size(coeffs(ii).modeCoeffs,2);
    scatter(coeffs(ii).modeCoeffs(2,:), coeffs(ii).modeCoeffs(3,:),...%coeffs(ii).modeCoeffs(2,:),... %
        'o', 'DisplayName', num2str(ii)); %'Color', colors(ii,:)
    hold on;
end
legend('show')

xlabel('a2')
ylabel('a3')
title('Phase space scatter plot')

saveas(gcf,['phaseSpacePlots' filesep 'scatterPlot_a2a3'],'fig')
saveas(gcf,['phaseSpacePlots' filesep 'scatterPlot_a2a3'],'tif')
saveas(gcf,['phaseSpacePlots' filesep 'scatterPlot_a2a3'],'epsc')

%%
% Now we bin them to create a coarse-grained phase space
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

counts = zeros(numel(a1Edges)-1,numel(a2Edges)-1);

for ii = chosenOnes
    a1 = coeffs(ii).modeCoeffs(1,:);
    a2 = coeffs(ii).modeCoeffs(2,:);
    a3 = coeffs(ii).modeCoeffs(3,:);

    counts = counts + histcounts2(a1,a2, a1Edges,a2Edges);
end
figure
pcolor(a1Edges,a2Edges,padarray(counts,[1,1],'post'))
title('counts')
xlabel('a1')
ylabel('a2')
colorbar
saveas(gcf,['phaseSpacePlots' filesep 'histogram'],'fig')
saveas(gcf,['phaseSpacePlots' filesep 'histogram'],'tif')
saveas(gcf,['phaseSpacePlots' filesep 'histogram'],'epsc')
figure
pcolor(a1Edges,a2Edges,padarray(log(counts),[1,1],'post'))
title('log(counts)')
xlabel('a1')
ylabel('a2')
colorbar
saveas(gcf,['phaseSpacePlots' filesep 'log_histogram'],'fig')
saveas(gcf,['phaseSpacePlots' filesep 'log_histogram'],'tif')
saveas(gcf,['phaseSpacePlots' filesep 'log_histogram'],'epsc')

