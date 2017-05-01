% Plot FIESTA output
close all

if ~exist('tif','dir')
    mkdir('tif')
end

if ~exist('eps','dir')
    mkdir('eps')
end

if ~exist('fig','dir')
    mkdir('fig')
end

nFrames = size(Objects,2);
colors = colormap(parula(nFrames));
hold on
for jj = 1:nFrames
    data = Objects{jj}.data;
    nFila = size(data,2);
    for ii = 1:nFila
        
        if ~isempty(data{ii})
            filaData = data{ii}(:,1:2);
            if size(filaData,1)<50
                continue
            end
        else
            continue
        end
        plot(filaData(:,1)./1000,filaData(:,2)./1000, 'Color', colors(jj,:));
    end

end
title('Objects');
saveas(gcf, ['fig' filesep 'objects'], 'fig');
saveas(gcf, ['tif' filesep 'objects'], 'tif');
saveas(gcf, ['eps' filesep 'objects'], 'epsc');


% This plots individual filaments separately
nFila = size(Filament,2);
for jj = 1:nFila
    data = Filament(jj).Data;
    nFrames = size(data,2);
    colors = colormap(parula(nFrames));
    disp(['Filament number: ', num2str(jj),'. Number of frames: ', num2str(nFrames)])
    figure, hold on
    title(['Filament ', num2str(jj)])
    for ii = 1:nFrames
        % if ~isempty(data{ii})
            filaData = data{ii}(:,1:2);
        % else
        %     continue
        % end
        plot(filaData(:,1)./1000,filaData(:,2)./1000, 'Color', colors(ii,:));
    end
    saveas(gcf,['fig' filesep 'filament_',num2str(jj)],'fig')
    saveas(gcf,['tif' filesep 'filament_',num2str(jj)],'tif')
    saveas(gcf,['eps' filesep 'filament_',num2str(jj)],'epsc')
end