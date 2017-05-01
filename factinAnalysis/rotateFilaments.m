% Rotate filaments and get the tangent angle data
% Filaments chosen by eye to be good for analysis
%chosenOnes = [3,4,7,32];%,43,39];

% Chosen ones for 560_lessthan_0.6tc(20170411T131818268)
%chosenOnes = [49, 48, 45, 42, 34, 25, 20, 13, 11, 7, 6, 5, 4, 35, 22, 15, 3];

% Chosen ones for 560(20170425T191553204)
% chosenOnes = [7, 9, 10, 16, 21, 28, 32, 34, 36, 40, 41, 44, 51, 52, 61, 63, 67, 70, 73, 83, 94, 101, 107, 118];

N = numel(Filament);

for jj = 1:N
    data = Filament(jj).Data;

    % Find angle from horizontal of line that joins first and last point in 20th frame
    theta = atan2(data{1}(end,2) - data{1}(1,2), data{1}(end,1) - data{1}(1,1));

    % rotate clockwise by angle theta
    R = [cos(theta), sin(theta);-sin(theta), cos(theta)]; 


    figure;


    colors = colormap(parula(size(data,2)));
    
    nframes = size(data,2);
    Filament(jj).theta = cell(1,nframes);

    for ii = 1:nframes
        temp = [data{ii}(:,1)';data{ii}(:,2)'];

        npoints = size(temp,2);

        L = Filament(jj).Results(ii,7);
        positions = -L/2:(L/(npoints-2)):L/2;
        
        % Find center of mass
        % com = sum(temp,2)/size(temp,2);
        
        % Rotate and subtract center of mass
        temp = R*(temp - [temp(1,1).*ones(1,size(temp,2));temp(2,1).*ones(1,size(temp,2))]);
        % [com(1).*ones(1,size(temp,2));com(2).*ones(1,size(temp,2))]);

        % Get tangent angle information
        displacements = diff(temp, 1,2); % take difference along columns
        tangents = atan2(displacements(2,:), displacements(1,:));
        Filament(jj).theta{ii} = tangents;

        % Plot
        subplot(2,1,1), hold on
        plot(temp(1,:), temp(2,:), 'Color', colors(ii,:));
        subplot(2,1,2), hold on
        plot(positions, tangents, 'Color', colors(ii,:))

    end

    subplot(2,1,1);
    xlabel('x (nm)')
    ylabel('y (nm)')
    title(['Rotated and anchored filament ', num2str(jj)])
    subplot(2,1,2);
    xlabel('s(nm)')
    ylabel('\theta')
    title(['\theta(s) of filament ', num2str(jj)])
%     if ismember(jj,chosenOnes)
%         saveas(gcf,['~/Desktop/angleInfo/fig' filesep 'filament', num2str(jj),'_angleInfo'],'fig')
%         saveas(gcf,['~/Desktop/angleInfo/tif' filesep 'filament', num2str(jj),'_angleInfo'],'tif')
%         saveas(gcf,['~/Desktop/angleInfo/eps' filesep 'filament', num2str(jj),'_angleInfo'],'epsc')
%     end
end

