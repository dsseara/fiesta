% Rotate filaments and get the tangent angle data
% Filaments chosen by eye to be good for analysis
%chosenOnes = [3,4,7,32];%,43,39];

% Chosen ones for 560_lessthan_0.6tc(20170411T131818268)
% chosenOnes = [3,4,5,6,7,11,13,15,20,22,25,34,35,42,45,48,49];

% Chosen ones for 560_lessThan0.6tc(20170419T170302783)
chosenOnes = [3, 8, 12, 13, 15, 16, 18, 21, 27, 36, 57, 69, 73,...
    ];

N = numel(Filament);

for jj = 1:N
    data = Filament(jj).Data;

    % Find angle from horizontal of line that joins first and last point in 20th frame
    theta = atan2(data{1}(end,2) - data{1}(1,2), data{1}(end,1) - data{1}(1,1));

    % rotate clockwise by angle theta
    R = [cos(theta), sin(theta);-sin(theta), cos(theta)]; 

    figure, hold on;
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
        
        % Rotate and subtract starting position
        temp = R*(temp - [temp(1,1).*ones(1,size(temp,2));temp(2,1).*ones(1,size(temp,2))]);
        % Replace 2nd term above with this to subtract center of mass
        % [com(1).*ones(1,size(temp,2));com(2).*ones(1,size(temp,2))]);

        % Get tangent angle information
        displacements = diff(temp, 1,2); % take difference along columns
        tangents = atan2(displacements(2,:), displacements(1,:));
        Filament(jj).theta{ii} = tangents;

        % Plot
        subplot(2,1,1), hold on
        plot(temp(1,:), temp(2,:), 'Color', colors(ii,:));
        subplot(2,1,2), hold on
        plot(temp(1,1:end-1), tangents, 'Color', colors(ii,:))

    end

    subplot(2,1,1);
    xlabel('x (nm)')
    ylabel('y (nm)')
    title(['Rotated and anchored filament ', num2str(jj)])
    subplot(2,1,2);
    xlabel('s(nm)')
    ylabel('\theta')
    title(['\theta(s) of filament ', num2str(jj)])
    if ismember(jj,chosenOnes)
        saveas(gcf,['fig', filesep, 'filament', num2str(jj),'_angleInfo'],'fig')
        saveas(gcf,['tif', filesep, 'filament', num2str(jj),'_angleInfo'],'tif')
        saveas(gcf,['eps', filesep, 'filament', num2str(jj),'_angleInfo'],'epsc')
    end
end

