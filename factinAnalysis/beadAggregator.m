% function [xall, yall] = beadAggregator(path,rg_cutoff, tc)
%
% This program aggregates all the x and y positions of beads tracked using 
% Maria Kilfoil's MATLAB particle tracking software into two arrays
%
% INPUTS      path - the base path for the experiment. Reads in the individual
%                    beads files and the correspondance matrix with Rg from the 
%                    "Bead_Tracking\ddposum_files\individual_beads\" subfolder.
%        rg_cutoff - [min max] of radius of gyration to use
%               tc - (optional) critical time (as a frame number) to cut time
%                    into.
%
% OUTPUTS:    xall - column vector of all mean subtracted x positions over time
%             yall - column vector of all mean subtracted y positions over time
%
% Created by Daniel Seara at 2017/04/25 17:42

function [xall, yall] = beadAggregator(path,rg_cutoff, tc)

%load([path 'Bead_Tracking\ddposum_files\individual_beads\correspondance'])
if ispc
    load([path 'Bead_Tracking\ddposum_files\individual_beads\correspondance_RG'])
elseif isunix
    load([path 'Bead_Tracking/ddposum_files/individual_beads/correspondance_RG'])
end
xall = [];
yall = [];

for ii = 1:length(correspondance(:,1))
    
    if ispc
        load([path 'Bead_Tracking\ddposum_files\individual_beads\bead_' num2str(ii)]);
    elseif isunix
        load([path 'Bead_Tracking/ddposum_files/individual_beads/bead_' num2str(ii)]);
    end

    if correspondance(ii,4) < rg_cutoff(2) && correspondance(ii,4)>rg_cutoff(1)
        if isempty(tc)
            bsecx=(bsec(:,1)-bsec(1,1));
            x = bsecx - mean(bsecx);
            
            bsecy=(bsec(:,2)-bsec(1,2));
            y = bsecy - mean(bsecy);

            xall = [xall;x];
            yall = [yall;y];
        else
            %%% Pre tc %%%
            pre  = bsec(bsec(:,3)<tc+1,:);
            
            if isempty(pre)
                % msdx.pre = 0;
                % msdy.pre = 0;
                % msd.pre  = 0;
                disp('empty pre')
                continue
            end

            pre_bsecx=(pre(:,1)-pre(1,1));
            pre_x = pre_bsecx - mean(pre_bsecx);
            
            pre_bsecy=(pre(:,2)-pre(1,2));
            pre_y = pre_bsecy - mean(pre_bsecy);

            xall = [xall;pre_x];
            yall = [yall;pre_y];

            %%% Post tc %%%
            % post  = bsec(bsec(:,3)>tc,:);
            
            % if isempty(post)
            %     % msdx.post = 0;
            %     % msdy.post = 0;
            %     % msd.post  = 0;
            %     disp('empty post')
            %     continue
            % end

            % post_bsecx=(post(:,1)-post(1,1));
            % post_x = post_bsecx - mean(post_bsecx);
            
            % post_bsecy=(post(:,2)-post(1,2));
            % post_y = post_bsecy - mean(post_bsecy);

            % [probMap.post, fluxField.post,xEdges.post, yEdges.post] =...
            %     probabilityFlux([post_x, post_y], dt, dbin, []);
        end % end if statement that parses whether you have a tc
    end % end if statement that gets beads with correct rg_cutoff

    if mod(ii,50) == 0
        disp(['Finished computing msd for bead number ' num2str(i)])
    end
end % end for loop over all particles