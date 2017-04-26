% function [probMap, fluxField] = beadMotionFluxLoop(path, nbins, dt,rg_cutoff, tc)
%
% This program calculates the flux field and probability
% map of particles tracked with Maria Kilfoil's particle
% tracking software. Just a hacked version of the
% many particle msd calculation
%
% INPUTS      path - the base path for the experiment. Reads in the individual
%                    beads files and the correspondance matrix with Rg from the 
%                    "Bead_Tracking\ddposum_files\individual_beads\" subfolder.
%            nbins - Either an integer umber of bins in both directions,
%                    or a 1x2 array with nbins in x and y direction
%               dt - the (average) delta-t between frames, in seconds
%        rg_cutoff - [min max] of radius of gyration to use
%               tc - (optional) critical time (as a frame number) to cut time
%                    into.
%
% OUTPUTS:    probMap - (nbiny)x(nbinx) histogram of probability distribution of 
%                       phase space, estimated as time spent in each coarse
%                       grained bin in phase space 
%           fluxField - (nbiny)x(nbinx)x2 matrix, containing x and y components 
%                       of the time averaged flux vector field
%
% Created by Daniel Seara at 2017/04/25 17:42

function [probMap, fluxField, xEdges, yEdges] = beadMotionFluxLoop(path, nbins, dt,rg_cutoff, tc)

%load([path 'Bead_Tracking\ddposum_files\individual_beads\correspondance'])
if ispc
    load([path 'Bead_Tracking\ddposum_files\individual_beads\correspondance_RG'])
elseif isunix
    load([path 'Bead_Tracking/ddposum_files/individual_beads/correspondance_RG'])
end

% if isempty(tc)
    
% else
%     %%% Pre tc %%%
%     pre_number_of_frames = tc;
%     probMap.pre=zeros(nbinx, nbiny);
%     %%% Post tc %%%
%     post_number_of_frames = number_of_frames - tc;
%     probMap.post=[zeros(nbinx, nbiny);
%     post_pointtracer=zeros(post_number_of_frames-1,1);

%     % Make two structs, one for pre and one for post tc
%     msd.pre=zeros(pre_number_of_frames-1,1);
%     msdx.pre=zeros(pre_number_of_frames-1,1);
%     msdy.pre=zeros(pre_number_of_frames-1,1);

%     msd.post=zeros(post_number_of_frames-1,1);
%     msdx.post=zeros(post_number_of_frames-1,1);
%     msdy.post=zeros(post_number_of_frames-1,1);
% end

for i = 1:length(correspondance(:,1))
    
    if ispc
        load([path 'Bead_Tracking\ddposum_files\individual_beads\bead_' num2str(i)]);
    elseif isunix
        load([path 'Bead_Tracking/ddposum_files/individual_beads/bead_' num2str(i)]);
    end
    
    if correspondance(i,4) < rg_cutoff(2) && correspondance(i,4)>rg_cutoff(1)
        if isempty(tc)
            bsecx=(bsec(:,1)-bsec(1,1));
            meanDX = mean(diff(bsecx));
            dx = diff(bsecx) - meanDX;
            
            bsecy=(bsec(:,2)-bsec(1,2));
            meanDY = mean(diff(bsecy));
            dy = diff(bsecy) - meanDY;

            [probMap, fluxField, xEdges, yEdges] = probabilityFlux([dx, dy], dt, nbins, []);
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
            pre_meanDX = mean(diff(pre_bsecx));
            pre_dx = diff(pre_bsecx) - pre_meanDX;
            
            pre_bsecy=(pre(:,2)-pre(1,2));
            pre_meanDY = mean(diff(pre_bsecy));
            pre_dy = diff(pre_bsecy) - pre_meanDY;

            [probMap.pre, fluxField.pre,xEdges.pre, yEdges.pre] = probabilityFlux([pre_dx, pre_dy], dt, nbins, []);

            %%% Post tc %%%
            post  = bsec(bsec(:,3)<tc+1,:);
            
            if isempty(post)
                % msdx.post = 0;
                % msdy.post = 0;
                % msd.post  = 0;
                disp('empty post')
                continue
            end

            post_bsecx=(post(:,1)-post(1,1));
            post_meanDX = mean(diff(post_bsecx));
            post_dx = diff(post_bsecx) - post_meanDX;
            
            post_bsecy=(post(:,2)-post(1,2));
            post_meanDY = mean(diff(post_bsecy));
            post_dy = diff(post_bsecy) - post_meanDY;

            [probMap.post, fluxField.post, xEdges.pre, yEdges.post] = probabilityFlux([post_dx, post_dy], dt, nbins, []);
        end
    end
    if mod(i,50) == 0
        disp(['Finished computing msd for bead number ' num2str(i)])
    end
end
