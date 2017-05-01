% Cleaning filament angle data for easy analysis
% Have two data sets in the following locations, followed by the filaments I want
% to keep as well as what needs to be done for each of them to clean everything up.
% The "cleaning" really just refers to taking out tracked filaments that differ
% significantly from the rest of the data for that filament.
% 
% Data set 1
%/Volumes/Storage/Danny/122311_5ulRED_2pt0uMact_phall_beads_nta_nochol_mc_1to300xlink_2/adding6pt2ulunspunskmuscmyo_5s/560/gaussSmoothed_2px/560_lessthan_0.6tc(20170411T131818268)
%
% Filaments:
% 3,4,5,6,7,8,11,12,13,15,20,22,25,30,34,36,42,45,48,49,55,62
% 
% Data set 2
% /Volumes/Storage/Danny/122311_5ulRED_2pt0uMact_phall_beads_nta_nochol_mc_1to300xlink/adding6pt2ulunspunskmuscmyo_5s/560/gaussSmoothed_2px/560_lessThan0.6tc(20170419T170302783)
%
% Filaments:
% 3,8,9,12,13,15,16,17,18,21,22,24,27,28,36,39,40,57,69,73,75,80,115,116,125

% Data set 3
% Danny/122311_5ulRED_2pt0uMact_phall_beads_nta_nochol_mc_1to300xlink_2/tc1_5s/560(20170425T191553204)/560(20170425T191553204).mat
%
% Filaments:
% 
% 7,9,10,16,21,28,32,34,36,40,41,44,51,52,61,63,67,70,73,83,94,101,107,118
%
% Set a cut-off size for the segments, and see how many we're left with!
filaList = {[3,4,5,6,7,8,11,12,13,15,20,22,25,30,34,36,42,45,48,49,55,62],...
			[3,8,9,12,13,15,16,17,18,21,22,24,27,28,36,39,40,57,69,73,75,80,115,116,125],...
			[7, 9, 10, 16, 21, 28, 32, 34, 36, 40, 41, 44, 51, 52, 61, 63, 67, 70, 73, 83, 94, 101, 107, 118]};

% 2nd to last column of aggregateData tells you which data set, last column tells
% which filament inside that data set.
aggregateData = [];

% % First have to cd to the server
% load('Danny/122311_5ulRED_2pt0uMact_phall_beads_nta_nochol_mc_1to300xlink_2/adding6pt2ulunspunskmuscmyo_5s/560/gaussSmoothed_2px/560_lessthan_0.6tc(20170411T131818268)/560_lessthan_0.6tc(20170411T131818268).mat');
% 
% for ii = filaList{1}
% 	thetas = Filament(ii).theta;
% 	thetas = thetas(cellfun(@(x) size(x,2)>cutoff,thetas));
% 	for jj = 1:numel(thetas)
% 		aggregateData = [aggregateData; [thetas{jj}(1:cutoff),1, ii]];
% 	end
% end
% 
% load('Danny/122311_5ulRED_2pt0uMact_phall_beads_nta_nochol_mc_1to300xlink/adding6pt2ulunspunskmuscmyo_5s/560/gaussSmoothed_2px/560_lessThan0.6tc(20170419T170302783)/560_lessThan0.6tc(20170419T170302783).mat');
% 
% for ii = filaList{2}
% 	thetas = Filament(ii).theta;
% 	thetas = thetas(cellfun(@(x) size(x,2)>cutoff,thetas));
% 	for jj = 1:numel(thetas)
% 		aggregateData = [aggregateData; [thetas{jj}(1:cutoff),2, ii]];
% 	end
% end

load('Danny/122311_5ulRED_2pt0uMact_phall_beads_nta_nochol_mc_1to300xlink_2/tc1_5s/560(20170425T191553204)/560(20170425T191553204).mat')

% Chosen ones for 560(20170425T191553204)
for ii = filaList{3}
	thetas = Filament(ii).theta;
	thetas = thetas(cellfun(@(x) size(x,2)>cutoff,thetas));
	for jj = 1:numel(thetas)
		aggregateData = [aggregateData; [thetas{jj}(1:cutoff),2, ii]];
	end
end