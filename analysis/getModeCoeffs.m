% Script to return a vector of coefficients of elastohydrodynamics modes of a
% slender rod for angle data

%chosenOnes = [3,4,7,32];

% Chosen ones for 560_lessthan_0.6tc(20170411T131818268)
chosenOnes = [49, 48, 45, 42, 34, 25, 20, 13, 11, 7, 6, 5, 4, 35, 22, 15, 3];


nmax = 3;
bcs = 'free'; % Only know free end stuff for now...
coeffs = struct([]);
N = numel(Filament);
for jj = 1:N
    %[~,ind] = ismember(jj, chosenOnes);
    coeffs(jj).filament = num2str(jj);

    allTheta = Filament(jj).theta;
    
    nframes = size(allTheta,2);
    a = zeros(3,nframes);
    
    for ii = 1:nframes
        theta = allTheta{ii};
        if max(size(theta))<50
            continue
        else
            L = Filament(jj).Results(ii,7);

            xy = Filament(jj).Data{ii}(:,1:2);
            dxy = diff(xy);
            ds = sqrt(dxy(:,1).^2 + dxy(:,2).^2);

            a(:,ii) = elastohydroModes(theta, ds, L, nmax, bcs);
        end
        
    end
    
    coeffs(jj).modeCoeffs = a;
end

% % Try with the Fourier coefficients
% % Say we have N+1 points
% % a_n = sqrt(2/L) sum(k=1 to N) theta_k * deltaS_k * cos(n pi s_k(mid) /L)
% % with s_k(mid) = sum(j=1 to k-1) deltaS_j + 1/2 deltaS_k
% % deltaS_k is the displacement from bead k to bead k+1

% for jj = 3
%     %coeffs(jj).filament = num2str(jj);

%     allTheta = Filament(jj).theta;
    
%     nframes = size(allTheta,2);
%     a = zeros(3,nframes);
    
%     for ii = 1:nframes
%         theta = allTheta{ii};
%         L = Filament(jj).Results(ii,7);
%         xy = Filament(jj).Data{ii}(:,1:2);
%         dxy = diff(xy);
%         ds = sqrt(dxy(:,1).^2 + dxy(:,2).^2);
%         for kk = 1:size(a,1)

%             a(kk,ii) = sqrt(2/L) * sum(theta.* ds .* cos(kk*pi*))
%         end
%     end

%     coeffs(jj).modeCoeffs = a;
% end