% Housekeeping
%==========================================================================
clear all
D           = ds_housekeeping;
Fbase       = D.Fbase;
Fscripts    = D.Fscripts;
Fanalysis   = D.Fanalysis;
Fdata       = D.Fdata;
Fdcm        = D.Fdcm;
p           = ds_definefiles(Fbase);
fs          = filesep;
load([Fscripts fs 'SubInfo']);


%% Load all subject DCMs into single file array
%==========================================================================
for s = 1:length(sub)
    TCM = load([Fdcm fs 'DCM_' sub{s}.name '.mat']);
    ACM{s} = TCM.DCM; 
    ecal(s) = sub{s}.ecal;
end

% Rank by alpha power
%--------------------------------------------------------------------------
[val ind] = sort(ecal);


%% Extract spectra and plot for each subjects 
%==========================================================================
spectra = zeros(length(ecal), 30);
modes = 3;

for m = 1:modes
    
cols = cbrewer('div', 'Spectral', length(ecal));
colsort = cols(ind,:);

for s = 1:length(sub)
    prd(s,:) = (abs(ACM{s}.Hc{1}(:,m,m)));
    obs(s,:) = (abs(ACM{s}.xY.y{1}(:,m,m))); 

    set(gcf, 'Color', 'w');
 
%--------------------------------------------------------------------------    
subplot(2,modes,m)
%--------------------------------------------------------------------------
    title('Predicted Values', 'Fontweight', 'bold');
    plot(prd(s,:), 'Color', colsort(s,:), 'Linewidth', 1.5); hold on
    
    % Titles and Labels
    ylabel('Power');
    xlabel('Frequency');
    
%--------------------------------------------------------------------------
subplot(2,modes,modes+m)
%--------------------------------------------------------------------------
    title('Empirical Values', 'Fontweight', 'bold');
    plot(obs(s,:), 'Color', colsort(s,:), 'Linewidth', 1.5); hold on
    
	% Titles and Labels
    ylabel('Power');
    xlabel('Frequency');
    
end
end