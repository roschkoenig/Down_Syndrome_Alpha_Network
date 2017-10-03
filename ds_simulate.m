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

% Load all files into single cell array
%--------------------------------------------------------------------------
RCM = load([Fdcm fs 'Reduced_Models.mat']);
ACM = RCM.RCM;

% Calculate Bayesian Parameter Average as starting point for simulations
%--------------------------------------------------------------------------
BMA = spm_dcm_bma(ACM');
TCM = ACM{1};

%% Parameter definitions for simulations
%--------------------------------------------------------------------------
Ps      = BMA.Ep;
steps   = 100;
rnge    = 5 * [-1 1];

% MFG
%--------------------------------------------------------------------------
    G{1} = [0 0 0 0 0 0 0 0;   % ltV1 
            0 0 0 0 0 0 0 0;   % rtV1 
            0 0 0 0 0 0 0 0;   % ltSPL 
            0 0 0 0 0 0 0 0;   % rtSPL 
            1 0 0 0 0 0 0 0;   % ltMFg 
            1 0 0 0 0 0 0 0];  % rtMFg 
        
% SPL
%--------------------------------------------------------------------------        
    G{2} = [0 0 0 0 0 0 0 0;  % ltV1 
            0 0 0 0 0 0 0 0;  % rtV1 
            1 0 0 0 0 0 0 0;  % ltSPL 
            1 0 0 0 0 0 0 0;  % rtSPL 
            0 0 0 0 0 0 0 0;  % ltMFG 
            0 0 0 0 0 0 0 0]; % rtMFg    


% Calculate linearly spaced parameter combinations
%--------------------------------------------------------------------------
g1       = linspace(rnge(1), rnge(2), steps); 
g2       = linspace(rnge(1), rnge(2), steps+1);
gname   = {'MFG self-inhibition', 'SPL self-inhibition'};

textprogressbar(['Simulating ' num2str(steps) ' steps:   ']);
clear Hc al al1 al2

% Simulate effects of parameter affecting all levels of hierarchy 
%--------------------------------------------------------------------------
count = 0;
for i1 = 1:steps
for i2 = 1:steps+1
    count = count + 1;
    textprogressbar(count * 100 / (steps * (steps+1)));
    TPs         = Ps;
    TPs.G       = Ps.G + g1(i1) * G{1} + g2(i2) * G{2};

    TempHc      = spm_csd_mtf(TPs, TCM.M, TCM.xU);
    Hc{i1,i2}   = TempHc{1};
end
end
textprogressbar(' Done');

%% Plotting routines
%==========================================================================
franges{1} = [4:7];
franges{2} = [8:13];
fnames     = {'theta', 'alpha'};
figure(1)

for f = 1:length(franges)
clear al
cols    = flip(cbrewer('div', 'Spectral', 100));
colormap(cols)

% Extract feature of interest (e.g. alpha peak)
%--------------------------------------------------------------------------
frange  = franges{f};
m       = 1;

for i1 = 1:steps
for i2 = 1:steps+1
    al(i1,i2)   =   log(abs(mean(Hc{i1,i2}(frange,m,m))));
end 
end

for i1 = 1:steps,   al1(i1,:) = (abs(Hc{i1, fix(steps/2)}(:,m,m))); end     % effect of dimension 1 - frontal inh
for i2 = 1:steps+1, al2(i2,:) = (abs(Hc{fix(steps/2), i2}(:,m,m))); end     % effect of dimension 2 - spl inh 


% Plot effect across all levels of hierarchy
%-------------------------------------------------------------------------- 
subplot(1,2,f)
    
% Plotting
    imagesc(g1, g2, al), hold on
    
    % Labels
    title(['Effects of ' FX ' on on mean ' fnames{f}], 'Fontsize', 15);
    xlabel(gname{2});
    ylabel(gname{1});
    
    % Settings
    axis square
    set(gca, 'YDir', 'normal');
    set(gcf, 'Position', [200 200 1400 800]);

set(gcf, 'color', 'w');
end


for a = 1:length(ACM)
    g(a,:) = ACM{a}.Ep.G(:,1);
end
figure(1)
scatter(g(:,3), g(:,5), 'filled');

ylim([-2 1]);
xlim([-4 4]);

%%
figure(2)
cols = flip(cbrewer('div', 'Spectral', size(al1,1)));
colormap(cols);
for a = 1:size(al1,1)
    subplot(1,2,1)
    plot(al1(a,:), 'color', cols(a,:)); hold on
    
	subplot(1,2,2)
    plot(al2(a,:), 'color', cols(a,:)); hold on
end
