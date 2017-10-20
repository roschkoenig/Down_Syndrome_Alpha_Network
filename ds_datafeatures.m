%% Frequency Analysis for Down Syndrome spont EEG data
%--------------------------------------------------------------------------
% This routine will call on a number of different function to perform
% simple analysis of spectral differences between different participants
% with Down syndrome. 
%==========================================================================
% (c) Richard Rosch 2016

%% Define files and subjects to use
%-------------------------------------------------------------------------
clear all

% Housekeeping
%==========================================================================
fs  	= filesep;
D       = ds_housekeeping;
Fbase   = D.Fbase;
FEO     = [Fbase fs 'Data' fs 'EO processed 2-sec'];
FEC     = [Fbase fs 'Data' fs 'EC processed 2-sec'];

addpath(FEO);
addpath(FEC);
spm('defaults', 'eeg');

p           = ds_definefiles(Fbase);
sub         = p;
filter      = 'all';

% Loop through individual subjects and calculate spectra
%==========================================================================
clear eo ec
EL{1} = {'E74', 'E75', 'E70', 'E71', 'E76', 'E83', 'E82'};     % occipital
EL{2} = {'E11', 'E12', 'E19', 'E18', 'E16', 'E10', 'E4', 'E5'}; % frontal
names = {'Occipital', 'Frontal'};

for E = 1:length(EL)
for s = 1:length(sub)
    
% Load and Filter data
%--------------------------------------------------------------------------
% Data loading into structure, e.g.: 
% S.eo_data(nChans x nSamples x nTrials)

clear eoft ecft EO_Oz EC_Oz 
disp(s);

Fs  = 250;
S   = ds_read(sub{s});
el  = EL{E};

for e = 1:length(el)
    idx             = find(strcmp(el{e}, S.ec_head.label));
    EC_Oz(e,:,:)    = squeeze(S.ec_data(idx, :, :));
end

% Estimate fourier transforms of of data from single electrodes
%--------------------------------------------------------------------------
count = 0;
for t = 1:size(EC_Oz, 3)
for e = 1:size(EC_Oz, 1)
    count   = count + 1;
    ft      = abs(fft(squeeze(EC_Oz(e,:,t)), Fs));
    ecft(count,:) = ft(1:fix(end/2));
    clear ft
end
end

% Generate mean spectral composition across time windows within subject
% --------------------------------------------------------------------------
ec_all{E}(s,:) = log(mean(ecft.^2));

end
end

% Plotting routines for Full Group Spectra
%==========================================================================
for E = 1:length(EL)
ec = ec_all{E};
% Calculate standard errors and plot Eyes Closed
%--------------------------------------------------------------------------
for e = 1:size(ec,2)
    sem(e) = std(ec(:,e))/sqrt(length(ec(:,e)));
end

m_ec      = mean(ec, 1);
ec_upper  = m_ec + 1.96*sem;
ec_lower  = m_ec - 1.96*sem;

% Plot overall group spectra and standard errors
%--------------------------------------------------------------------------
figure(1)
subplot(1,length(EL), E)
    % Plots
    pec = plot(1:length(m_ec), m_ec, 'r'); hold on
    plotshaded(1:length(m_ec), [ec_upper; ec_lower], 'r')
    plot(1:length(m_ec), ec_upper, 'r');
    plot(1:length(m_ec), ec_lower, 'r');
    % Labe;s
    ylabel('Power');
    xlabel('Frequency');
    title(names{E});
    % Settings
    xlim([1 30]);
    
% Median Split by K-bit scores and compare spectra 
%==========================================================================
clear k a c
for s = 1:length(sub)
    works(s) = 1;
    try k(s) = sub{s}.Kbit; catch works(s) = 0; end
    try a(s) = sub{s}.age;  catch works(s) = 0; end
end

wi  = find(works);
mk  = median(k);
hi  = k > mk;
lo  = k <= mk;

c{1} = ec(find(hi), :);
c{2} = ec(find(lo), :);

% Estimate means and standard errors by condition
%--------------------------------------------------------------------------
for i = 1:length(c)
    for e = 1:size(c{i},2)
        sem(e) = std(c{i}(:,e))/sqrt(length(c{i}(:,e)));
    end
    m      = mean(c{i}, 1);
    upper  = m + sem;
    lower  = m - sem;
    
    % Define colors for conditions
    %----------------------------------------------------------------------
    cols = cbrewer('qual','Set1', 9);    
    if mod(i,2) == 1, col = cols(1,:);
    else col = cols(2,:);
    end
    
    % Plot
    %----------------------------------------------------------------------    
    figure(2)
    subplot(1,length(EL), E);
    % Plots
    pkbit(i) = plot(1:length(m), m, 'Color', col); hold on;
%     plot(1:length(m), upper, 'color', col);
%     plot(1:length(m), lower, 'color', col);
    plotshaded(1:length(m), [upper; lower], col) 
    % Labels
    xlabel('Frequency', 'fontweight', 'bold');
    ylabel('Power', 'fontweight', 'bold')
    % Settings
    xlim([1 30]);
%     ylim([0 1300]);
    axis square
end
end

le1 = legend([pkbit(1) pkbit(2)], {'High K-bit', 'Low K-bit'}); 
set(gcf, 'color', 'w')


%% Linear Regression Model (Peak Power)
%==========================================================================
for E = 1:length(ec_all)
fqr = [4 13];
for s = 1:length(sub)
    works(s) = 1;
    try a(s) = sub{s}.age;  catch, works(s) = 0;    end
    try k(s) = sub{s}.Kbit; catch, works(s) = 0;    end
    
    ecal(s) = max( (ec_all{E}(s, fqr(1):fqr(end))) );
    sub{s}.ecal = ecal(s);
end
wi      = find(works);
a       = a(wi);
k       = k(wi);
ecal    = exp(ecal(wi));

% Orthogonalise Age and Kbit
%--------------------------------------------------------------------------
P       = polyfit(k,a,1);
ares    = a - P(1)*k - P(2);
tbl     = table(ares', k', log(ecal)', 'VariableNames', {'Age', 'Kbit', 'EC_Peak'});
lmec 	= fitlm(tbl, 'EC_Peak ~ Kbit')

figure
plot(lmec)
axis square
title(names{E});
set(gcf, 'color', 'w');
save('SubInfo', 'sub');

end

% %% Median Split by K-bit scores and compare spectra 
% %==========================================================================
% clear k a c
% for s = 1:length(sub)
%     try k(s) = sub{s}.Kbit; end
%     try a(s) = sub{s}.age; end
% end
% 
% % Median split by K-bit score
% %--------------------------------------------------------------------------
% mk  = median(k);
% hi  = k > mk;
% lo  = k <= mk;
% 
% ma  = median(a);
% oa  = a > ma;
% ya  = a <= ma;
% 
% c{1} = eo(intersect(find(hi), find(ya)), :);
% c{2} = eo(intersect(find(lo), find(ya)), :);
% c{3} = eo(intersect(find(hi), find(oa)), :);
% c{4} = eo(intersect(find(lo), find(oa)), :);
% 
% c{5} = ec(intersect(find(hi), find(ya)), :);
% c{6} = ec(intersect(find(lo), find(ya)), :);
% c{7} = ec(intersect(find(hi), find(oa)), :);
% c{8} = ec(intersect(find(lo), find(oa)), :);
% 
% % Estimate means and standard errors by condition
% %--------------------------------------------------------------------------
% for i = 1:length(c)
%     for e = 1:size(c{i},2)
%         sem(e) = std(c{i}(:,e))/sqrt(length(c{i}(:,e)));
%     end
%     m      = mean(c{i}, 1);
%     upper  = m + sem;
%     lower  = m - sem;
% 
%     % Define subplot index for conditions
%     %----------------------------------------------------------------------
%     if      i <= 2, plidx = 1;
%     elseif  i <= 4, plidx = 2;
%     elseif  i <= 6, plidx = 3;
%     elseif  i <= 8, plidx = 4;
%     end
%     
%     % Define colors for conditions
%     %----------------------------------------------------------------------
%     cols = cool(5);    
%     if mod(i,2) == 1, col = cols(2,:);
%     else col = cols(5,:);
%     end
%     
%     % Plot
%     %----------------------------------------------------------------------    
%     subplot(2,2,plidx)
%     pkbit(i) = plot(1:length(m), m, 'Color', col); hold on;
%     plotshaded(1:length(m), [upper; lower], col) 
%     xlabel('Frequency');
%     ylabel('Power')
%     xlim([1 30]);
%     axis square
% end
% 
% subplot(2,2,1), legend([pkbit(1) pkbit(2)], {'High K-bit', 'Low K-bit'}); title('Eyes open, young age');
% subplot(2,2,2), legend([pkbit(3) pkbit(4)], {'High K-bit', 'Low K-bit'}); title('Eyes open, old age');
% subplot(2,2,3), legend([pkbit(5) pkbit(6)], {'High K-bit', 'Low K-bit'}); title('Eyes closed, young age');
% subplot(2,2,4), legend([pkbit(7) pkbit(8)], {'High K-bit', 'Low K-bit'}); title('Eyes closed, old age');
% 
% set(gcf, 'color', 'w')
% 
% % Median Split by K-bit - alpha power and peak frequency
% %==========================================================================
% clear alpha
% fqr = [8 13];
% 
% cond{1}.i     = intersect(find(hi), find(oa));    cond{1}.name = 'OA, HI';
% cond{2}.i     = intersect(find(lo), find(oa));    cond{2}.name = 'OA, LO';
% cond{3}.i     = intersect(find(hi), find(ya));    cond{3}.name = 'YA, HI';
% cond{4}.i     = intersect(find(lo), find(ya));    cond{4}.name = 'YA, LO';
% 
% for c = 1:length(cond)
% i = 0;
% for h = cond{c}.i
%     i = i+1;
%     alpha{c}.d(i) = mean(ec(h,fqr(1):fqr(end))); alpha{c}.name = cond{c}.name;
% end
% end
% subplot(2,1,1)
%     p = ranksum([alpha{1}.d, alpha{3}.d], [alpha{2}.d, alpha{4}.d]);
%     ds_dotplot(alpha);
%     title('Eyes closed, grouped by Kbit (median split) and Age (median split)');
%     xlabel(['Wilcoxon between K-bit: ' num2str(p)]);
%     axis square
% 
% 
% for c = 1:length(cond)
% i = 0;
% for h = cond{c}.i
%     i = i+1;
%     alpha{c}.d(i) = mean(eo(h,fqr(1):fqr(end))); alpha{c}.name = cond{c}.name;
% end
% end
% 
% subplot(2,1,2)
%     p = ranksum([alpha{1}.d, alpha{3}.d], [alpha{2}.d, alpha{4}.d]);
%     ds_dotplot(alpha);
%     xlabel(['Wilcoxon between K-bit: ' num2str(p)]);
%     title('Eyes open, grouped by Kbit (median split) and Age (median split)');
%     axis square
%     set(gcf, 'color', 'w');
% 
% 
% 
% 
% 
