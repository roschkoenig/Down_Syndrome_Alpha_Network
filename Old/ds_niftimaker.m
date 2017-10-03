% Make Niftis
%==========================================================================
% (c) Richard Rosch 2016

% Define files and subjects to use
%--------------------------------------------------------------------------
clear all
p           = ds_definefiles;
sub         = p;
filter      = 'all';

% Housekeeping
%==========================================================================
fs  	= filesep;
if strcmp(computer, 'MACI64')
elseif strcmp(computer, 'PCWIN64');
    scripts = 'D:\Research_Data\1608 DS alpha';
    f_spm   = 'D:\My Documents\MATLAB\tools\spm';
	Fbase   = 'D:\Research_Data\1608 DS alpha';
    FEO     = [Fbase fs 'Data' fs 'EO processed 2-sec'];
    FEC     = [Fbase fs 'Data' fs 'EC processed 2-sec'];
    Fanalysis = [Fbase fs 'Matlab Files'];
end

addpath(f_spm);
addpath(scripts);
addpath(FEO);
addpath(FEC);
spm('defaults', 'eeg');


%%
clear sm fmap tmap map eo_map ec_map
for s = 1:length(sub)

disp(s);
S = ds_read(sub{s});

lbl     = S.eo_head.label;
tim_ax  = linspace(0,S.eo_head.nSamples/S.eo_head.Fs, S.eo_head.nSamples);
i       = 0;

% Extract time series for two conditions
%--------------------------------------------------------------------------
for o = 1:size(S.eo_data,3)
    i = i + 1;
    cond{i}         = 'Eyes Open';
    ftdata.trial{i} = squeeze(S.eo_data(:,:,o));
    ftdata.time{i}  = tim_ax;
end
for c = 1:size(S.ec_data,3)
    i = i + 1;
 	cond{i} = 'Eyes Closed';
    ftdata.trial{i} = squeeze(S.ec_data(:,:,c));
    ftdata.time{i}  = tim_ax;    
end    


% Calculate alpha bandpower maps for eyes open condition
% -------------------------------------------------------------------------
Fs  	= S.eo_head.Fs;
sens    = S.eo_head.elec;
xy      = spm_eeg_project3D(sens, 'EEG');
xpos    = xy(1,:); xpos = fix(xpos * 30);
ypos    = xy(2,:); ypos = fix(ypos * 30);

for r = 1 :size(S.eo_data,3)
for c = 1 :size(S.eo_data,1)
    tmap(r, xpos(c),ypos(c)) = bandpower(squeeze(S.eo_data(c,:,r)), Fs, [8 13]); 
end
end

mmap = squeeze(mean(tmap,1));
sm   = fspecial('gaussian', 10, 10);
eo_map(s,:,:)  = filter2(sm, mmap);

nii = make_nii(squeeze(eo_map(s,:,:)));
save_nii(nii, [sub{s}.name '_EO']);    clear nii;


% Calculate alpha bandpower maps for eyes closed condition
% -------------------------------------------------------------------------
sens    = S.ec_head.elec;
xy      = spm_eeg_project3D(sens, 'EEG');
xpos    = xy(1,:); xpos = fix(xpos * 30);
ypos    = xy(2,:); ypos = fix(ypos * 30);

clear map tmap fmap
Fs = S.eo_head.Fs;
for r = 1 :size(S.ec_data,3)
for c = 1 :size(S.ec_data,1)
    tmap(r, xpos(c),ypos(c)) = bandpower(squeeze(S.ec_data(c,:,r)), Fs, [8 13]); 
end
end

mmap = squeeze(mean(tmap,1));
sm        = fspecial('gaussian', 10, 10);
ec_map(s,:,:)  = filter2(sm, mmap);

nii = make_nii(squeeze(ec_map(s,:,:)));
save_nii(nii, [sub{s}.name '_EC']);    clear nii;
end



%% Plotting Routines
colormap gray
pl_range = [0 6];

subplot(2,1,1)
imagesc(squeeze(mean(ec_map,1)), pl_range);
set(gca, 'YDir', 'normal');
axis square

subplot(2,1,2)
imagesc(squeeze((eo_map(12,:,:))), pl_range); hold on
set(gca, 'YDir', 'normal');
axis square