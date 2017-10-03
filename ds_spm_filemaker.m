% Make SPM files
%==========================================================================
% (c) Richard Rosch 2016

% Housekeeping
%==========================================================================
clear all
fs  	= filesep;
D       = ds_housekeeping;

Fbase       = D.Fbase;
Fscripts    = D.Fscripts;
Fanalysis   = D.Fanalysis;
Fdata       = D.Fdata;

% Define files and subjects to use
%--------------------------------------------------------------------------
p           = ds_definefiles(Fbase);
sub         = p;
filter      = 'all';

FEO         = [Fbase fs 'Data' fs 'EO processed 2-sec'];
FEC         = [Fbase fs 'Data' fs 'EC processed 2-sec'];
   
addpath(FEO);
addpath(FEC);
spm('defaults', 'eeg');


%% Loop through individual subjects and generate SPM files
%==========================================================================
for s = 1:length(sub)

% Load and Filter data
%--------------------------------------------------------------------------
% Data loading into structure, e.g.: 
% S.eo_data(nChans x nSamples x nTrials)

disp(s);
[S no_open] = ds_read(sub{s});


if no_open,     length_eo   = 0; 
else            length_eo   = length(S.eo_head.label); 
end

length_ec = length(S.ec_head.label);


ls      = [length_eo, length_ec];
[v i]   = min(ls);

% Identify label indices that are common between the two conditions
%--------------------------------------------------------------------------
if no_open
    lbl = S.ec_head.label;
    eoidx = [];
    ecidx = [];
    for l = 1:length(lbl)
        ecidx = [ecidx find(strcmp(S.ec_head.label, lbl{l}))];
    end
    
else
    if i == 1; 
        lbl = S.eo_head.label;
        eoidx = [];
        ecidx = [];
        for l = 1:length(lbl)
            ecidx = [ecidx find(strcmp(S.ec_head.label, lbl{l}))];
        end
        for l = ecidx
            eoidx = [eoidx find(strcmp(S.eo_head.label, S.ec_head.label{l}))];
        end
    else
        lbl = S.ec_head.label;
        ecidx = [];
        eoidx = [];
        for l = 1:length(lbl)
            eoidx = [eoidx find(strcmp(S.eo_head.label, lbl{l}))];
        end
        for l = eoidx
            ecidx = [ecidx find(strcmp(S.ec_head.label, S.eo_head.label{l}))];
        end
    end
end

tim_ax  = linspace(0,S.ec_head.nSamples/S.ec_head.Fs, S.ec_head.nSamples);
i       = 0;

% Label conditions
%--------------------------------------------------------------------------
clear cond ftdata 
if ~no_open
for o = 1:size(S.eo_data,3)
    i = i + 1;
    cond{i}         = 'Eyes Open';
    ftdata.trial{i} = squeeze(S.eo_data(eoidx,:,o));
    ftdata.time{i}  = tim_ax;
end
end

for c = 1:size(S.ec_data,3)
    i = i + 1;
 	cond{i} = 'Eyes Closed';
    ftdata.trial{i} = squeeze(S.ec_data(ecidx,:,c));
    ftdata.time{i}  = tim_ax;    
end

cd([Fanalysis fs 'MEEGs']);
ftdata.label    = lbl;
ftdata.label    = ftdata.label(ecidx);
pos             = S.ec_head.elec.chanpos;

D = spm_eeg_ft2spm(ftdata, [S.ec_head.name]);
D = type(D, 'single');

for c = 1:length(cond)
    D = conditions(D, c, cond{c});
end


% Locate EEG templates in SPM folder
%--------------------------------------------------------------------------
spm_path    = which('spm');
seppos      = find(spm_path == fs);
spm_path    = spm_path(1:seppos(end));
posfile     = [spm_path 'EEGtemplates\egi128_GSN_HydroCel.sfp'];

% Load sensor locations
%--------------------------------------------------------------------------
SP.task     = 'loadeegsens';
SP.source   = 'locfile';
SP.sensfile = posfile;
SP.D        = D;
D           = spm_eeg_prep(SP);

% Save MEEG object
%--------------------------------------------------------------------------
save(D);
cd(Fbase);
end