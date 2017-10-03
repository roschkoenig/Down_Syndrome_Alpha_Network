clear all
% Housekeeping
%==========================================================================
D           = ds_housekeeping;
Fscripts 	= D.Fscripts;
Fbase       = D.Fbase;
Fanalysis   = D.Fanalysis;
Fdata       = D.Fdata;
fs          = filesep;
sub         = ds_definefiles(Fbase);

%% Project 3D sensor locations to 2D
%==========================================================================
omeegs = cellstr(spm_select('FPList', [Fanalysis fs 'MEEGs'], '^O.*.mat$'));
ymeegs = cellstr(spm_select('FPList', [Fanalysis fs 'MEEGs'], '^Y.*.mat$'));
meegs  = {omeegs{:}, ymeegs{:}}';

load(meegs{1});
chanxy 	= spm_eeg_project3D(D.sensors.eeg, 'EEG');
chanlbl = D.sensors.eeg.label;

for m = 1:length(meegs)
    load(meegs{m});
    for c = 1:length(D.channels)
        chid    = find(strcmp(D.channels(c).label, chanlbl));
        D.channels(c).X_plot2D = chanxy(1,chid);
        D.channels(c).Y_plot2D = chanxy(2,chid);
    end
    save(meegs{m}, 'D');
    clear D
end

%% Calculate time frequency analysis for each individual
%==========================================================================
omeegs = cellstr(spm_select('FPList', [Fanalysis fs 'MEEGs'], '^O.*.mat$'));
ymeegs = cellstr(spm_select('FPList', [Fanalysis fs 'MEEGs'], '^Y.*.mat$'));
meegs  = {omeegs{:}, ymeegs{:}}';
%%
for m = 1:length(meegs)
    clear job
    job{1}.spm.meeg.tf.tf.D = meegs(m);
    job{1}.spm.meeg.tf.tf.channels{1}.all = 'all';
    job{1}.spm.meeg.tf.tf.frequencies = [1:30];
    job{1}.spm.meeg.tf.tf.timewin = [-Inf Inf];
    job{1}.spm.meeg.tf.tf.method.mtmspec.timeres = 400;
    job{1}.spm.meeg.tf.tf.method.mtmspec.timestep = 50;
    job{1}.spm.meeg.tf.tf.method.mtmspec.bandwidth = 3;
    job{1}.spm.meeg.tf.tf.phase = 0;
    job{1}.spm.meeg.tf.tf.prefix = '';

    spm_jobman('run', job);
end

%% Average across repeated trials
%--------------------------------------------------------------------------
tfs = cellstr(spm_select('FPList', [Fanalysis fs 'MEEGs'], '^tf.*\.mat'));

for t = 1:length(tfs)
    clear job
    job{1}.spm.meeg.averaging.average.D = tfs(t);
    job{1}.spm.meeg.averaging.average.userobust.standard = false;
    job{1}.spm.meeg.averaging.average.plv = false;
    job{1}.spm.meeg.averaging.average.prefix = 'm';

    spm_jobman('run', job);
end

%% Generate scalp * frequency maps
%--------------------------------------------------------------------------
mtfs = cellstr(spm_select('FPList', [Fanalysis fs 'MEEGs'], '^mtf.*\.mat'));
for m = 1:length(mtfs)
    clear job
    job{1}.spm.meeg.images.convert2images.D = mtfs(m);
    job{1}.spm.meeg.images.convert2images.mode = 'scalp x frequency';
    job{1}.spm.meeg.images.convert2images.conditions = {'Eyes Open','Eyes Closed'}';
    job{1}.spm.meeg.images.convert2images.channels{1}.all = 'all';
    job{1}.spm.meeg.images.convert2images.timewin = [-Inf Inf];
    job{1}.spm.meeg.images.convert2images.freqwin = [4 13];
    job{1}.spm.meeg.images.convert2images.prefix = '';
    spm_jobman('run', job);
end

%% Smooth TF Data to allow averaging across individual differences
%--------------------------------------------------------------------------
dirs = cellstr(spm_select('FPList', [Fanalysis fs 'MEEGs'], 'dir')); 
filtr = '^condition_EyesClosed';

clear job
for d = 1:length(dirs)
    fname = spm_select('ExtFPList', dirs{d}, filtr);
    job{1}.spm.spatial.smooth.data = {fname};
    job{1}.spm.spatial.smooth.fwhm = [2 2 2];
    job{1}.spm.spatial.smooth.dtype = 0;
    job{1}.spm.spatial.smooth.im = 0;
    job{1}.spm.spatial.smooth.prefix = 's'; 
    spm_jobman('run', job);
end

%% Average first level scalp maps across individuals
%--------------------------------------------------------------------------
clear job flist fname
flist = cellstr(spm_select('FPList', Fdata, 'dir', '.'));

for f = 1:length(flist)
    fname{f} = [flist{f} fs 'scondition_EyesClosed.nii'];
    if f == 1,  cstring = '(i1';
    else        cstring = [cstring '+i' num2str(f)]; 
    end
end
fname   = fname';
cstring = [cstring ')/' num2str(length(flist))];

job{1}.spm.util.imcalc.input = fname;
job{1}.spm.util.imcalc.output = 'mean_sEC';
job{1}.spm.util.imcalc.outdir = {Fdata};
job{1}.spm.util.imcalc.expression = cstring;
job{1}.spm.util.imcalc.var = struct('name', {}, 'value', {});
job{1}.spm.util.imcalc.options.dmtx = 0;
job{1}.spm.util.imcalc.options.mask = 0;
job{1}.spm.util.imcalc.options.interp = 1;
job{1}.spm.util.imcalc.options.dtype = 4;

spm_jobman('run', job)

%% Second level - Kbit and age effects on theta-alpha networks
%==========================================================================
% Subject specific folders containing smoothed, mean images
%--------------------------------------------------------------------------
dirs = cellstr(spm_select('FPList', [Fanalysis fs 'MEEGs'], 'dir')); 
filtr = '^scondition_EyesClosed';

% Loop through folders to define filepaths and regressors
%--------------------------------------------------------------------------
count = 1;
fcnt  = 1;

clear job lvl1 K age pal vf split fullset

for d = 1:length(dirs)
    lvl1(d,:)   = spm_select('ExtFPList', dirs{d}, filtr);  
    name{d}     = dirs{d}(end-4:end);
    nah         = 0;
    
    try K(d)	= sub{d}.Kbit;      catch, nah = 1; end
    try age(d) 	= sub{d}.age;       catch, nah = 1; end
	try spl(d)  = sub{d}.split;     catch, nah = 1; end
    try sess(d) = sub{d}.sess;      catch, nah = 1; end
    
    if nah == 0
        fullset(fcnt)   = d;
        fcnt            = fcnt + 1;
    else 
        nah = 0;
    end
    
end

% Only include subjects will full datasets
%--------------------------------------------------------------------------
scans   = cellstr(lvl1(fullset,:));
K       = K(fullset)'; 
age     = age(fullset)';
spl     = spl(fullset)';
sess    = sess(fullset)';

clear job
Fspm = [Fanalysis fs 'SPM'];

job{1}.spm.stats.factorial_design.dir = {Fspm};
job{1}.spm.stats.factorial_design.des.mreg.scans    = scans;
job{1}.spm.stats.factorial_design.des.mreg.mcov     = struct('c', {}, 'cname', {}, 'iCFI', {}, 'iCC', {});
job{1}.spm.stats.factorial_design.des.mreg.incint   = 1;

job{1}.spm.stats.factorial_design.cov(1).c      = K;
job{1}.spm.stats.factorial_design.cov(1).cname  = 'K-bit';
job{1}.spm.stats.factorial_design.cov(1).iCFI   = 1;
job{1}.spm.stats.factorial_design.cov(1).iCC    = 1;

job{1}.spm.stats.factorial_design.cov(2).c      = age;
job{1}.spm.stats.factorial_design.cov(2).cname  = 'Age';
job{1}.spm.stats.factorial_design.cov(2).iCFI   = 1;
job{1}.spm.stats.factorial_design.cov(2).iCC    = 1;

bage    = age - min(age);
mage    = 2*bage / max(bage) - 1; 
bkbit   = K - min(K);
mkbit   = 2*bkbit / max(bkbit) - 1;
itx     = mkbit .* mage;

job{1}.spm.stats.factorial_design.cov(3).c      = itx;
job{1}.spm.stats.factorial_design.cov(3).cname  = 'Interaction';
job{1}.spm.stats.factorial_design.cov(3).iCFI   = 1;
job{1}.spm.stats.factorial_design.cov(3).iCC    = 1;

job{1}.spm.stats.factorial_design.cov(4).c      = spl;
job{1}.spm.stats.factorial_design.cov(4).cname  = 'Split';
job{1}.spm.stats.factorial_design.cov(4).iCFI   = 1;
job{1}.spm.stats.factorial_design.cov(4).iCC    = 1;

job{1}.spm.stats.factorial_design.cov(5).c      = sess;
job{1}.spm.stats.factorial_design.cov(5).cname  = 'Sess';
job{1}.spm.stats.factorial_design.cov(5).iCFI   = 1;
job{1}.spm.stats.factorial_design.cov(5).iCC    = 1;

job{1}.spm.stats.factorial_design.multi_cov = struct('files', {}, 'iCFI', {}, 'iCC', {});
job{1}.spm.stats.factorial_design.masking.tm.tm_none = 1;
job{1}.spm.stats.factorial_design.masking.im = 1;
job{1}.spm.stats.factorial_design.masking.em = {''};
job{1}.spm.stats.factorial_design.globalc.g_omit = 1;
job{1}.spm.stats.factorial_design.globalm.gmsca.gmsca_no = 1;
job{1}.spm.stats.factorial_design.globalm.glonorm = 1;

job{2}.spm.stats.fmri_est.spmmat = {[Fspm fs 'SPM.mat']};
job{2}.spm.stats.fmri_est.write_residuals = 0;
job{2}.spm.stats.fmri_est.method.Classical = 1;
spm_jobman('run', job);
clear job


% Define Contrasts
%--------------------------------------------------------------------------
job{1}.spm.stats.con.spmmat = {[Fspm fs 'SPM.mat']};

job{1}.spm.stats.con.consess{1}.tcon.name = 'Kbit';
job{1}.spm.stats.con.consess{1}.tcon.weights = [0 1 0 0];
job{1}.spm.stats.con.consess{1}.tcon.sessrep = 'none';

job{1}.spm.stats.con.consess{2}.tcon.name = 'Age';
job{1}.spm.stats.con.consess{2}.tcon.weights = [0 0 1 0];
job{1}.spm.stats.con.consess{2}.tcon.sessrep = 'none';

job{1}.spm.stats.con.consess{3}.tcon.name = 'Interaction';
job{1}.spm.stats.con.consess{3}.tcon.weights = [0 0 0 1];
job{1}.spm.stats.con.consess{3}.tcon.sessrep = 'none';

job{1}.spm.stats.con.delete = 1;
spm_jobman('run', job); 

d = job{1}.spm.stats.con.consess;
for dd = 1:length(d)
    subplot(1,length(d),dd)
    volfile = [Fspm fs 'spmT_000' num2str(dd) '.nii'];
    ds_plotscalp(volfile)
    title(d{dd}.tcon.name)
end

clear job;

%% Second level - PAL and VF correlates, post-hoc analysis
%==========================================================================
% Subject specific folders containing smoothed, mean images
%--------------------------------------------------------------------------
dirs = cellstr(spm_select('FPList', [Fanalysis fs 'MEEGs'], 'dir')); 
filtr = '^scondition_EyesClosed';

% Loop through folders to define filepaths and regressors
%--------------------------------------------------------------------------
count = 1;
fcnt  = 1;

clear job lvl1 K age pal vf split fullset

for d = 1:length(dirs)
    lvl1(d,:)   = spm_select('ExtFPList', dirs{d}, filtr);  
    name{d}     = dirs{d}(end-4:end);
    
    try pal(d)  = sub{d}.pal;       catch, nah = 1; end
    try vf(d) 	= sub{d}.verbalflu; catch, nah = 1; end
    try spl(d) 	= sub{d}.split;     catch, nah = 1; end
    try sess(d) = sub{d}.sess;      catch, nah = 1; end
    
    if nah == 0
        fullset(fcnt)   = d;
        fcnt            = fcnt + 1;
    else 
        nah = 0;
    end
    
end

% Only include subjects will full datasets
%--------------------------------------------------------------------------
scans   = cellstr(lvl1(fullset,:));
pal     = pal(fullset)';
vf      = vf(fullset)';
spl     = spl(fullset)';
sess    = sess(fullset)';

clear job
Fspm = [Fanalysis fs 'SPM'];
p   = sub;

job{1}.spm.stats.factorial_design.dir = {Fspm};
job{1}.spm.stats.factorial_design.des.mreg.scans    = scans;
job{1}.spm.stats.factorial_design.des.mreg.mcov     = struct('c', {}, 'cname', {}, 'iCFI', {}, 'iCC', {});
job{1}.spm.stats.factorial_design.des.mreg.incint   = 1;

% Define regressors
%--------------------------------------------------------------------------
job{1}.spm.stats.factorial_design.cov(1).c      = pal;
job{1}.spm.stats.factorial_design.cov(1).cname  = 'PAL';
job{1}.spm.stats.factorial_design.cov(1).iCFI   = 1;
job{1}.spm.stats.factorial_design.cov(1).iCC    = 1;

job{1}.spm.stats.factorial_design.cov(2).c      = vf;
job{1}.spm.stats.factorial_design.cov(2).cname  = 'Verbal Fluency';
job{1}.spm.stats.factorial_design.cov(2).iCFI   = 1;
job{1}.spm.stats.factorial_design.cov(2).iCC    = 1;

job{1}.spm.stats.factorial_design.cov(3).c      = spl;
job{1}.spm.stats.factorial_design.cov(3).cname  = 'Split';
job{1}.spm.stats.factorial_design.cov(3).iCFI   = 1;
job{1}.spm.stats.factorial_design.cov(3).iCC    = 1;

job{1}.spm.stats.factorial_design.cov(4).c      = sess;
job{1}.spm.stats.factorial_design.cov(4).cname  = 'Sess';
job{1}.spm.stats.factorial_design.cov(4).iCFI   = 1;
job{1}.spm.stats.factorial_design.cov(4).iCC    = 1;

% Set SPM design parameters
%--------------------------------------------------------------------------
job{1}.spm.stats.factorial_design.multi_cov = struct('files', {}, 'iCFI', {}, 'iCC', {});
job{1}.spm.stats.factorial_design.masking.tm.tm_none = 1;
job{1}.spm.stats.factorial_design.masking.im = 1;
job{1}.spm.stats.factorial_design.masking.em = {''};
job{1}.spm.stats.factorial_design.globalc.g_omit = 1;
job{1}.spm.stats.factorial_design.globalm.gmsca.gmsca_no = 1;
job{1}.spm.stats.factorial_design.globalm.glonorm = 1;

% Estimate SPM
%--------------------------------------------------------------------------
job{2}.spm.stats.fmri_est.spmmat = {[Fspm fs 'SPM.mat']};
job{2}.spm.stats.fmri_est.write_residuals = 0;
job{2}.spm.stats.fmri_est.method.Classical = 1;

spm_jobman('run', job); 
clear job

% Define Contrasts
%--------------------------------------------------------------------------
job{1}.spm.stats.con.spmmat = {[Fspm fs 'SPM.mat']};

job{1}.spm.stats.con.consess{1}.tcon.name = 'PAL';
job{1}.spm.stats.con.consess{1}.tcon.weights = [0 1 0];
job{1}.spm.stats.con.consess{1}.tcon.sessrep = 'none';

job{1}.spm.stats.con.consess{2}.tcon.name = 'VF';
job{1}.spm.stats.con.consess{2}.tcon.weights = [0 0 1];
job{1}.spm.stats.con.consess{2}.tcon.sessrep = 'none';

job{1}.spm.stats.con.delete = 1;
spm_jobman('run', job);

d = job{1}.spm.stats.con.consess;
for dd = 1:length(d)
    subplot(1,length(d),dd)
    volfile = [Fspm fs 'spmT_000' num2str(dd) '.nii'];
    ds_plotscalp(volfile)
    title(d{dd}.tcon.name)
end

clear job


