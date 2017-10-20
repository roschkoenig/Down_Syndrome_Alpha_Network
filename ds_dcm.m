function [DCM] = ds_dcm(sub, Fanalysis) 

% This function sets up the DCM for cross spectral densities to invert EEG
% alpha oscillation data from resting state EEG in people living with Down
% syndrome

% Housekeeping
%==========================================================================
clear DCM
fs          = filesep;
Fdata       = [Fanalysis fs 'MEEGs'];
Fdcm        = [Fanalysis fs 'DCM'];
spm('defaults', 'EEG');

% Set up DCM structure and invert baseline
%==========================================================================
DCM = [];

% Fix directory of canonical forward matrix
%--------------------------------------------------------------------------
DCM.xY.Dfile                            = [Fdata fs sub];
fixvol                                  = load(DCM.xY.Dfile);
if isfield(fixvol.D.other, 'inv')
    fullspmpath                             = which('spm');
    fixvol.D.other.inv{end}.forward.vol     = [fullspmpath(1:end-6) fs 'canonical' fs 'single_subj_T1_EEG_BEM.mat'];
    D = fixvol.D;
    save(DCM.xY.Dfile, 'D');
end

% Load MEEG object and extract sampling rate and info
%--------------------------------------------------------------------------
MEEG                = spm_eeg_load(DCM.xY.Dfile);
Fs                  = fsample(MEEG);
smpls               = size(MEEG,2);
timax               = linspace(0, smpls/Fs, smpls);

% Set up DCM details
%--------------------------------------------------------------------------
DCM.options.analysis = 'CSD';       % cross-spectral density 
DCM.options.model    = 'CMC';      	% structure canonical microcircuit (for now)
DCM.options.spatial  = 'IMG';           
DCM.options.Tdcm     = [1 40000];    % 1-30k ms 

DCM.options.Fdcm    = [1 30];     	% frequency range  
DCM.options.D       = 1;         	% frequency bin, 1 = no downsampling
DCM.options.Nmodes  = 8;          	% cosine reduction components used 
DCM.options.han     = 0;         	% no hanning 

DCM.options.location = 1;           % optmise location 
DCM.options.trials   = length(condlist(MEEG));    % index of ERPs within file

DCM.Sname           = {'lVI1', 'rVI1', 'lSPL', 'rSPL', 'lMFG', 'rMFG'};
DCM.Lpos            = [[-16; -92; 0] [12; -92; 21] [-48; -56; 52] [34; -51; 39] [-46; 37; 16] [46; 18; 21]];
DCM.M.Hz            = DCM.options.Fdcm(1):DCM.options.D:DCM.options.Fdcm(2);
DCM.xY.Hz           = DCM.M.Hz;
 
% Define different connection types
%==========================================================================
% Forward
%--------------------------------------------------------------------------
F      = zeros(6); 
F(3,1) = 1;     F(4,2) = 1;
F(5,3) = 1;     F(6,4) = 1;

% Backward
%--------------------------------------------------------------------------
B      = zeros(6);
B(1,3) = 1;     B(2,4) = 1;
B(3,5) = 1;     B(4,6) = 1;
   
% Lateral
%--------------------------------------------------------------------------
L      = zeros(6);
L(1,2) = 1;     L(2,1) = 1;
L(3,4) = 1;     L(4,3) = 1;
L(5,6) = 1;     L(6,5) = 1;
   
% Self
%--------------------------------------------------------------------------
S      = zeros(6);
S(1,1) = 1;     S(2,2) = 1;     S(3,3) = 1;
S(4,4) = 1;     S(5,5) = 1;     S(6,6) = 1;

 
% Define model arcitecture (A), conditional effects (B) and input (C) 
%==========================================================================
DCM.A{1}    =   F + L;
DCM.A{2}    =   B + L;
DCM.A{3}    =   S;

DCM.B           = {};
DCM.C           = sparse(length(DCM.A{1}),0); 

% Reorganise model parameters in specific structure
%==========================================================================
DCM.M.dipfit.Nm    = DCM.options.Nmodes;
DCM.M.dipfit.model = DCM.options.model;
DCM.M.dipfit.type  = DCM.options.spatial;

DCM.M.dipfit.Lpos  = DCM.Lpos;
DCM.M.dipfit.Nc    = size(MEEG,1);
DCM.M.dipfit.Ns    = length(DCM.A{1});
DCM.M.dipfit.silent_source = {};


% Define priors
%==========================================================================
% Load standard neural priors
%--------------------------------------------------------------------------
[pE,pC]  = ds_spm_dcm_neural_priors(DCM.A,DCM.B,DCM.C,DCM.options.model);

% EITHER    Manually define set (tested on YA008)
%--------------------------------------------------------------------------
% pE.T(1) = 1;
% pE.T(2) = 1;
% pE.T(3) = 3;
% pE.S(1) = -1;

% OR        Load set of empirical priors (from inv of YA008)
%--------------------------------------------------------------------------
load([Fdcm fs 'EmpPriors.mat']);
eE.G(:, 4:8)    = pE.G(:, 4:8);
pE              = eE;
pCvec           = spm_vec(pC);
pCvec           = pCvec * 2;
pC              = spm_unvec(pCvec,pC);

DCM.M.pE   = pE;
DCM.M.pC   = pC;

% Compute leadfields during first model inversion
%==========================================================================
if ~isfield(MEEG, 'inv')
    job{1}.spm.meeg.source.headmodel.D = {[Fdata fs sub '.mat']};
    job{1}.spm.meeg.source.headmodel.val = 1;
    job{1}.spm.meeg.source.headmodel.comment = '';
    job{1}.spm.meeg.source.headmodel.meshing.meshes.template = 1;
    job{1}.spm.meeg.source.headmodel.meshing.meshres = 3;
    job{1}.spm.meeg.source.headmodel.coregistration.coregspecify.fiducial(1).fidname = 'spmnas';
    job{1}.spm.meeg.source.headmodel.coregistration.coregspecify.fiducial(1).specification.select = 'nas';
    job{1}.spm.meeg.source.headmodel.coregistration.coregspecify.fiducial(2).fidname = 'spmlpa';
    job{1}.spm.meeg.source.headmodel.coregistration.coregspecify.fiducial(2).specification.select = 'lpa';
    job{1}.spm.meeg.source.headmodel.coregistration.coregspecify.fiducial(3).fidname = 'spmrpa';
    job{1}.spm.meeg.source.headmodel.coregistration.coregspecify.fiducial(3).specification.select = 'rpa';
    job{1}.spm.meeg.source.headmodel.coregistration.coregspecify.useheadshape = 0;
    job{1}.spm.meeg.source.headmodel.forward.eeg = 'EEG BEM';
    job{1}.spm.meeg.source.headmodel.forward.meg = 'Single Shell';
    spm_jobman('run', job);

    MEEG = spm_eeg_inv_forward([Fdata fs fname(MEEG)]);
end

% Invert DCM and log inversion output
%==========================================================================
DCM.name = [Fdcm fs 'DCM_' sub];
save(DCM.name, 'DCM');
