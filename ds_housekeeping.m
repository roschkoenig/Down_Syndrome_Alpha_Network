function D = ds_housekeeping
fs          = filesep;
Fbase     	= 'D:\Research_Data\1608 DS alpha\1708 Simplified'; 
Fscripts    = [Fbase fs 'Scripts'];
Fanalysis   = [Fbase fs 'Matlab Files'];
Fdata       = [Fanalysis fs 'MEEGs'];
Fdcm        = [Fanalysis fs 'DCM'];
addpath(genpath(Fscripts));

% Pack variables into export variable D
%--------------------------------------------------------------------------
D.Fbase     = Fbase;
D.Fscripts  = Fscripts;
D.Fanalysis = Fanalysis;
D.Fdata     = Fdata;
D.Fdcm      = Fdcm;
