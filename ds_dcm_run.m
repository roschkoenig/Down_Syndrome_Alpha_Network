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

for i = 1:length(p)
    sub     = p{i}.name;
    DCM{i}  = ds_dcm(sub, Fanalysis);
end