Cbase = ['/data/rrosch/ds_cluster'];
Cdata = ['/data/rrosch/ds_cluster'];
Fdata = ['D:\Research_Data\1608 DS alpha\1708 Simplified\Cluster\Data'];
%%
old     = cellstr(spm_select('FPList', Fdata, '^O.*.mat$'));
young   = cellstr(spm_select('FPList', Fdata, '^Y.*.mat$'));
datafiles = {old{:}, young{:}}';

for d = 1:length(datafiles)
    datafile    = datafiles{d};
    load(datafile);
    D.path      = Cdata;
    save([Fdata fs D.fname], 'D');
    clear D
end

%%
dcms = cellstr(spm_select('FPList', Fdata, '^D'));
for d = 1:length(dcms)
    dcm = dcms{d};
    load(dcm);
    DCM.xY.Dfile = [Cdata '/' DCM.xY.Dfile(end-4:end) '.mat'];
    DCM.name     = [Cdata '/' DCM.name(end-8:end)]
    save([Fdata fs DCM.name(end-8:end)], 'DCM');
end