function [S no_open] = ds_read(sub)
%--------------------------------------------------------------------------
% This function reads the data file from a subject specifier and returns
% either the header, or header and data
% [hdr, [data]] = ee_read(subject, [window length in minutes]);
%--------------------------------------------------------------------------

% Define optional variables
%--------------------------------------------------------------------------
% switch nargin 
% end
fs      = filesep;
Fbase   = 'D:\Research_Data\1608 DS alpha\Data';
FEO     = [Fbase fs 'EO processed 2-sec'];
FEC     = [Fbase fs 'EC processed 2-sec'];

% Load files 
%--------------------------------------------------------------------------
disp('closed')
D.ec_data         = ft_read_data(sub.ecfile);
D.ec_head         = ft_read_header(sub.ecfile);
D.ec_head.name    = sub.name;
D.ec_head.orig.filepath = FEC; 

if D.ec_head.Fs > 250
	cfg.resamplefs  = 250;
    cfg.detrend     = 'no';
    D.ec_data = ft_resampledata(cfg, D.ec_data);
end

try                     % problems with subject OA026 EO condition
disp('open')
D.eo_data         = ft_read_data(sub.eofile);
D.eo_head         = ft_read_header(sub.eofile);
D.eo_head.name          = sub.name;
D.eo_head.orig.filepath = FEO; 

if D.eo_head.Fs > 250
    cfg.resamplefs  = 250;
    cfg.detrend     = 'no';
    D.eo_data = ft_resampledata(cfg, D.eo_data);
end
no_open = 0;
catch
    no_open = 1;
end

S = D;