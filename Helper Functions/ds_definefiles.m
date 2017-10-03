function [p] = ds_definefiles(Fbase)
%%=========================================================================
% This function defines the locations of the EDF files, as well as the
% start and duration of the relevant artefact-free time windows.
%
% p = patients (p{x}.name, p{x}.start, p{x}.length, p{x}.path)
% c = controls (c{x}.name, c{x}.start, c{x}.length, c{x}.path)

fs      = filesep;
FEO     = [Fbase fs 'Data' fs 'EO processed 2-sec'];
FEC     = [Fbase fs 'Data' fs 'EC processed 2-sec'];   

yeo = cellstr(spm_select('FPList', FEO, '^Y.*\.set$'));
oeo = cellstr(spm_select('FPList', FEO, '^O.*\.set$'));
yec = cellstr(spm_select('FPList', FEC, '^Y.*\.set$'));
oec = cellstr(spm_select('FPList', FEC, '^O.*\.set$'));

count = 1;
for y = 1:length(yeo)
    p{count}.name = yeo{y}(end-18:end-14);
    p{count}.eofile = yeo{y};
    p{count}.ecfile = yec{y}; 
    count = count+1;
end

for o = 1:length(oeo)
    p{count}.name = oeo{o}(end-18:end-14);
    p{count}.eofile = oeo{o};
    p{count}.ecfile = oec{o}; 
    count = count+1;    
end

for i = 1:length(p)
    if p{i}.name == 'A058 ', p{i}.name = 'OA058'; end;
    if p{i}.name == 'A057 ', p{i}.name = 'YA057'; end;
    switch p{i}.name
        case 'OA007', p{i}.Kbit = 74; p{i}.age = 44; p{i}.sess = -1; p{i}.verbalflu = 11; p{i}.pal = 14; p{i}.split = -1;
        case 'OA010', p{i}.Kbit = 78; p{i}.age = 41; p{i}.sess = 1; p{i}.verbalflu = 27; p{i}.pal = 14; p{i}.split = -1;
        case 'OA019', p{i}.Kbit = 12; p{i}.age = 56; p{i}.sess = -1; p{i}.verbalflu = 11; p{i}.pal = 4; p{i}.split = -1;
        case 'OA027', p{i}.Kbit = 40; p{i}.age = 37; p{i}.sess = -1; p{i}.verbalflu = 5; p{i}.pal = 14; p{i}.split = -1;
        case 'OA058', p{i}.Kbit = 63; p{i}.age = 49; p{i}.sess = 1; p{i}.verbalflu = 11; p{i}.pal = 14; p{i}.split = 1;  
        case 'OA060', p{i}.Kbit = 73; p{i}.age = 38; p{i}.sess = 1; p{i}.verbalflu = 14; p{i}.pal = 11; p{i}.split = 1;
        case 'OA075', p{i}.Kbit = 63; p{i}.age = 51; p{i}.sess = -1; p{i}.verbalflu = 13; p{i}.pal = 10; p{i}.split = 1;    
        case 'OA103', p{i}.Kbit = 66; p{i}.age = 46; p{i}.sess = -1; p{i}.verbalflu = 20; p{i}.pal = 12; p{i}.split = 1;
        case 'OA138', p{i}.Kbit = 56; p{i}.age = 50; p{i}.sess = 1; p{i}.verbalflu = 18; p{i}.pal = 6; p{i}.split = 1;    
        case 'OA195', p{i}.Kbit = 55; p{i}.age = 43; p{i}.sess = -1; p{i}.split = 1;
        case 'YA001', p{i}.Kbit = 83; p{i}.age = 31; p{i}.sess = -1; p{i}.verbalflu = 16; p{i}.pal = 15; p{i}.split = -1;    
        case 'YA003', p{i}.Kbit = 66; p{i}.age = 26; p{i}.sess = -1; p{i}.verbalflu = 13; p{i}.pal = 16; p{i}.split = -1;
        case 'YA008', p{i}.Kbit = 58; p{i}.age = 32; p{i}.sess = 1; p{i}.verbalflu = 15; p{i}.pal = 10; p{i}.split = -1;
        case 'YA009', p{i}.Kbit = 58; p{i}.age = 17; p{i}.sess = -1; p{i}.verbalflu = 6; p{i}.pal = 16; p{i}.split = -1;   
        case 'YA011', p{i}.Kbit = 38; p{i}.age = 35; p{i}.sess = -1; p{i}.verbalflu = 13; p{i}.pal = 13; p{i}.split = -1;       
        case 'YA014', p{i}.Kbit = 51; p{i}.age = 27; p{i}.sess = 1; p{i}.verbalflu = 12; p{i}.pal = 7; p{i}.split = -1;
        case 'YA019', p{i}.Kbit = 53; p{i}.age = 26; p{i}.sess = -1; p{i}.verbalflu = 10; p{i}.pal = 14; p{i}.split = 1;
        case 'YA021', p{i}.Kbit = 38; p{i}.age = 33; p{i}.sess = -1; p{i}.verbalflu = 13; p{i}.pal = 2; p{i}.split = -1;    
        case 'YA022', p{i}.Kbit = 86; p{i}.age = 33; p{i}.sess = 1; p{i}.verbalflu = 12; p{i}.pal = 14; p{i}.split = -1;
        case 'YA023', p{i}.Kbit = 10; p{i}.age = 23; p{i}.sess = -1; p{i}.split = -1;
        case 'YA025', p{i}.Kbit = 51; p{i}.age = 34; p{i}.sess = -1; p{i}.verbalflu = 14; p{i}.pal = 11; p{i}.split = -1;
        case 'YA028', p{i}.Kbit = 61; p{i}.age = 25; p{i}.sess = 1; p{i}.verbalflu = 12; p{i}.pal = 13; p{i}.split = -1; 
        case 'YA032', p{i}.Kbit = 33; p{i}.age = 22; p{i}.sess = 1; p{i}.verbalflu = 9; p{i}.pal = 10; p{i}.split = -1;    
        case 'YA038', p{i}.Kbit = 44; p{i}.age = 17; p{i}.sess = 1; p{i}.verbalflu = 8; p{i}.pal = 13; p{i}.split = -1;
        case 'YA050', p{i}.Kbit = 56; p{i}.age = 16; p{i}.sess = 1; p{i}.verbalflu = 8; p{i}.pal = 9; p{i}.split = 1;    
        case 'YA051', p{i}.Kbit = 27; p{i}.age = 31; p{i}.sess = -1; p{i}.verbalflu = 12; p{i}.pal = 4; p{i}.split = 1; 
        case 'YA057', p{i}.Kbit = 32; p{i}.age = 21; p{i}.sess = 1; p{i}.split = 1;
        case 'YA066', p{i}.Kbit =102; p{i}.age = 27; p{i}.sess = 1; p{i}.verbalflu = 24; p{i}.pal = 19; p{i}.split = 1;      
        case 'YA067', p{i}.Kbit = 64; p{i}.age = 18; p{i}.sess = -1; p{i}.verbalflu = 8; p{i}.pal = 13; p{i}.split = 1;   
        case 'YA072', p{i}.Kbit = 58; p{i}.age = 26; p{i}.sess = 1; p{i}.verbalflu = 14; p{i}.pal = 10; p{i}.split = 1;
        case 'YA073', p{i}.Kbit = 54; p{i}.age = 35; p{i}.sess = 1; p{i}.verbalflu = 17; p{i}.pal = 9; p{i}.split = 1;  
        case 'YA088', p{i}.Kbit = 27; p{i}.age = 19; p{i}.sess = 1; p{i}.verbalflu = 6; p{i}.pal = 9; p{i}.split = 1; 
        case 'YA090', p{i}.Kbit = 70; p{i}.age = 20; p{i}.sess = 1; p{i}.verbalflu = 9; p{i}.pal = 17; p{i}.split = 1; 
        case 'YA091', p{i}.Kbit = 53; p{i}.age = 20; p{i}.sess = 1; p{i}.verbalflu = 14; p{i}.pal = 10; p{i}.split = 1; 
        case 'YA096', p{i}.Kbit = 68; p{i}.age = 27; p{i}.sess = 1; p{i}.verbalflu = 22; p{i}.split = 1;      
        case 'YA117', p{i}.Kbit = 53; p{i}.age = 17; p{i}.sess = 1; p{i}.verbalflu = 16; p{i}.pal = 13; p{i}.split = 1;
    end
    
end