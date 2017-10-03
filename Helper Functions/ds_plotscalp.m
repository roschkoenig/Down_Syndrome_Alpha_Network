function ds_plotscalp(volfile)

V   = spm_vol(volfile);
Y   = spm_read_vols(V);

Tmin = 3.66;
Tmap = Y * 0;
Tmap(Y > Tmin) = Y(Y > Tmin);

for r = 1:size(Tmap,1)
for c = 1:size(Tmap,2)
    Smap(r,c) = mean(find(Tmap(r,c,:)));
    if isnan(Smap(r,c)), Smap(r,c) = 0; end
    
    rl  = size(Tmap,1);     cl = size(Tmap,2);
    r2  = floor(rl/2);      c2 = floor(rl/2);
    
    if sqrt((r - r2)^2 + (c - c2)^2) > r2
        mask(r,c) = 0; 
    else mask(r,c) = 1;
    end
end
end

% Smap    = Smap .* mask;

imagesc(Smap); axis square; hold on
t = linspace(0,2*pi);
plot(r2*cos(t)+r2+.5,r2*sin(t)+r2+.5)