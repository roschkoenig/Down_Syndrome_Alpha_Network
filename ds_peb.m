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
fs          = filesep;

% Load all files into single cell array
%--------------------------------------------------------------------------
for pp = 1:length(p)
    TCM = load([Fdcm fs 'DCM_' p{pp}.name '.mat']);
    ACM{pp} = TCM.DCM; 
end

% PEB Analysis
%==========================================================================
clear X Kbit Age 

% Effects of age and kbit and session number
%--------------------------------------------------------------------------
for aa = 1:length(ACM)  
    dotpos      = find(ACM{aa}.name == '.');
    seppos      = find(ACM{aa}.name == '/');
    if isempty(seppos), seppos = find(ACM{aa}.name == '\'); end
    if isempty(dotpos), dotpos = length(ACM{aa}.name) + 1; end
    
    astring{aa} = ACM{aa}.name(seppos(end)+5 : dotpos(end)-1);
end

for pp = 1:length(p)
    pstring = p{pp}.name;
    aindx   = find(strcmp(pstring, astring));
    kbit(aindx) = p{pp}.Kbit;
    age(aindx)  = p{pp}.age;
    sess(aindx) = p{pp}.sess;
    split(aindx) = p{pp}.split;
end

% Mean-center and normalise kbit scores
%--------------------------------------------------------------------------
kbit = 2*(kbit - mean(kbit)) / range(age);

% Orthogonalise Age
%--------------------------------------------------------------------------
age         = (age - mean(age)) / range(age);

lin_coeff   = polyfit(kbit, age, 1);
agefit      = polyval(lin_coeff, kbit);
age         = 2*(age-agefit);

clear X
X(:,1) = kbit;      Xnames{1}   = ('Kbit');
X(:,2) = age;       Xnames{2}   = ('Age');
X(:,3) = sess;      Xnames{3}   = ('Session');
X(:,4) = split;     Xnames{4}   = ('Split');

% Main Group Effect
%--------------------------------------------------------------------------
X(:,end+1)  = ones(1,length(ACM));
Xnames{end+1}   = ('Group Mean');

clear C F PEB

C{1} = {'A{1}(3,1)', 'A{1}(4,2)', 'A{1}(5,3)', 'A{1}(6,4)', 'A{2}(3,1)', 'A{2}(4,2)', 'A{2}(5,3)', 'A{2}(6,4)'}; % Forward
C{2} = {'A{3}(1,3)', 'A{3}(2,4)', 'A{3}(3,5)', 'A{3}(4,6)', 'A{4}(1,3)', 'A{4}(2,4)', 'A{4}(3,5)', 'A{4}(4,6)'}; % Backward
C{3} = {'A'};   % All extrinsic 

C{4} = {'G(1,1)', 'G(1,2)', 'G(1,3)', 'G(1,4)', 'G(1,5)', 'G(1,6)'};      % Modulatory
C{5} = {C{1}{:}, C{4}{:}};  % Forward Mod
C{6} = {C{2}{:}, C{4}{:}};  % Backward Mod
C{7} = {C{3}{:}, C{4}{:}};  % Extrinsic and Mod

Cnames = {'F', 'B', 'FB', 'i', 'Fi', 'Bi', 'FBi'};

% Run PEB
%--------------------------------------------------------------------------
M.X         = X;
M.Xnames    = Xnames;
for c = 1:length(C)
    PEB{c} = spm_dcm_peb(ACM', M, C{c});
end

%% Plot Bayesian model comparison at the second level
%--------------------------------------------------------------------------
for b = 1:length(PEB)
   F(b) = PEB{b}.F; 
end
Fsorted     = sort(F);
dF          = Fsorted(end) - Fsorted(end-1);
[val iF]   = max(F);

figure
subplot(2,1,1)
    bar(F - min(F));
    title('Bayesian Model Comparison')
    xlabel(['dF = ' num2str(dF)]);
    ylabel('Free Energy');
    set(gca, 'XTickLabel', Cnames);
    hold on
    
    bar(iF, F(iF) - min(F), 'r');
    
subplot(2,1,2)
    bar(spm_softmax(F'));
%% Perform Bayesian model average at second level and extract reduced models
%==========================================================================
clear M
M.X         = X;
[REB RCM]   = spm_dcm_peb(ACM', M, C{iF});
BMA         = spm_dcm_peb_bmc(PEB{iF});
save([Fdcm fs 'Reduced_Models'], 'RCM');
%%
clear g
p = ds_definefiles(Fbase);
for r = 1:length(RCM)
    g(r,:)      = RCM{r}.Ep.G(:,1);
    kbit(r)     = p{r}.Kbit;
    age(r)      = p{r}.age;
end

% scatter(g(:,5), kbit, 'filled'); hold on
% scatter(g(:,6), kbit, 'filled')
subplot(2,1,1)
scatter(kbit, mean(g(:,1:2),2), [], 'filled'); hold on
subplot(2,1,2)
scatter(kbit, mean(g(:,5:6),2), [], 'filled');
%% Plotting parameter changes per region
%==========================================================================
% Extract parameter estimates and covariances from reduced model
%--------------------------------------------------------------------------
Cp = diag(BMA.Cp);
Ep = BMA.Ep;

% Plot parameter changes for each region
%--------------------------------------------------------------------------
for r = 1:6
rid = [];

% Identify indices from specified region from PEB
%--------------------------------------------------------------------------
for p = 1:length(PEB{iF}.Pnames)
    if ~isempty(findstr(PEB{iF}.Pnames{p}, ['G(' num2str(r) ',']));
        rid = [rid, p];
    end
end

R{r} = rid;
end

% Define order for plotting of parameters
%--------------------------------------------------------------------------
pid = {[1 2 5], [3 4 6]};   % mod(sp, ii, ss); exc(ss>ii, dp>ii, ss>sp); 

for e = 1:2
for r = 1:length(R)

Np  = length(PEB{iF}.Pnames);
mod = R{r}(pid{1});
exc = R{r}(pid{2});
id = [mod + (e-1)*Np, exc + (e-1)*Np];

figure(e)
subplot(3,2,r)
spm_plot_ci(Ep(id), Cp(id));
title([BMA.Xnames{e} ' effect in ' ACM{1}.Sname{r}]);
set(gca, 'XTickLabel', {BMA.Pnames{mod}, BMA.Pnames{exc}});

end
end

for e = 1:2
for sp = 1:3
    
figure(e)
subplot(3,2, 2*sp-1),   yl(1,:) = ylim;
subplot(3,2, 2*sp),     yl(2,:) = ylim;

subplot(3,2, 2*sp-1),   ylim([min(yl(:,1)), max(yl(:,2))]);
subplot(3,2, 2*sp),     ylim([min(yl(:,1)), max(yl(:,2))]);

end
end


%% This is where code goes to die
%--------------------------------------------------------------------------
% Define PEBs for 2nd level model comparison
%--------------------------------------------------------------------------
% Generate field names for intrinsic connections
%--------------------------------------------------------------------------
% Ns      = size(ACM{1}.A{1},1);   % Number of sources in the model
% types   = {'mod', 'inh', 'exc'};
% tid     = {[7 4 1], [2 3], [5 6 8]};
% 
% for t = 1:length(types);
% count   = 0;  
% for s = 1:Ns;
% for ti = 1:length(tid{t});
%     count           = count + 1;
%     i               = tid{t}(ti);
%     Gtyp{t}{count} = ['G(' num2str(s) ',' num2str(i) ')'];
% end
% end
% end


% for g = 1:length(Gtyp)
%     C{end+1} = Gtyp{g};
% end
% C{end+1} = {'M', Gtyp{3}{:}};
% C{end+1} = {Gtyp{1}{:}, Gtyp{3}{:}};
% C{end+1} = {Gtyp{2}{:}, Gtyp{3}{:}};
% 
% Cnames = {'Ext', 'F', 'B', 'M', 'T', 'Mod', 'Inh', 'Exc', 'M/Exc', 'Mod/Exc', 'Exc/Inh', 'all'};