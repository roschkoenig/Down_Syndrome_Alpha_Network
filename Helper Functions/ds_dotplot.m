function ds_dotplot(M)
%==========================================================================
% This function plots a dot-plot with mean and standard errors for grouped
% variables (here applied for the K-bit separated alpha bandpowers
%==========================================================================

% Loop through individual groups
%==========================================================================
names = [];
for m = 1:length(M)
  
    % Define slightly jittered x-locations 
    %----------------------------------------------------------------------
    o       = ones(1, length(M{m}.d)) * (m-1);
    jit     = randn(1, length(o))/15;
    
    % Plot individual data points (here mean power per subject)
    %----------------------------------------------------------------------
    scatter(o + jit, M{m}.d, 'filled'); hold on
    mm = median(M{m}.d);
    plot([m-1.25 m-0.75], [mm mm], 'k', 'Linewidth', 2);
    names = [names; M{m}.name];
    
    % Calculate and plot means, and standard errors
    %----------------------------------------------------------------------
% 	s = std(mean(M{m}.d,2));
%     sem = s / sqrt(length(mean(M{m}.d,2)));
%     u = mm + sem;
%     l = mm - sem;
%     plot([m-1 m-1], [l u], 'r', 'Linewidth', 1.5);    
%     
%     
%     
end

set(gca, 'XTick', [1:length(M)]-1);
set(gca, 'XTickLabel', names);