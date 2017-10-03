function ee_eegplot(data, head, sub, varargin)
%--------------------------------------------------------------------------
% This function plots a ten second EEG section as average montage.
% ee_eegplot(Data, Header, Subject Structure, [Starting time in seconds])
%--------------------------------------------------------------------------

if nargin > 3, eeg_start = varargin{1}
else eeg_start = 10; end

figure;

% Plot EEG trace across all channels
%--------------------------------------------------------------------------
for c = 1:size(data,1)
    if any(c==head.left), cl = [0.7 0 0]; end
    if any(c==head.centre), cl = [0 0 0]; end
    if any(c==head.right), cl = [0.4 0 0.7]; end
    timeaxis = linspace(0,10,10*head.Fs);
    subplot(4,1,1:3) 
        plot(timeaxis, data(c, eeg_start*head.Fs + (1:10*head.Fs)) - 100*c, 'Color', cl); hold on 

end
title([sub.cond '- aged: ' num2str(sub.age) ' months. Filter: ' head.filter]);
eeg = gca;
eeg.YTick = flip([-100*((1:size(data,1)))]);
eeg.YTickLabel = flip(head.label);
ylim([-Inf 100]);
xlabel('time in seconds');
set(gcf, 'Color', 'white');

% Plot single trace across whole time window to check for artefact areas
%--------------------------------------------------------------------------
subplot(4,1,4)
    timeaxis = linspace(0,length(data)/(head.Fs*60), length(data));
	section = zeros(1,length(data));
    section(1,floor(eeg_start*head.Fs):ceil((eeg_start+10)*head.Fs)) = ones*1000;
    section = section;
    area(timeaxis, section, 'FaceColor', [0.8, 0.8, 0.8], 'LineStyle', ':'), hold on;
    plot(timeaxis, data(head.centre(1),:)+500); hold on;
    xlabel('time in minutes');
    xlim([0 Inf]);
    ylim([min(data(head.centre(1),:))+500, max(data(head.centre(1),:))+500]);
end
