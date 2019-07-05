%% StimulusPreparationMice

%% Visualization of the raw recording %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

data = data((600*fs)+1:end);

%%

trial = 1;

tmpData = data(((trial-1)*10*fs)+1:trial*10*fs);

figure('Position',get(groot,'Screensize'))
subplot(2,1,1)
window = hamming(512);
noverlap = 384;
nfft = 1024;
[~,F,T,P] = spectrogram(tmpData,window,noverlap,nfft,fs,'yaxis');
imagesc(T,F,10*log10(P));
ylim([0000 100000])
set(gca,'YTick',[0 10000 20000 30000 40000 50000 60000 70000 80000 90000 100000])
set(gca,'YTickLabels',{'0','10','20','30','40','50','60','70','80','90','100'})
ylabel('Fq (kHz)')
set(gca,'XTick',[0 1 2 3 4 5 6 7 8 9 10 11 12])
% set(gca,'XTickLabel',{'-4','-3','-2','-1','0','1','2','3','4','5','6','7','8'})
xlabel('Time (sec)')
colormap(gray)
colormap(flipud(colormap))
set(gca,'clim',[-105 -80]);
grid off
box on
set(gca,'YDir','Normal')

subplot(2,1,2)
plot(tmpData)

%% Save the separated recording

audiowrite('_Squeak1.wav',tmpData,fs)

%% Trimming %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure('Position',get(groot,'Screensize'))
subplot(2,1,1)
window = hamming(512);
noverlap = 384;
nfft = 1024;
[~,F,T,P] = spectrogram(data,window,noverlap,nfft,fs,'yaxis');
imagesc(T,F,10*log10(P));
ylim([0000 100000])
set(gca,'YTick',[0 10000 20000 30000 40000 50000 60000 70000 80000 90000 100000])
set(gca,'YTickLabels',{'0','10','20','30','40','50','60','70','80','90','100'})
ylabel('Fq (kHz)')
% set(gca,'XTick',[0 1 2 3 4 5 6 7 8 9 10 11 12])
% set(gca,'XTickLabel',{'-4','-3','-2','-1','0','1','2','3','4','5','6','7','8'})
xlabel('Time (sec)')
colormap(gray)
colormap(flipud(colormap))
set(gca,'clim',[-105 -80]);
grid off
box on
set(gca,'YDir','Normal')

subplot(2,1,2)
plot(data)
% plot((1:length(data))*(1/fs),data)
% xlabel('Time (sec)')

%% Get the onset sample

[x1,~] = ginput(1);

%% Get the offset sample

[x2,~] = ginput(1);

%% Cut the relevant section

periStimDur = 0.005; % in sec

data2 = data(round(x1-(fs*periStimDur)):round(x2+(fs*periStimDur)));

%% Save the raw stimulus

audiowrite('_Squeak1_Raw.wav',data2,fs)

%% Tapering %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

taperLength = 0.01; % in sec, the summation of beggining and the ending

w = tukeywin(length(data2),(taperLength*fs)/length(data2));
data3 = data2.*w;

%% Save the tapered stimulus

audiowrite('_Squeak1_Tapered.wav',data3,fs)
