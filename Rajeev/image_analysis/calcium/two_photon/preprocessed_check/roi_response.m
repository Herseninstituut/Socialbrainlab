%% created by Rajeev Rajendran
% 2018_10_03
%% split quad data
pn = uigetdir();
cd(pn)
quad = load(uigetfile('*_quadrature.mat'));
quad = quad.quad_data;
quad_split1=NaN(1,floor((size(quad,2))/3));
quad_split2=NaN(1,floor((size(quad,2))/3));
quad_split3=NaN(1,floor((size(quad,2))/3));

for i=1:(floor(size(quad,2)))
    if floor((i-1)/3)==(i-1)/3
    quad_split1(1,floor(i/3)+1)= quad(1,i);
    elseif floor((i-2)/3)==(i-2)/3
    quad_split2(1,floor(i/3)+1)= quad(1,i);
    elseif floor((i-3)/3)==(i-3)/3
    quad_split3(1,floor(i/3)+1)= quad(1,i);
    end 
end
save('quad_split1','quad_split1')
save('quad_split2','quad_split2')
save('quad_split3','quad_split3')


%% load responses
% current folder have to be the correct one
% select the same split files for signals and running
SPSIG_file = load(uigetfile('*_SPSIG.mat'));
info = load(uigetfile('*normcorr.mat')); % 2p log file
r = load(uigetfile('quad_split*.mat')); % running log file
sig = SPSIG_file.sig;
den = SPSIG_file.den;
decon = SPSIG_file.decon;
sigraw = SPSIG_file.sigraw;

%% convert run data
r=r.quad_split1;
ra = 10;
dist = (2*pi*ra)*r/1000; %in cms
freqx = SPSIG_file.freq;
speed = (dist')*freqx; % in cms/s
%%  check stim timings
stims=info.info.frame;
stims=stims/3;
%% plot figures
% figure
% Spec=[3,22, 40, 60, 67, 70, 86, 89]; % specific ROIs
% for sr=1:size(Spec,2)
%     i=Spec(sr);
for i=1:size(sig,2)
    plot (sig(:,i),'color','r');
    hold on
    plot((1+(speed(:,1))/100),'color','g');
    legend ('activity df/F','speed 100 cm/s')
    title(['ROI ',num2str(i)])
    hold on
    for k=1:22
        y=(0:1:5);
        x=zeros(1,size(y,2));
        for j=1:size(y,2)
            x(j)=stims(k);
        end
        plot (x,y,'color','b','Linewidth',1)
        hold on
    end
    w = waitforbuttonpress;
    if isequal(w,0)
    end
    clf('reset')
    % plot (sig(:,i),'color','r');
    % plot (den(:,1),'color','g'), hold on
    % plot ((1+decon(:,1)),'color','b')
    % plot ((sigraw(:,1))/3200,'color','k')
end

% hold on
% Con = PP.Con;
% for k = 1:PP.Cnt
%     vx = Con(k).x;
%     vy = Con(k).y;
%     plot(vx, vy, 'r')
% end 
%% deconvolved data
for i=1:size(den,2)
    plot (sig(:,i),'color','r');
    hold on
    plot (den(:,i),'color','b');
    hold on
    plot (decon(:,i),'color','c');
    hold on
    plot((1+(speed(:,1))/100),'color','g');
    legend ('activity df/F','speed 100 cm/s')
    title(['ROI ',num2str(i)])
    hold on
    for k=1:22
        y=(0:1:5);
        x=zeros(1,size(y,2));
        for j=1:size(y,2)
            x(j)=stims(k);
        end
        plot (x,y,'color','b','Linewidth',1)
        hold on
    end
    w = waitforbuttonpress;
    if isequal(w,0)
    end
    clf('reset')
    % plot (sig(:,i),'color','r');
    % plot (den(:,1),'color','g'), hold on
    % plot ((1+decon(:,1)),'color','b')
    % plot ((sigraw(:,1))/3200,'color','k')
end
%% import stim
% have to be in the correct folder
stim = load(uigetfile('*.mat'));
contrast = stim.log(:,2);
angles = stim.log(:,1);
stim_duration = stim.Parameters.time;
ISI = stim.Parameters.wait;
stim_onset=round((info.info.frame)/4);
frame_rate=SPSIG_file.freq;
pre_stim_onset_frames=round(stim_onset-(frame_rate));
pre_stim_offset_frames=round(stim_onset-1);
stim_onset_frames = stim_onset;
stim_offset_frames = floor(stim_onset+(((stim_duration))*frame_rate));
a = unique(angles); % angles
c = unique(contrast); % contrast
x = size(angles,1)/((size(a,1))*(size(c,1)));% repeats
x=(1:x);
x=x';
rois=size(sig,2);
rois=(1:rois);
rois=rois';
%% stim details
for n=1:size(a) % unique angles
    for m=1:size(c) % unique contrasts
        k=1;
        for i=1:size(angles,1) % all stimuli
            if stim.log(i,1)==a(n) && stim.log(i,2)==c(m)
                stim.log(i,3)=k;
                k=k+1;
                %     else
                %        stim.log(i,3)=nan;
            end
        end
    end
end
%% create psth
resp_array = NaN(round((((stim_duration))+1)*frame_rate),size(angles,1),size(sig,2),2);
for i=1:size(sig,2)     %  i = rois
    for n=1:size(angles,1) % n = all stimuli
        k=1;
        for j=pre_stim_onset_frames(n):pre_stim_offset_frames(n)  % j = pre_stim_signals
            resp_array(k,n,i,1)=sig(j,i);
            k=k+1;
        end
        k=1;
        for j=stim_onset_frames(n):stim_offset_frames(n)  % j = pre_stim_signals
            resp_array(k,n,i,2)=sig(j,i);
            k=k+1;
        end
    end
end
%% responsive rois
[f,s,roi,b]=size(resp_array);
psth_array = NaN(f,size(x,1),size(a,1),size(c,1),roi,b);
% for i=1:roi
%     for j=1:size(c,1)
%         for m=1:size(a,1)
%             
% i=1:f; % frames
% j=1:size(x,1); % repeats
% k=1:size(a,1); % angles
% m=1:size(c,1); % contrasts
% n=1:roi;
% p=1:b; % prestim (1) or stim (2)
% resp_array=(i,all_stims,roi,p)
for p=1:b % prestim (1) or stim (2)
    for n=1:roi
        for m=1:size(c,1)% contrasts
            for k=1:size(a,1)% angles
                for j=1:size(x,1)% repeats
                    for q=1:s % all_stimuli
                        for i=1:f % frames
                            if stim.log(q,1)==a(k) && stim.log(q,2)==c(m) && stim.log(q,3)==x(j)
                               psth_array(i,j,k,m,n,p)=resp_array(i,q,n,p);
                               % psth_array(frame,repeat,angle,contrast,roi,prestim_or_stim)
                            end
                        end
                    end
                end
            end
        end
    end
end
save('quad_split1','quad_split1')
save('quad_split2','quad_split2')
save('quad_split3','quad_split3')
save('quad_split4','quad_split4')
save('psth_array','psth_array')
save('resp_array','resp_array')
save('psth.details.stim.mat','angles','stim','a','c','contrast','frame_rate',...
    'info','ISI','speed','dist','rois','s','f','b')


%% QC
% psth_array(i,j,k,m,n,p)=resp_array(i,q,n,p);
% psth_array(frame,repeat,angle,contrast,roi,prestim_or_stim)
% resp_array(frame,all_stim,roi,prestim_or_stim)

%make sure the value is not NaN

rng('shuffle')
x=s/(size(a,1)*size(c,1));
roi=size(rois,1);
i=randperm(f,1); % i=1:f; frames
n=randperm(roi,1);% n=1:roi;
p=randperm(b,1);% p=1:b; prestim (1) or stim (2)
q=randperm(s,1);% all stimuli
k=find(a==(stim.log(q,1)),1);% k=1:size(a,1); angles
m=find(c==(stim.log(q,2)),1);% m=1:size(c,1); contrasts
j=stim.log(q,3);% j=1:size(x,1); repeats
psth_array(i,j,k,m,n,p)==resp_array(i,q,n,p)
resp_array(i,q,n,p)
psth_array(i,j,k,m,n,p)


