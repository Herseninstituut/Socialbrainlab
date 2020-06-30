%% created by Rajeev Rajendran
% 2019_09_24
% Finds responsive cells based on the given stimuli, calculated by Wilcoxon
% sign rank test between paired trials of averaged signals between pre-
% onset and post- onset time windows. Correct stim indices and appropriate
% stim duration details are needed as input arguments.


function responsive_cells_sigraw(stim_first,stim_last,pre_stim_dur,stim_dur,...
    post_stim_dur)

% select the correct directory and make it current
[fn , pn] = uigetfile('*_SPSIG.mat');
cd(pn)

%% load data files
SPSIG_file = load(uigetfile('*_SPSIG.mat')); % signal file
info = load(uigetfile('*_normcorr.mat', 'info')); % normcorr mat file
quad = load(uigetfile('*quad_split*.mat')); % run file
stim=floor((info.info.frame)/3); % get stim onset and offset frames

%% calculate speed
run=['quad.quad_',pn(end-6:end-1)];
run=eval(run);
ra = 10;
dist = (2*pi*ra)*run/1000; %in cms
speed = (dist)*SPSIG_file.freq; % in cms/s

n=size(SPSIG_file.sig,1); % total frames available

%% select stim on and off frames
Ts=size((stim_first:2:stim_last),2);
for i=1:Ts
    StimOn(i,1)=stim((stim_first+(2*i-2)),1);
    StimOff(i,1)=stim((stim_first+(2*i-2))+1,1);
end

%% distribute the signals into pre and stim groups for all stims
% no. of frames before stim onset
Tf_pre=ceil((pre_stim_dur)*SPSIG_file.freq);
% no. of frames from stim onset
Tf_st=ceil((stim_dur+post_stim_dur)*SPSIG_file.freq);
% no. of ROIs
Tr=size(SPSIG_file.sig,2);
% make blank array for stats
stat_arr_raw=nan(Tf_st,2,Ts,Tr);
for i=1:Tr % for each ROI
    for j=1:Ts % for each stim
        for k=1:Tf_pre % for each frame
            if (StimOn(j,1)+k-(Tf_pre+1))>=1
                stat_arr_raw(k,1,j,i)=SPSIG_file.sigraw((StimOn(j,1)+k-(Tf_pre+1)),i); % baseline
            end
        end
        for k=1:Tf_st % for each frame
            stat_arr_raw(k,2,j,i)=SPSIG_file.sigraw((StimOn(j,1)+k-1+2),i); % post stimon
        end
    end
end
%%
stat_arr_dF=nan(Tf_st,2,Ts,Tr);
for i=1:Tr % for each ROI
    for j=1:Ts % for each stim
        for k=1:Tf_pre % for each frame
            if nanmean(stat_arr_raw(:,1,j,i))<0
                stat_arr_dF(k,1,j,i)=(stat_arr_raw(k,1,j,i)-nanmean(stat_arr_raw(:,1,j,i)))/((-1)*nanmean(stat_arr_raw(:,1,j,i)));
            else
                stat_arr_dF(k,1,j,i)=(stat_arr_raw(k,1,j,i)-nanmean(stat_arr_raw(:,1,j,i)))/nanmean(stat_arr_raw(:,1,j,i)); % baseline
            end
        end
        for k=1:Tf_st % for each frame
            if nanmean(stat_arr_raw(:,1,j,i))<0
                stat_arr_dF(k,2,j,i)=(stat_arr_raw(k,2,j,i)-nanmean(stat_arr_raw(:,1,j,i)))/((-1)*nanmean(stat_arr_raw(:,1,j,i)));
            else
                stat_arr_dF(k,2,j,i)=(stat_arr_raw(k,2,j,i)-nanmean(stat_arr_raw(:,1,j,i)))/nanmean(stat_arr_raw(:,1,j,i)); % post stimon
            end
        end
    end
end

%%
stat_arr=nan(Ts,2,Tr);
% enter correct values (:,1,:) for baseline, (:,2,:) for post stim onset
for i=1:Tr % for each ROI
    for j=1:Ts % for each stim
%         if (StimOn(j)-Tf_pre)>=1
            stat_arr(j,1,i)=nanmean(stat_arr_dF(:,1,j,i));
%         end
        stat_arr(j,2,i)=nanmean(stat_arr_dF(:,2,j,i));
        % 2 added for ~ half a second shift into the stim
    end
end



%%
% Wilcoxon sign rank test (for paired non uniform samples) since all roi
% signal distribution of either baseline or post stim was not normal
% make blank signrank stat cell variable
sr.val={};
% binary inputs for all rois, 1= responsive, else 0
resp=zeros(1,Tr);
% index of responsive rois
resp_rois=[];
% incrementing count for responsive rois
c=1;
for i=1:Tr % for each roi
    % tail is left since signrank(x,y) y>x, we expect the signal to
    % increase after the stim onset. p = p-value, h = hypothesis (1
    % indicates null hypothesis is rejected else 0).
    [p,h,stats] = signrank(stat_arr(:,1,i),stat_arr(:,2,i),'alpha',0.05,...
        'tail','left');
    sr.val.p(i)=p;
    if sr.val.p(i)<0.05
        resp(1,i)=1;
        resp_rois(c)=i;
        c=c+1;
    end
    sr.val.h(i)=h;
    sr.val.stat(i)=stats;
end

%
disp(['No. of responsive rois: ', num2str(size(resp_rois,2))]);
disp(['Percent of responsive rois: ',num2str(round(100*size(resp_rois,2)...
    /Tr),2)]);
%%
save(fn, 'resp','resp_rois','sr','StimOn','stat_arr','speed', 'Tf_pre',...
    'Tf_st','stat_arr_dF','-append')
clear


