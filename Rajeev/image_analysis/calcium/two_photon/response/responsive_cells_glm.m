%% created by Rajeev Rajendran
% 2019_06_27 
% Open directory
% split folders
pn = uigetdir();
cd(pn)

%% Extract stim frames
% Prestim = 4s, stim 2s frames, post stim 4s
prompt='how many seconds before stim onset?';
pre_stim_dur=input(prompt); % in seconds
prompt='how many seconds of stim?';
stim_dur=input(prompt); % in seconds
prompt='how many seconds after stim offset?';
post_stim_dur=input(prompt); % in seconds

%% Load data files

load('glm_speed_regr.mat')
n=size(sig,1); % total frames available

% subtract the glm predicted values from sig
sig_regr=[];
Tr=size(sig,2); % number of ROIs
Ts=size(stim_on,1);
for i=1:Tr % for each ROI
    for j=1:n % for each frame
        if j < (fr_shift*fr_sh_interval+1)
            sig_regr(j,i)=nan;
        elseif j > size(pred(1).beta,1)+(fr_shift*fr_sh_interval)
            sig_regr(j,i)=nan;
        else
            sig_regr(j,i)=sig(j,i)-pred(i).beta(j-(fr_shift*fr_sh_interval),1);
        end
    end
end

% distribute the signals into pre and stim groups for all stims
Tf_pre=ceil((pre_stim_dur)*freq); % no. of frames before stim onset
Tf_st=ceil((stim_dur+post_stim_dur)*freq); % no. of frames from stim onset
stat_arr=nan(Ts,2,Tr);
for i=1:Tr % for each ROI
    for j=1:Ts % for each stim
        if (stim_on(j)-Tf_pre)>=1
            stat_arr(j,1,i)=mean(sig_regr((stim_on(j)-Tf_pre):(stim_on(j)-1),i));
        end
            stat_arr(j,2,i)=mean(sig_regr(stim_on(j)+2:(stim_on(j)+Tf_st),i));
            % +2 for ~ half a second shift into the stim
    end
end

% Wilcoxon sign rank test (for paired non uniform samples) 
sr.val={};
resp=zeros(1,Tr);
resp_rois=[];
c=1;
for i=1:Tr % for each ROI
    [p,h,stats] = signrank(stat_arr(:,1,i),stat_arr(:,2,i),'alpha',0.05,'tail','left');
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
save('glm_speed_regr', 'resp','resp_rois','sr','stat_arr', '-append')
clear
clc

