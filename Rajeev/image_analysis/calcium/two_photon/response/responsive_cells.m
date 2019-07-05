%% created by Rajeev Rajendran
% 2019_06_27 
% Open directory
% split folders
[fn , pn] = uigetfile('*_SPSIG.mat');
cd(pn)

%% Load data files
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

%% manually check by opening stim
prompt='Which is the first stim onset index?';
stim1=input(prompt); % first stim onset
prompt='Which is the last stim onset index?';
stim2=input(prompt); % last stim onset
n=size(SPSIG_file.sig,1); % total frames available

% Selecting the stim on and off frames
Ts=size((stim1:2:stim2),2);
for i=1:Ts
    StimOn(i,1)=stim((stim1+(2*i-2)),1);
    StimOff(i,1)=stim((stim1+(2*i-2))+1,1);
end

% label frames when stim is on
m=1;
for i=1:n
    if i >= StimOn(m,1) && i < StimOff(m,1)
        st_fr(i,1)=1;
    elseif i == StimOff(m,1)
        st_fr(i,1)=1;
        if m < Ts
            m=m+1;
        end
    else
        st_fr(i,1)=0;
    end
end
%% Extract stim frames
% Prestim = 4s, stim 2s frames, post stim 8s
prompt='how many seconds before stim onset?';
pre_stim_dur=input(prompt); % in seconds
prompt='how many seconds of stim?';
stim_dur=input(prompt); % in seconds
prompt='how many seconds after stim offset?';
post_stim_dur=input(prompt); % in seconds
%% distribute the signals into pre and stim groups for all stims
Tf_pre=ceil((pre_stim_dur)*SPSIG_file.freq); % no. of frames before stim onset
Tf_st=ceil((stim_dur+post_stim_dur)*SPSIG_file.freq); % no. of frames from stim onset
Tr=size(SPSIG_file.sig,2);
stat_arr=nan(Ts,2,Tr);
for i=1:Tr % for each ROI
    for j=1:Ts % for each stim
        if (StimOn(j)-Tf_pre)>=1
            stat_arr(j,1,i)=mean(SPSIG_file.sig((StimOn(j)-Tf_pre):(StimOn(j)-1),i));
        end
            stat_arr(j,2,i)=mean(SPSIG_file.sig(StimOn(j):(StimOn(j)+Tf_st),i));  
    end
end

%% shapiro wilk
% sw.stat={};
% for i=1:Tr % for each ROI
%     [H, pValue, SWstatistic]=swtest(stat_arr(:,1,i),0.05);
%     sw.stat.H(i,1)=H;
%     sw.stat.p(i,1)=pValue;
%     sw.stat.SW(i,1)=SWstatistic;
%     [H, pValue, SWstatistic]=swtest(stat_arr(:,2,i),0.05);
%     sw.stat.H(i,2)=H;
%     sw.stat.p(i,2)=pValue;
%     sw.stat.SW(i,2)=SWstatistic;
% end
%% Wilcoxon rank sum test (for non-paired non uniform samples) 
% % equivalent to the Mann-Whitney U test
% 
% rs.val={};
% resp=zeros(1,Tr);
% resp_rois=[];
% c=1;%count
% for i=1:Tr % for each ROI
%     [p,h,stats]=ranksum(stat_arr(:,1,i),stat_arr(:,2,i),'alpha',0.05,'tail','left');
%     rs.val.p(i)=p;
%     if rs.val.p(i)<0.05
%         resp(1,i)=1;
%         resp_rois(c)=i;
%         c=c+1;
%     end
%     rs.val.h(i)=h;
%     rs.val.stat(i)=stats;
% end
% tot_roi=sum(rs.val.h);

%% Wilcoxon sign rank test (for paired non uniform samples) 
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

%%
save(fn, 'resp','resp_rois','sr','StimOn','stat_arr','speed', 'Tf_pre', 'Tf_st', '-append')
clear
clc
