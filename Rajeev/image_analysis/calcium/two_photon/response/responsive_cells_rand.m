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
% Prestim = 4s, stim 2s frames, post stim 4s
prompt='how many seconds before stim onset?';
pre_stim_dur=input(prompt); % in seconds
prompt='how many seconds of stim?';
stim_dur=input(prompt); % in seconds
prompt='how many seconds after stim offset?';
post_stim_dur=input(prompt); % in seconds

%%
random_roi={};
tic
for rpt=1:1000
    % define variables
    Tf_pre=ceil((pre_stim_dur)*SPSIG_file.freq); % no. of frames before stim onset
    Tf_st=ceil((stim_dur+post_stim_dur)*SPSIG_file.freq); % no. of frames from stim onset
    Tr=size(SPSIG_file.sig,2);
    % randomize stims
    a=(floor(stim(stim1,1)/100))*100;
    if a < 1
       a=1; 
    end
    b=size(speed,2)-(1+Tf_st);
    r = sort(round(a + (b-a).*rand(size((stim1:2:stim2),2),1)));
    
    % Selecting the stim on and off frames
    Ts=size((stim1:2:stim2),2);
    for i=1:Ts
        StimOn(i,1)=r(i,1);
        StimOff(i,1)=StimOn(i,1)-10;
    end
    
    pre_stim_dur=4;
    stim_dur=2;
    post_stim_dur=4;
    
    % distribute the signals into pre and stim groups for all stims
    
    stat_arr=nan(Ts,2,Tr);
    for i=1:Tr % for each ROI
        for j=1:Ts % for each stim
            if (StimOn(j)-Tf_pre)>=1
                stat_arr(j,1,i)=mean(SPSIG_file.sig((StimOn(j)-Tf_pre):(StimOn(j)-1),i));
            end
            stat_arr(j,2,i)=mean(SPSIG_file.sig(StimOn(j):(StimOn(j)+Tf_st),i));
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
    %     random_roi.resp.rpt=resp_rois;
    random_roi.n(rpt)=size(resp_rois,2);
end

%
prompt='how many rois are responsive to actual stim?';
true_resp=input(prompt);
c=0;
for i=1:1000
    if random_roi.n(1,i) < true_resp
        c=c+1;
    end
end
random_roi.percentile=c/10;
random_roi.mean=mean(random_roi.n(1,:));
random_roi.std=std(random_roi.n(1,:));
random_roi.med=median(random_roi.n(1,:));
histogram(random_roi.n(1,:))
toc
save('rand_test','random_roi','speed','run','stim')

