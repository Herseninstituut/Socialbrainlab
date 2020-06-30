%% created by Steven Voges based on script by Rajeev Rajendran
% 2020_06_29
% First finds selection of responsive rois to any stimulus using Wilcoxon
% Signed Rank one Tailed 
% Then using this selection compares these stimuli to each other to see
% which are significant over their controls
% Stim onset times should be stored in stims_split.stimname using extractstims.mat  (Full list: shock, shockCTRL,
% laserHI, laserLO, laserCTRL, squeak, squeakCTRL)
% Data should be from caiSig file extracted using extractcaiman.mat

%% load caiman data file and stim time file


%%


n=size(caiSig,1); % total frames available

% select stim on and off frames
Ts=10

pre_stim_dur=2 %Number of seconds before stim
stim_dur = 2 %stim dueration
post_stim_dur = 2 %number of seconds after stim
freq = 5.1666
% distribute the signals into pre and stim groups for all stims
% no. of frames before stim onset
Tf_pre=ceil((pre_stim_dur)*freq);
% no. of frames from stim onset
Tf_st=ceil((stim_dur+post_stim_dur)*freq);
% no. of ROIs
Tr=size(caiSig,2)
% make blank array for stats
stat_arr.shock=nan(Ts,2,Tr);
% enter correct values (:,1,:) for baseline, (:,2,:) for post stim onset
for i=1:Tr % for each ROI
    for j=1:Ts % for each stim
        if (stims_split.shock(j)-Tf_pre)>=1
            stat_arr.shock(j,1,i)=mean(caiSig((stims_split.shock(j)-Tf_pre):...
                (stims_split.shock(j)-1),i));
        end
        stat_arr.shock(j,2,i)=mean(caiSig(stims_split.shock(j)+2:(stims_split.shock(j)...
            +Tf_st),i));
        % 2 added for ~ half a second shift into the stim
    end
end


% Wilcoxon sign rank test (for paired non uniform samples) since all roi
% signal distribution of either baseline or post stim was not normal
% make blank signrank stat cell variable
sr.val={};
% binary inputs for all rois, 1= responsive, else 0
resp.shock=zeros(1,Tr);
% index of responsive rois
resp_rois.shock=[];
% incrementing count for responsive rois
c=1;
for i=1:Tr % for each roi
    % tail is left since signrank(x,y) y>x, we expect the signal to
    % increase after the stim onset. p = p-value, h = hypothesis (1
    % indicates null hypothesis is rejected else 0).
    [p,h,stats] = signrank(stat_arr.shock(:,1,i),stat_arr.shock(:,2,i),'alpha',0.05,...
        'tail','left');
    sr.val.p(i)=p;
    if sr.val.p(i)<0.05
        resp.shock(1,i)=1;
        resp_rois.shock(c)=i;
        c=c+1;
    end
    sr.val.h(i)=h;
    sr.val.stat(i)=stats;
end

disp('Shock ')
disp(['No. of responsive rois: ', num2str(size(resp_rois.shock,2))]);
disp(['Percent of responsive rois: ',num2str(round(100*size(resp_rois.shock,2)...
    /Tr),2)]);
%%

stat_arr.shockCTRL=nan(Ts,2,Tr);
% enter correct values (:,1,:) for baseline, (:,2,:) for post stim onset
for i=1:Tr % for each ROI
    for j=1:Ts % for each stim
        if (stims_split.shockCTRL(j)-Tf_pre)>=1
            stat_arr.shockCTRL(j,1,i)=mean(caiSig((stims_split.shockCTRL(j)-Tf_pre):...
                (stims_split.shockCTRL(j)-1),i));
        end
        stat_arr.shockCTRL(j,2,i)=mean(caiSig(stims_split.shockCTRL(j)+2:(stims_split.shockCTRL(j)...
            +Tf_st),i));
        % 2 added for ~ half a second shift into the stim
    end
end


% Wilcoxon sign rank test (for paired non uniform samples) since all roi
% signal distribution of either baseline or post stim was not normal
% make blank signrank stat cell variable
sr.val={};
% binary inputs for all rois, 1= responsive, else 0
resp.shockCTRL=zeros(1,Tr);
% index of responsive rois
resp_rois.shockCTRL=[];
% incrementing count for responsive rois
c=1;
for i=1:Tr % for each roi
    % tail is left since signrank(x,y) y>x, we expect the signal to
    % increase after the stim onset. p = p-value, h = hypothesis (1
    % indicates null hypothesis is rejected else 0).
    [p,h,stats] = signrank(stat_arr.shockCTRL(:,1,i),stat_arr.shockCTRL(:,2,i),'alpha',0.05,...
        'tail','left');
    sr.val.p(i)=p;
    if sr.val.p(i)<0.05
        resp.shockCTRL(1,i)=1;
        resp_rois.shockCTRL(c)=i;
        c=c+1;
    end
    sr.val.h(i)=h;
    sr.val.stat(i)=stats;
end


%
disp('Shock Ctrl')
disp(['No. of responsive rois: ', num2str(size(resp_rois.shockCTRL,2))]);
disp(['Percent of responsive rois: ',num2str(round(100*size(resp_rois.shockCTRL,2)...
    /Tr),2)]);


%%
pre_stim_dur=2
stim_dur = 2
post_stim_dur = 2
freq = 5.1666
% distribute the signals into pre and stim groups for all stims
% no. of frames before stim onset
Tf_pre=ceil((pre_stim_dur)*freq);
% no. of frames from stim onset
Tf_st=ceil((stim_dur+post_stim_dur)*freq);

Ts=9
stat_arr.laserHI=nan(Ts,2,Tr);
% enter correct values (:,1,:) for baseline, (:,2,:) for post stim onset
for i=1:Tr % for each ROI
    for j=1:Ts % for each stim
        if (stims_split.laserHI(j)-Tf_pre)>=1
            stat_arr.laserHI(j,1,i)=mean(caiSig((stims_split.laserHI(j)-Tf_pre):...
                (stims_split.laserHI(j)-1),i));
        end
        stat_arr.laserHI(j,2,i)=mean(caiSig(stims_split.laserHI(j)+2:(stims_split.laserHI(j)...
            +Tf_st),i));
        % 2 added for ~ half a second shift into the stim
    end
end


% Wilcoxon sign rank test (for paired non uniform samples) since all roi
% signal distribution of either baseline or post stim was not normal
% make blank signrank stat cell variable
sr.val={};
% binary inputs for all rois, 1= responsive, else 0
resp.laserHI=zeros(1,Tr);
% index of responsive rois
resp_rois.laserHI=[];
% incrementing count for responsive rois
c=1;
for i=1:Tr % for each roi
    % tail is left since signrank(x,y) y>x, we expect the signal to
    % increase after the stim onset. p = p-value, h = hypothesis (1
    % indicates null hypothesis is rejected else 0).
    [p,h,stats] = signrank(stat_arr.laserHI(:,1,i),stat_arr.laserHI(:,2,i),'alpha',0.05,...
        'tail','left');
    sr.val.p(i)=p;
    if sr.val.p(i)<0.05
        resp.laserHI(1,i)=1;
        resp_rois.laserHI(c)=i;
        c=c+1;
    end
    sr.val.h(i)=h;
    sr.val.stat(i)=stats;
end


%
disp('LaserHI')
disp(['No. of responsive rois: ', num2str(size(resp_rois.laserHI,2))]);
disp(['Percent of responsive rois: ',num2str(round(100*size(resp_rois.laserHI,2)...
    /Tr),2)]);
%save(fn, 'resp','resp_rois','sr','StimOn','stat_arr','speed', 'Tf_pre',...
   % 'Tf_st', '-append')
%clear


%%


Ts=10
stat_arr.laserLO=nan(Ts,2,Tr);
% enter correct values (:,1,:) for baseline, (:,2,:) for post stim onset
for i=1:Tr % for each ROI
    for j=1:Ts % for each stim
        if (stims_split.laserLO(j)-Tf_pre)>=1
            stat_arr.laserLO(j,1,i)=mean(caiSig((stims_split.laserLO(j)-Tf_pre):...
                (stims_split.laserLO(j)-1),i));
        end
        stat_arr.laserLO(j,2,i)=mean(caiSig(stims_split.laserLO(j)+2:(stims_split.laserLO(j)...
            +Tf_st),i));
        % 2 added for ~ half a second shift into the stim
    end
end


% Wilcoxon sign rank test (for paired non uniform samples) since all roi
% signal distribution of either baseline or post stim was not normal
% make blank signrank stat cell variable
sr.val={};
% binary inputs for all rois, 1= responsive, else 0
resp.laserLO=zeros(1,Tr);
% index of responsive rois
resp_rois.laserLO=[];
% incrementing count for responsive rois
c=1;
for i=1:Tr % for each roi
    % tail is left since signrank(x,y) y>x, we expect the signal to
    % increase after the stim onset. p = p-value, h = hypothesis (1
    % indicates null hypothesis is rejected else 0).
    [p,h,stats] = signrank(stat_arr.laserLO(:,1,i),stat_arr.laserLO(:,2,i),'alpha',0.05,...
        'tail','left');
    sr.val.p(i)=p;
    if sr.val.p(i)<0.05
        resp.laserLO(1,i)=1;
        resp_rois.laserLO(c)=i;
        c=c+1;
    end
    sr.val.h(i)=h;
    sr.val.stat(i)=stats;
end


%
disp('LaserLO')

disp(['No. of responsive rois: ', num2str(size(resp_rois.laserLO,2))]);
disp(['Percent of responsive rois: ',num2str(round(100*size(resp_rois.laserLO,2)...
    /Tr),2)]);



%%
Ts=10
stat_arr.laserCTRL=nan(Ts,2,Tr);
% enter correct values (:,1,:) for baseline, (:,2,:) for post stim onset
for i=1:Tr % for each ROI
    for j=1:Ts % for each stim
        if (stims_split.laserCTRL(j)-Tf_pre)>=1
            stat_arr.laserCTRL(j,1,i)=mean(caiSig((stims_split.laserCTRL(j)-Tf_pre):...
                (stims_split.laserCTRL(j)-1),i));
        end
        stat_arr.laserCTRL(j,2,i)=mean(caiSig(stims_split.laserCTRL(j)+2:(stims_split.laserCTRL(j)...
            +Tf_st),i));
        % 2 added for ~ half a second shift into the stim
    end
end


% Wilcoxon sign rank test (for paired non uniform samples) since all roi
% signal distribution of either baseline or post stim was not normal
% make blank signrank stat cell variable
sr.val={};
% binary inputs for all rois, 1= responsive, else 0
resp.laserCTRL=zeros(1,Tr);
% index of responsive rois
resp_rois.laserCTRL=[];
% incrementing count for responsive rois
c=1;
for i=1:Tr % for each roi
    % tail is left since signrank(x,y) y>x, we expect the signal to
    % increase after the stim onset. p = p-value, h = hypothesis (1
    % indicates null hypothesis is rejected else 0).
    [p,h,stats] = signrank(stat_arr.laserCTRL(:,1,i),stat_arr.laserCTRL(:,2,i),'alpha',0.05,...
        'tail','left');
    sr.val.p(i)=p;
    if sr.val.p(i)<0.05
        resp.laserCTRL(1,i)=1;
        resp_rois.laserCTRL(c)=i;
        c=c+1;
    end
    sr.val.h(i)=h;
    sr.val.stat(i)=stats;
end


%
disp('LaserCTRL')
disp(['No. of responsive rois: ', num2str(size(resp_rois.laserCTRL,2))]);
disp(['Percent of responsive rois: ',num2str(round(100*size(resp_rois.laserCTRL,2)...
    /Tr),2)]);


%%
n=size(caiSig,1); % total frames available

% select stim on and off frames
Ts=10

pre_stim_dur=2 %Number of seconds before stim
stim_dur = 2 %stim dueration
post_stim_dur = 2 %number of seconds after stim
freq = 5.1666
% distribute the signals into pre and stim groups for all stims
% no. of frames before stim onset
Tf_pre=ceil((pre_stim_dur)*freq);
% no. of frames from stim onset
Tf_st=ceil((stim_dur+post_stim_dur)*freq);
% no. of ROIs
Tr=size(caiSig,2)
% make blank array for stats
stat_arr.squeak=nan(Ts,2,Tr);
% enter correct values (:,1,:) for baseline, (:,2,:) for post stim onset
for i=1:Tr % for each ROI
    for j=1:Ts % for each stim
        if (stims_split.squeak(j)-Tf_pre)>=1
            stat_arr.squeak(j,1,i)=mean(caiSig((stims_split.squeak(j)-Tf_pre):...
                (stims_split.squeak(j)-1),i));
        end
        stat_arr.squeak(j,2,i)=mean(caiSig(stims_split.squeak(j)+2:(stims_split.squeak(j)...
            +Tf_st),i));
        % 2 added for ~ half a second shift into the stim
    end
end


% Wilcoxon sign rank test (for paired non uniform samples) since all roi
% signal distribution of either baseline or post stim was not normal
% make blank signrank stat cell variable
sr.val={};
% binary inputs for all rois, 1= responsive, else 0
resp.squeak=zeros(1,Tr);
% index of responsive rois
resp_rois.squeak=[];
% incrementing count for responsive rois
c=1;
for i=1:Tr % for each roi
    % tail is left since signrank(x,y) y>x, we expect the signal to
    % increase after the stim onset. p = p-value, h = hypothesis (1
    % indicates null hypothesis is rejected else 0).
    [p,h,stats] = signrank(stat_arr.squeak(:,1,i),stat_arr.squeak(:,2,i),'alpha',0.05,...
        'tail','left');
    sr.val.p(i)=p;
    if sr.val.p(i)<0.05
        resp.squeak(1,i)=1;
        resp_rois.squeak(c)=i;
        c=c+1;
    end
    sr.val.h(i)=h;
    sr.val.stat(i)=stats;
end

disp('squeak ')
disp(['No. of responsive rois: ', num2str(size(resp_rois.squeak,2))]);
disp(['Percent of responsive rois: ',num2str(round(100*size(resp_rois.squeak,2)...
    /Tr),2)]);

%%
n=size(caiSig,1); % total frames available

% select stim on and off frames
Ts=10

pre_stim_dur=2 %Number of seconds before stim
stim_dur = 2 %stim dueration
post_stim_dur = 2 %number of seconds after stim
freq = 5.1666
% distribute the signals into pre and stim groups for all stims
% no. of frames before stim onset
Tf_pre=ceil((pre_stim_dur)*freq);
% no. of frames from stim onset
Tf_st=ceil((stim_dur+post_stim_dur)*freq);
% no. of ROIs
Tr=size(caiSig,2)
% make blank array for stats
stat_arr.squeakCTRL=nan(Ts,2,Tr);
% enter correct values (:,1,:) for baseline, (:,2,:) for post stim onset
for i=1:Tr % for each ROI
    for j=1:Ts % for each stim
        if (stims_split.squeakCTRL(j)-Tf_pre)>=1
            stat_arr.squeakCTRL(j,1,i)=mean(caiSig((stims_split.squeakCTRL(j)-Tf_pre):...
                (stims_split.squeakCTRL(j)-1),i));
        end
        stat_arr.squeakCTRL(j,2,i)=mean(caiSig(stims_split.squeakCTRL(j)+2:(stims_split.squeakCTRL(j)...
            +Tf_st),i));
        % 2 added for ~ half a second shift into the stim
    end
end


% Wilcoxon sign rank test (for paired non uniform samples) since all roi
% signal distribution of either baseline or post stim was not normal
% make blank signrank stat cell variable
sr.val={};
% binary inputs for all rois, 1= responsive, else 0
resp.squeakCTRL=zeros(1,Tr);
% index of responsive rois
resp_rois.squeakCTRL=[];
% incrementing count for responsive rois
c=1;
for i=1:Tr % for each roi
    % tail is left since signrank(x,y) y>x, we expect the signal to
    % increase after the stim onset. p = p-value, h = hypothesis (1
    % indicates null hypothesis is rejected else 0).
    [p,h,stats] = signrank(stat_arr.squeakCTRL(:,1,i),stat_arr.squeakCTRL(:,2,i),'alpha',0.05,...
        'tail','left');
    sr.val.p(i)=p;
    if sr.val.p(i)<0.05
        resp.squeakCTRL(1,i)=1;
        resp_rois.squeakCTRL(c)=i;
        c=c+1;
    end
    sr.val.h(i)=h;
    sr.val.stat(i)=stats;
end

disp('squeakCTRL ')
disp(['No. of responsive rois: ', num2str(size(resp_rois.squeakCTRL,2))]);
disp(['Percent of responsive rois: ',num2str(round(100*size(resp_rois.squeakCTRL,2)...
    /Tr),2)]);


%% all responsives

resp_rois.all =  unique([resp_rois.shock resp_rois.shockCTRL resp_rois.laserHI resp_rois.laserLO resp_rois.laserCTRL resp_rois.squeak resp_rois.squeakCTRL])

%%Test if shock is higher than shock ctrl for all resp rois


sr.val={};
sig_rois=[]
sig_rois_ind=[];
c=1;

st_pre.shock = stat_arr.shock(:,2,:)-stat_arr.shock(:,1,:);

st_pre.shockCTRL = stat_arr.shockCTRL(:,2,:)-stat_arr.shock(:,1,:);

for i=1:length(resp_rois.all)% for each ROI
    %     if isfinite(st_pre_nan(:,2,i)) && isfinite(st_pre_nan(:,1,i))
    [p,h,stats] = signrank(st_pre.shock(:,resp_rois.all(i)),st_pre.shockCTRL(:,resp_rois.all(i)),'alpha',0.05,'tail','right');
    sr.val.p(i)=p;
    if sr.val.p(i)<0.05
        sig_rois(1,i)=1;
        sig_rois_ind(c,1)=i;
        c=c+1;
    end
    sr.val.h(i)=h;
    sr.val.stat(i)=stats;
    %     end
end

resp_rois.shockOVERctrl = resp_rois.all(nonzeros(sig_rois_ind))

%% Test if laserHI is higher than laserCTRL for all resp rois
sr.val={};
sig_rois=[]
sig_rois_ind=[];
c=1;
st_pre.laserHI = stat_arr.laserHI(:,2,:)-stat_arr.laserHI(:,1,:);

st_pre.laserCTRL = stat_arr.laserCTRL(:,2,:)-stat_arr.laserCTRL(:,1,:);
c=1 
for i=1:length(resp_rois.all)% for each ROI
    %     if isfinite(st_pre_nan(:,2,i)) && isfinite(st_pre_nan(:,1,i))
    [p,h,stats] = signrank(st_pre.laserHI(:,resp_rois.all(i)),st_pre.laserCTRL(1:9,resp_rois.all(i)),'alpha',0.05,'tail','right');
    sr.val.p(i)=p;
    if sr.val.p(i)<0.05
        sig_rois(1,i)=1;
        sig_rois_ind(c,1)=i;
        c=c+1;
    end
    sr.val.h(i)=h;
    sr.val.stat(i)=stats;
    %     end
end

resp_rois.laserHIOVERctrl = resp_rois.all(nonzeros(sig_rois_ind))

%% Test if LaserLO is higher than laser control for all resp_rois
sr.val={};
sig_rois=[]
sig_rois_ind=[];
c=1;
st_pre.laserLO = stat_arr.laserLO(:,2,:)-stat_arr.laserLO(:,1,:);

st_pre.laserCTRL = stat_arr.laserCTRL(:,2,:)-stat_arr.laserCTRL(:,1,:);
c=1; 
for i=1:length(resp_rois.all)% for each ROI
    %     if isfinite(st_pre_nan(:,2,i)) && isfinite(st_pre_nan(:,1,i))
    [p,h,stats] = signrank(st_pre.laserLO(:,resp_rois.all(i)),st_pre.laserCTRL(:,resp_rois.all(i)),'alpha',0.05,'tail','right');
    sr.val.p(i)=p;
    if sr.val.p(i)<0.05
        sig_rois(1,i)=1;
        sig_rois_ind(c,1)=i;
        c=c+1;
    end
    sr.val.h(i)=h;
    sr.val.stat(i)=stats;
    %     end
end

resp_rois.laserLOOVERctrl = resp_rois.all(nonzeros(sig_rois_ind))



%% Test if laserHI is higher than laserLO for all resp_rois
sr.val={};
sig_rois=[]
sig_rois_ind=[];

c=1; 
for i=1:length(resp_rois.all)% for each ROI
    %     if isfinite(st_pre_nan(:,2,i)) && isfinite(st_pre_nan(:,1,i))
    [p,h,stats] = signrank(st_pre.laserHI(:,resp_rois.all(i)),st_pre.laserLO(1:9,resp_rois.all(i)),'alpha',0.05,'tail','right');
    sr.val.p(i)=p;
    if sr.val.p(i)<0.05
        sig_rois(1,i)=1;
        sig_rois_ind(c,1)=i;
        c=c+1;
    end
    sr.val.h(i)=h;
    sr.val.stat(i)=stats;
    %     end
end

resp_rois.laserHIOVERLO = resp_rois.all(sig_rois_ind)

%% test if laserLO is higher than laser hi for all resp rois
sr.val={};
sig_rois=[]
sig_rois_ind=[];
c=1; 

for i=1:length(resp_rois.all)% for each ROI
    %     if isfinite(st_pre_nan(:,2,i)) && isfinite(st_pre_nan(:,1,i))
    [p,h,stats] = signrank(st_pre.laserLO(1:9,resp_rois.all(i)),st_pre.laserHI(:,resp_rois.all(i)),'alpha',0.05,'tail','right');
    sr.val.p(i)=p;
    if sr.val.p(i)<0.05
        sig_rois(1,i)=1;
        sig_rois_ind(c,1)=i;
        c=c+1;
    end
    sr.val.h(i)=h;
    sr.val.stat(i)=stats;
    %     end
end

resp_rois.laserLOOVERHI = resp_rois.all(nonzeros(sig_rois_ind))

%% Test if laser ctrl is higher than laserHI for all resp rois
sr.val={};
sig_rois=[]
sig_rois_ind=[];
c=1; 

for i=1:length(resp_rois.all)% for each ROI
    %     if isfinite(st_pre_nan(:,2,i)) && isfinite(st_pre_nan(:,1,i))
    [p,h,stats] = signrank(st_pre.laserCTRL(1:9,resp_rois.all(i)),st_pre.laserHI(:,resp_rois.all(i)),'alpha',0.05,'tail','right');
    sr.val.p(i)=p;
    if sr.val.p(i)<0.05
        sig_rois(1,i)=1;
        sig_rois_ind(c,1)=i;
        c=c+1;
    end
    sr.val.h(i)=h;
    sr.val.stat(i)=stats;
    %     end
end

resp_rois.laserCTRLOVERHI = resp_rois.all(nonzeros(sig_rois_ind))
%% Test if squeak is higher than squeak controll for all resp rois
sr.val={};
sig_rois=[]
sig_rois_ind=[];
c=1;

st_pre.squeak = stat_arr.squeak(:,2,:)-stat_arr.squeak(:,1,:);

st_pre.squeakCTRL = stat_arr.squeakCTRL(:,2,:)-stat_arr.squeakCTRL(:,1,:);

for i=1:length(resp_rois.all)% for each ROI
    %     if isfinite(st_pre_nan(:,2,i)) && isfinite(st_pre_nan(:,1,i))
    [p,h,stats] = signrank(st_pre.squeak(:,resp_rois.all(i)),st_pre.squeakCTRL(:,resp_rois.all(i)),'alpha',0.05,'tail','right');
    sr.val.p(i)=p;
    if sr.val.p(i)<0.05
        sig_rois(1,i)=1;
        sig_rois_ind(c,1)=i;
        c=c+1;
    end
    sr.val.h(i)=h;
    sr.val.stat(i)=stats;
    %     end
end

resp_rois.squeakOVERctrl = resp_rois.all(nonzeros(sig_rois_ind))





