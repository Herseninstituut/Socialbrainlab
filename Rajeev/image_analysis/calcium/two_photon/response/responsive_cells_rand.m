%% created by Rajeev Rajendran
% 2019_06_27
% Shuffles the stim onset amonst the available frames and tests
% resposnsiveness of all the detected rois, and repeats it 1000 times to
% get an average idea of number of random responsive cells. Also see
% 'responsive_cells.m'. Correct stim indices and appropriate
% stim duration details are needed as input arguments.

function responsive_cells_rand(stim_first,stim_last,pre_stim_dur,...
    stim_dur,post_stim_dur)
format compact
% select the correct directory and make it current
[~ , pn] = uigetfile('*_SPSIG.mat');
cd(pn)

%% load data files
SPSIG_file = load(uigetfile('*_SPSIG.mat')); % signal file
info = load(uigetfile('*_normcorr.mat', 'info')); % normcorr mat file
quad = load(uigetfile('*quad_split*.mat')); % run file
stim=floor((info.info.frame)/3); % get stim onset and offset frames

% calculate speed
run=['quad.quad_',pn(end-6:end-1)];
run=eval(run);
ra = 10;
dist = (2*pi*ra)*run/1000; %in cms
speed = (dist)*SPSIG_file.freq; % in cms/s

n=size(SPSIG_file.sig,1); % total frames available

% empty cell to store the iteration data
random_roi={};
tic
for rpt=1:1000
    % define variables
    % no. of frames before stim onset
    Tf_pre=ceil((pre_stim_dur)*SPSIG_file.freq);
    % no. of frames from stim onset
    Tf_st=ceil((stim_dur+post_stim_dur)*SPSIG_file.freq);
    % no. of rois detected in preprocessing
    Tr=size(SPSIG_file.sig,2);
    % randomize stims. a is the starting frame and b the last frame for
    % stim onset
    a=(floor(stim(stim_first,1)/100))*100;
    if a < 1
        a=1;
    end
    b=size(speed,2)-(1+Tf_st);
    r = sort(round(a + (b-a).*rand(size((stim_first:2:stim_last),2),1)));
    
    % Selecting the stim on and off frames
    Ts=size((stim_first:2:stim_last),2);
    for i=1:Ts
        StimOn(i,1)=r(i,1);
        StimOff(i,1)=StimOn(i,1)-10;
    end
    
    % distribute the signals into pre and post stim onset groups for all
    % stims. enter correct values (:,1,:) for baseline, (:,2,:) for post
    % stim onset
    stat_arr=nan(Ts,2,Tr);
    for i=1:Tr % for each ROI
        for j=1:Ts % for each stim
            if (StimOn(j)-Tf_pre)>=1
                stat_arr(j,1,i)=mean(SPSIG_file.sig((StimOn(j)-Tf_pre):...
                    (StimOn(j)-1),i));
            end
            % 2 added for ~ half a second shift into the stim
            % condition added in case not enough frames are left after the
            % random last stim onset.
            
            if StimOn(j)+Tf_st > size(SPSIG_file.sig,1)
                stat_arr(j,2,i)=mean(SPSIG_file.sig(StimOn(j)+2:end,i));
            else
                stat_arr(j,2,i)=mean(SPSIG_file.sig(StimOn(j)+2:...
                    (StimOn(j)+Tf_st),i));
            end
        end
    end
    
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
        [p,h,stats] = signrank(stat_arr(:,1,i),stat_arr(:,2,i),'alpha',...
            0.05,'tail','left');
        sr.val.p(i)=p;
        if sr.val.p(i)<0.05
            resp(1,i)=1;
            resp_rois(c)=i;
            c=c+1;
        end
        sr.val.h(i)=h;
        sr.val.stat(i)=stats;
    end
    random_roi.n(rpt)=size(resp_rois,2);
end
toc
% number of actual responsive rois, see responsive_cells.m
true_resp=size(SPSIG_file.resp_rois,2);
c=0;
for i=1:1000
    % incrementing count for calculating percentile.
    if random_roi.n(1,i) < true_resp
        c=c+1;
    end
end
%
random_roi.percentile=c/10;
random_roi.mean=mean(random_roi.n(1,:));
random_roi.std=std(random_roi.n(1,:));
random_roi.med=median(random_roi.n(1,:));
save('rand_test','random_roi','speed','run','stim')
figure;
histogram(random_roi.n(1,:))
hold on
% to draw vertical line at actual responsive roi number.
y=ylim;
y=(1:y(2)/10:y(2));
x=ones(1,size(y,2));
x=x + true_resp;
plot(x,y,'color','r')
hold on
title('Randomized shuffle of stim onset to check chance of random responsive cell')
xlabel('number of responsive cells, red line = true responsive cell')
ylabel('number of shuffles')
% save png file in the current folder
saveas(gcf,'Random responsive cells.png')
disp(['Random probability of number of responsive cells are ',num2str(random_roi.mean),' +/- ',num2str(random_roi.std)]);
disp(['Percentile of actual response is  ',num2str(random_roi.percentile)]);
clear

