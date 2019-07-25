%% created by Rajeev Rajendran
% 2019_07_02 
% Open directory
pn = uigetdir();
cd(pn)

%%
cd ([pn '\mult_sessions\'])
load('chr_matched_rois_split1','linkMat2');
cd ([pn '\shock_observation\analysis\split1'])
sig{1}=load('glm_speed_regr');
cd ([pn '\shock_control\analysis\split1'])
sig{2}=load('glm_speed_regr');
cd ([pn '\laser_self_select\analysis\split1'])
sig{3}=load('glm_speed_regr');
cd ([pn '\squeak_playback\analysis\split1'])
sig{4}=load('glm_speed_regr');
% decon{1}=decon1;
% decon{2}=decon2;
% decon{3}=decon3;
% decon{4}=decon4;

%% find number of cells common to all sessions
for i = 1:size(linkMat2,1)
    for j=1:size(linkMat2,2)
        if linkMat2(i,j)~=0
            continue
        else
            num=i-1; % number of cells common to all sessions
            break
        end
    end
    if linkMat2(i,j)==0
        break
    end
end

%% select rois responsive in atleast one session
sess_rois=zeros(num,1); % logical
sess_rois_num=[]; % actual numbers only
c=1; % count
for i=1:num % no. of rois
    for j=1:size(linkMat2,2) % number of sessions
        if sig{1,j}.resp(1,linkMat2(i,j))==1
            sess_rois(i,1)=1;
            sess_rois_num(c,1)=i;
            c=c+1;
            break
        end
    end
end

%% collate stim-pre values of resp rois
x=0;
for j=1:size(linkMat2,2) % index of sessions
    if size(sig{1, j}.stat_arr,1)>x
        x=j;
    end
end
st_pre=nan(x,size(linkMat2,2),size(sess_rois_num,1));
for i=1:size(sess_rois_num,1) % no. of responsive rois
    for j=1:size(linkMat2,2) % index of sessions
        for k=1:size(sig{1, j}.stat_arr,1) % number of trials
            st_pre(k,j,i)=sig{1, j}.stat_arr(k,2,linkMat2(sess_rois_num(i,1),j))-sig{1, j}.stat_arr(k,1,linkMat2(sess_rois_num(i,1),j));
        end
    end
end

% replace zeros (not actual measurements) with nans
st_pre_nan=nan(size(st_pre,1),size(st_pre,2),size(st_pre,3));
for i=1:size(st_pre,3) % no. of responsive rois
    for j=1:size(st_pre,2) % index of sessions
        for k=1:size(st_pre,1) % number of trials
            if st_pre(k,j,i)~=0 || isnan(st_pre(k,j,i))
                st_pre_nan(k,j,i)=st_pre(k,j,i);
            end
        end
    end
end

%% compare rois shock observation > shock control
sr.val={};
sig_rois=zeros(1,size(st_pre_nan,3));
sig_rois_ind=[];
c=1;
for i=1:size(st_pre_nan,3) % for each ROI
%     if isfinite(st_pre_nan(:,2,i)) && isfinite(st_pre_nan(:,1,i))
        [p,h,stats] = signrank(st_pre_nan(:,2,i),st_pre_nan(:,1,i),'alpha',0.05,'tail','left');
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

% check whether the shock observation > shock control rois are also 
% responsive for laser or squeak

sig_rois_ind_ori=[];
for i=1:size(sig_rois_ind,1)
    sig_rois_ind_ori(i,1)=linkMat2(sess_rois_num(sig_rois_ind(i,1),1),1);
    if any(sig{1, 3}.resp_rois(1,:)==linkMat2(sess_rois_num(sig_rois_ind(i,1),1),3))
        sig_rois_ind_ori(i,2)=linkMat2(sess_rois_num(sig_rois_ind(i,1),1),3);
    end
    if any(sig{1, 4}.resp_rois(1,:)==linkMat2(sess_rois_num(sig_rois_ind(i,1),1),4))
        sig_rois_ind_ori(i,3)=linkMat2(sess_rois_num(sig_rois_ind(i,1),1),4);
    end
end

%% save values

save('comm_rois','sig_rois','sig_rois_ind','sig_rois_ind_ori','sess_rois','sess_rois_num')
clear
clc



