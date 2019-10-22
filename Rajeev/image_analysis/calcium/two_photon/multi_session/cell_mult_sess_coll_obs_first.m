%% created by Rajeev Rajendran
% 2019_07_02
% Open directory
function cell_mult_sess_coll()
pn = uigetdir();
cd(pn);
if exist([pn '\mult_sessions\split1\'])
else
    mkdir([pn '\mult_sessions\split1\'])
end
if exist([pn '\mult_sessions\split2\'])
else
    mkdir([pn '\mult_sessions\split2\'])
end
if exist([pn '\mult_sessions\split3\'])
else
    mkdir([pn '\mult_sessions\split3\'])
end
% go to correct split folder in mult_sessions

%% select mult_session split folder
pn1 = uigetdir();
cd(pn1);
% run=['quad.quad_',pn(end-6:end-1)];
fn1=(['chr_matched_rois_split',pn1(end)]);
load(fn1,'linkMat2');
% shock observation folder and file
pn_obs=([pn,'\shock_observation\analysis\split',pn1(end)]);
cd (pn_obs);
sig{1}=load(uigetfile('*_SPSIG.mat'));
% shock control folder and file
pn_con=([pn,'\shock_control\analysis\split',pn1(end)]);
cd (pn_con);
sig{2}=load(uigetfile('*_SPSIG.mat'));
% laser folder and file
if exist([pn,'\laser_self_select\analysis\split',pn1(end)])
    pn_las=([pn,'\laser_self_select\analysis\split',pn1(end)]);
else
    pn_las=([pn,'\laser_self\analysis\split',pn1(end)]);
end
cd (pn_las);
sig{3}=load(uigetfile('*_SPSIG.mat'));
% laser folder and file
pn_sq=([pn,'\squeak_playback\analysis\split',pn1(end)]);
cd (pn_sq);
sig{4}=load(uigetfile('*_SPSIG.mat'));

% decon{1}=decon1;
% decon{2}=decon2;
% decon{3}=decon3;
% decon{4}=decon4;
%% obs responsive cells also present in control
% collate dF values of responsive obs rois and same control

A=(sig{1, 1}.resp_rois)';
B=linkMat2(:,1);
[LIA,LOCB]=ismember(A,B);

%%
x=0;
for j=1:size(linkMat2,2) % index of sessions
    if size(sig{1, j}.stat_arr,1)>x
        x=size(sig{1, j}.stat_arr,1);
    end
end
c=1;
df.obs_df=[];
for i=1:size(sig{1,1}.resp_rois,2) % responsive rois
    if LOCB(i,1)~=0
        if linkMat2(LOCB(i,1),2)~=0
            for j=1:x % trials
                if j<=size(sig{1, 1}.stat_arr,1)
                    df.obs_df(j,2,c)=sig{1, 1}.stat_arr(j,2,i);
                else
                    df.obs_df(j,2,c)=NaN;
                end
                if j<=size(sig{1, 2}.stat_arr,1)
                    df.obs_df(j,1,c)=sig{1, 2}.stat_arr(j,2,i);
                else
                    df.obs_df(j,1,c)=NaN;
                end
            end
            df.roi(c,1)=A(i,1);
            c=c+1;
        end
    end
end
%
%
%
%
%


%% select rois responsive for observation compared to control

sr.val={};
sig_obs_rois=zeros(1,size(df.obs_df,3));
sig_obs_rois_ind=[];
c=1;
for i=1:size(df.obs_df,3) % for each ROI
    [p,h,stats] = signrank(df.obs_df(:,1,i),df.obs_df(:,2,i),'alpha',0.05,'tail','left');
    sr.val.p(i)=p;
    if sr.val.p(i)<0.05
        sig_obs_rois(1,i)=1;
        sig_obs_rois_ind(c,1)=i;
        c=c+1;
    end
    sr.obs.val.h(i)=h;
    sr.obs.val.stat(i)=stats;
end

% %%
% obs_rois=zeros(num,1); % logical
% obs_rois_num=[]; % actual numbers only
% c=1; % count
% for i=1:num % no. of rois
%     if sig{1,1}.resp(1,linkMat2(i,1))==1
%         obs_rois(i,1)=1;
%         obs_rois_num(c,1)=i;
%         c=c+1;
%         break
%     end
% end
% 
% 
% 
% 
% %% find number of cells common to all sessions
% for i = 1:size(linkMat2,1)
%     for j=1:size(linkMat2,2)
%         if linkMat2(i,j)~=0
%             continue
%         else
%             num=i-1; % number of cells common to all sessions
%             break
%         end
%     end
%     if linkMat2(i,j)==0
%         break
%     end
% end
% %% collate normalized (df/F) stim values of all common rois
% x=0;
% for j=1:size(linkMat2,2) % index of sessions
%     if size(sig{1, j}.stat_arr,1)>x
%         x=size(sig{1, j}.stat_arr,1);
%     end
% end
% %%
% st_dF=nan(x,size(linkMat2,2),num);
% for i=1:num % no. of common rois
%     for j=1:size(linkMat2,2) % index of sessions
%         for k=1:size(sig{1, j}.stat_arr,1) % number of trials
%             st_dF(k,j,i)=sig{1, j}.stat_arr(k,2,linkMat2(i,j));
%         end
%     end
% end
% 
% %% replace zeros (not actual measurements) with nans
% % st_dF_nan=nan(size(st_dF,1),size(st_dF,2),size(st_dF,3));
% % for i=1:size(st_dF,3) % no. of responsive rois
% %     for j=1:size(st_dF,2) % index of sessions
% % %         for k=1:size(dF,1) % number of trials
% %             if st_pre(k,j,i)~=0 | isnan(st_pre(k,j,i))
% %                 st_pre_nan(k,j,i)=st_pre(k,j,i);
% %             end
% %         end
% %     end
% % end
% 
% %
% %
% %
% %
% %
% 
% 
% 
% %% select rois responsive in atleast one session
% sess_rois=zeros(num,1); % logical
% sess_rois_num=[]; % actual numbers only
% c=1; % count
% for i=1:num % no. of rois
%     for j=1:size(linkMat2,2) % number of sessions
%         if sig{1,j}.resp(1,linkMat2(i,j))==1
%             sess_rois(i,1)=1;
%             sess_rois_num(c,1)=i;
%             c=c+1;
%             break
%         end
%     end
% end
% 
% %% collate stim-pre values of resp rois
% x=0;
% for j=1:size(linkMat2,2) % index of sessions
%     if size(sig{1, j}.stat_arr,1)>x
%         x=j;
%     end
% end
% st_pre=nan(x,size(linkMat2,2),size(sess_rois_num,1));
% for i=1:size(sess_rois_num,1) % no. of responsive rois
%     for j=1:size(linkMat2,2) % index of sessions
%         for k=1:size(sig{1, j}.stat_arr,1) % number of trials
%             st_pre(k,j,i)=sig{1, j}.stat_arr(k,2,linkMat2(sess_rois_num(i,1),j))-sig{1, j}.stat_arr(k,1,linkMat2(sess_rois_num(i,1),j));
%         end
%     end
% end
% 
% % replace zeros (not actual measurements) with nans
% st_pre_nan=nan(size(st_pre,1),size(st_pre,2),size(st_pre,3));
% for i=1:size(st_pre,3) % no. of responsive rois
%     for j=1:size(st_pre,2) % index of sessions
%         for k=1:size(st_pre,1) % number of trials
%             if st_pre(k,j,i)~=0 | isnan(st_pre(k,j,i))
%                 st_pre_nan(k,j,i)=st_pre(k,j,i);
%             end
%         end
%     end
% end
% 
% %% compare rois shock observation > shock control
% sr.val={};
% sig_rois=zeros(1,size(st_pre_nan,3));
% sig_rois_ind=[];
% c=1;
% for i=1:size(st_pre_nan,3) % for each ROI
%     %     if isfinite(st_pre_nan(:,2,i)) && isfinite(st_pre_nan(:,1,i))
%     [p,h,stats] = signrank(st_pre_nan(:,2,i),st_pre_nan(:,1,i),'alpha',0.05,'tail','left');
%     sr.val.p(i)=p;
%     if sr.val.p(i)<0.05
%         sig_rois(1,i)=1;
%         sig_rois_ind(c,1)=i;
%         c=c+1;
%     end
%     sr.val.h(i)=h;
%     sr.val.stat(i)=stats;
%     %     end
% end
% 
% % check whether the shock observation > shock control rois are also
% % responsive for laser or squeak
% 
% sig_rois_ind_ori=[];
% for i=1:size(sig_rois_ind,1)
%     sig_rois_ind_ori(i,1)=linkMat2(sess_rois_num(sig_rois_ind(i,1),1),1);
%     if any(sig{1, 3}.resp_rois(1,:)==linkMat2(sess_rois_num(sig_rois_ind(i,1),1),3))
%         sig_rois_ind_ori(i,2)=linkMat2(sess_rois_num(sig_rois_ind(i,1),1),3);
%     end
%     if any(sig{1, 4}.resp_rois(1,:)==linkMat2(sess_rois_num(sig_rois_ind(i,1),1),4))
%         sig_rois_ind_ori(i,3)=linkMat2(sess_rois_num(sig_rois_ind(i,1),1),4);
%     end
% end

%% save values
cd(pn1)
save('comm_rois','sig_obs_rois','sig_obs_rois_ind')
clear
clc



