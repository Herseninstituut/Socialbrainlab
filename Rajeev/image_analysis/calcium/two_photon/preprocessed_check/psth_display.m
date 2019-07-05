%% created by Rajeev Rajendran
% 2018_10_22
% load data from analysis folder
% psth_array(frame,repeat,angle,contrast,roi,prestim_or_stim)
% psth_array(i,j,k,m,n,p)
load('psth_array.mat')
load('resp_array.mat')
load('psth.details.stim.mat')
% array should have mean_value with dimensions of angles,time,roi,angles,contrast
psth_mean=nanmean(psth_array,2); % 6D(frame,repeat=1,angle,contrast,roi,prestim_or_stim)
cmap = jet(12);
roi=size(rois,1);
figure;
for g=1:roi
    for q=1:size(c,1)
        for p=1:size(a)
            pre_psth=psth_mean(:,1,p,q,g,1);
            post_psth=psth_mean(:,1,p,q,g,2);
            all_psth=[pre_psth;post_psth];
            all_psth2=all_psth;
            all_psth2(isnan(all_psth2))=[];
            psth(:,p)=all_psth2;
        end
        subplot(1,2,q)
        for p=1:size(a)
            plot(psth(:,p),'color',cmap(p,:))
            hold on
            legend ('0','30','60','90','120','150','180','210','240','270','300','330','Location', 'NorthWest')
            title(['ROI ',num2str(rois(g))])
        end
    end
    w = waitforbuttonpress;
    if isequal(w,0)
    end
    clf('reset')
end
save('psth.details.stim.mat','cmap')
