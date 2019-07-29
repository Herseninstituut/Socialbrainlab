%Ca image processing pipeline as of Octber 2017
%This document is under constant revision as the scripts we use are also
%changing.
%save a local copy of this file to your local matlab path so you can change
%whatever you like without changing the GIT version

%start parallel pool if neccessary
%g = gcp;

%use dbBrowser to select datasets, and view metadata
% dbBrowser

%change working directory
pn = uigetdir();
cd(pn)

%have a look at your sbx file, 
%split data in separate slices if neccessary
%average over selected number of trials, at a selected position and save the image
%remove a frame with an artifact

Showsbx % split data

%concatenate sbx files that were recorded consecutively at the same position
%to do registration and roi segmentation over several stacks simultaneously
%sbxconcatenate


%% using simons foundation; normcorre  alternative for (python) CaIman from Simons foundation
% simonalign  %or 
simonsplitalign % select all splits
% after splitting move files (except original data) to analysis folder and 
% each split in a different folder

%% Directly retrieve Rois for the red channel (structural based, no cnmf or spectral processing)
% getroisfromredchannel

%% functions for roi selection and signal extraction %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% spectral correlation of adjacent pixels and roi segmentation, and deconvolve signals
Spectralprocess_rr %please open this file and run step by step ; gets rois and signals

%% or do ca source extration %please open first, and run step by step
% NMFmemmap %gets rois and signals

%% Chronic analysis
% Match ROIs of multiple recordings together
chronic_matching_multiple

%% Displayrois(Mask, PP, BImg): show rois, input: mask, contours, background image (input is the outcome of
% Spectralprocess or NMFmemmap, compare with other datasets)
Displayrois()
%Check_overlap : to see outcome of the GLTbl 

%% extract stimulus ordered trials, save with all info to '*_RES.mat' file
% extract_response();
% 
% extract();

%% compare behaviour with signals and stimuli
% plot_Res()

%% do reverse correlation for receptive field plotting, 
% reversecor


%% orientation and direction tuning (Koen) for Grating, Cirgratsound, Movingbar
% oridirtuning()

%% extract speed from quad data
% file in Socialbrainlab\Rajeev\image_analysis\calcium\two_photon\running
% file can be run directly
quad_split3 % for 3 z-slices only

%% check for responsive cells
% file in Socialbrainlab\Rajeev\image_analysis\calcium\two_photon\response
% just run the file 
% inputs needed are 
% 1. stim_first = index of first stim
% 2. stim_last = index of last stim
% 3. pre_stim_dur = duration of pre_stim for F calculation, in seconds
% 4. stim_dur = duration of stim, in seconds
% 5. post_stim_dur = duration of post_stim for dF calculation, in seconds

responsive_cells

% also same inputs are needed for randomly active cells

responsive_cells_rand

%% collate and compare responsive cells across sessions
% file in Socialbrainlab\Rajeev\image_analysis\calcium\two_photon\multi_session 
% just run the file 
% inputs needed are 
% 1. parent folder of all the sessions 
% 2. specific split subfolder in mult_sessions folder
% 3. load SPSIG file from the following sessions in order
%       shock_observation
%       shock_control
%       laser_self or laser_self_select
%       squeak_playback

cell_mult_sess_coll

%% for regression of run speed related activity
RunGLMFitCaResFullNeg
responsive_cells_glm
cell_mult_sess_coll_glm



%% Align your image stack using scanbox routines 
% 
%  %scanbox version of aligning
% [fn, pn] = uigetfile(' *.sbx', 'Get Sbx file', 'MultiSelect', 'on');
% cd(pn);
% if iscell(fn)
%     sbxaligndir(fn)
% else
%     sbxaligndir({fn})
% end
% 
% sbxaligntool
% 
%  % Segment rois with scanbox tool
% sbxsegmenttool


%% locate and open segmentation file for mask data after using dario's segmentation
    % [fn pn] = uigetfile(' *.segment', 'segment');
    % Rm = matfile([pn fn]);
    % mask = Rm.mask;
    % figure, imagesc(mask)
    % rois = unique(mask);
    % rois(1) = [];
    % 
    % %get median positions for each roi
    % for i = 1:length(rois)
    %     [iy,ix] = find(mask == rois(i));
    %     mx = median(ix);
    %     my = median(iy);
    %     text(mx, my, num2str(i), 'color', 'w')
    %     Rois(i).mx = mx;
    %     Rois(i).my = my;
    %     Rois(i).msdelay = Tline * my * 1000;
    % end


%% locate and load signals file (imaging PC)
    % [fn pn] = uigetfile(' *.signals', 'signals');
    % Rm = matfile([pn fn]);
    % sig = Rm.sig;
    % sig = sig(:,rois);
    % figure, plot(sig)






%% DONE rest is obsolete or examples 

%% Subtract baseline

%  Yp = bsxfun(@minus, sig, median(sig));
%  Yp = double(Yp./(nanmedian(max(Yp))));


%%   Deconvolve the signals with Dario's method to spike records (see scanbox blog: extremely simple deconvolution algorithm )

% shat = [];
% Est = [];
% a = lpc(Yp,10); % order=10
% for i = 1:size(Yp,2)
%     Est(:,i) = filter([0 -a(i,2:end)],1,Yp(:,i));
%     uhat = Yp(:,i)-Est(:,i); % estimated spike train
%     shat(:,i) = uhat > graythresh(uhat); 
% end


%%

% figure, 
% for i = 1:size(Yp,2)
%     hold off
%     %plot(FrameTimes + Rois(i).msdelay/1000, Yp(:,i).*10)
%     plot(FrameTimes + Rois(i).msdelay/1000, Yp(:,i))
%     hold on, plot(StimTimes, -0.1, '+')
%     %plot(FrameTimes + Rois(i).msdelay/1000, shat(:,i))
%     plot(FrameTimes + Rois(i).msdelay/1000, Est(:,i), 'g')
% 
%     spike = shat(:,i) > 0; 
%     plot(FrameTimes(spike) + Rois(i).msdelay/1000, shat(spike,i)*0.1-0.2, '.')
% 
%    % plot(RunTimes, RunSpeed - 1, 'r')
%     title(num2str(rois(i)))
%     pause
% end


%% transform data in piecewise bicubic spline fits, taking the delay of each
%roi into account
% P.pp = [];
% for i = 1:size(sig,2)
%     P(i).pp = spline(FrameTimes + Rois(i).msdelay/1000,Yp(:,i));
% end


%% example averaging over stimuli 
% [fn pn] = uigetfile(' *.mat', 'Get Log');
% load([pn fn]);
% 
% 
% Angles = unique(log(:,1));
% win = [-1.0 2];
% %we can now define a sampling rate to subsample from the spline fits; say
% %at 20 hrz
% xx = (win(1):0.05:win(2));
% numofsamples = length(xx);
% numofRois = size(sig,2);
% numofreps = length(log)/length(Angles);
% 
% %samples, numoftrials, rois, angles
% Data = zeros(numofsamples, numofreps, numofRois, length(Angles));
% for i = 1:length(Angles)
%     Adx = find(log == Angles(i));
%     for j = 1:length(Adx)
%         Stimonset = StimTimes(Adx(j));
%         for k = 1:numofRois
%             Data(:,j,k,i) = ppval(P(k).pp,xx+Stimonset);
%         end      
%     end
% end



%%  testing and plotting the outcome

%Average responses over all rois
% figure, 
% for i = 1:length(Angles)
%     Res = mean(mean(Data(:,:,:,i),3),2);
%     SE = mean(squeeze(std(Data(:,:,:,i), [] , 2)),2)/sqrt(size(Data(:,:,:,i),2));
%     plot(xx, Res);
%     hold on, plot(xx, Res+mean(SE,2), 'm')
%     title(num2str(Angles(i)))
%     pause
%     hold off
% end
% 
% %responses for each ROI separately
% h = figure; 
% for j = 1:numofRois
%     set(h, 'name', ['Roi: ' num2str(j)])
%     Res = squeeze(mean(Data(:,:,j,:),2));
%     SE = squeeze(std(Data(:,:,j,:), [] , 2)/sqrt(size(Data(:,:,j,:),2)));
%     for k = 1:length(Angles)
%         subplot(2,length(Angles)/2,k)
%         plot(xx, Res(:,k))
%         hold on, plot(xx, Res(:,k)+SE(:,k), 'm')
%         axis tight
%         hold off
%     end
%     pause 
% end
% 
