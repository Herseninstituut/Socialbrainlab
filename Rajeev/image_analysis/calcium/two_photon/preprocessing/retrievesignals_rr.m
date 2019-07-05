%retrieve signals using rois
% 
% %first get rois from spectral components
% [fn, pn] = uigetfile('*_ROI.mat');
% filename = [pn fn];
% load(filename)
% 
% %Can be skipped : to display and find overlapping rois plot mask or image
% %from one session and rois(PP) from other session
% %load spectral images and take lowest spectral component
% %  load([filename(1:end-8) '.mat'])
% %  BImg = log(SPic(:,:,2));
% %  clear SPic
% % Sax(1) = [];
% [fn, pn] = uigetfile('*.tif');
% BImg = imread([pn fn]);
% 
% %show rois, mask, contours, 1st spectral image
% Displayrois(Mask, PP, BImg)

%% get calcium data associated with rois from transposed block 
% there are small differences with the rois taken from the image stack!!!
  
   [fn, fp] = uigetfile('*Trans*.dat');
   if fn ~= 0 %if not cancelled
        Sigfn = [fp fn(1:end-10) '_SPSIG.mat'];
        load(Sigfn, 'Mask', 'PP');
        %load transposed data
        %[dim, freq] = Zgetinfo([fp fn]);
        %mfn = strsplit(fn, '.');

        %sbxt = memmapfile([fp fn], 'Format', 'uint16', 'Offset', 500);
        [sbxt, dim, freq] = transmemap([fp fn]);

        Mt = Mask'; %the mask needs to be transposed
        indices = find(Mt(:)> 0); %get the image indices for all rois
        %indices = indices - double(dim(2)); %add a line?! 

        %Sigrois = XYgetZ(indices, sbxt, dim(1));
        Sigrois = sbxt.Data.y(:,indices);   
        Sigrois = double(Sigrois);
        %sig = mean(Sigrois,2);

        sig = zeros(dim(1),PP.Cnt);
        for i = 1:PP.Cnt
            ichan = find(Mt(:)== i);
            [~, ia] = intersect(indices, ichan); %get indices corresponding to this roi
            sig(:,i) = mean(Sigrois(:,ia),2); %and average this selection 
        end

        xas = (1:length(sig))./freq;

        window            = round(5000*freq/15); %adapted for framerate
        %percentile        = 10;  % default

        % baseline estimation due to Pnevmatikakis et al
        sigraw = sig;
        sig = basecorrect(sig, window);
        figure, plot(xas, sig)

        save(Sigfn, 'sig','sigraw', 'freq', '-append')
   end

%% Deconvolving using MLspike


% parameters
parallel_process = true; % default
%load('MLspike_paramvals_PLinterneurons', 'MLspike_params'); % default for Prem
% MLspike_params = load('chrisMLpest', 'MLpest'); %these are better parameters, but pleas check with G = spk_demoGUI
% MLspike_params = load('MLpest20180927', 'MLpest'); %these are better parameters, but pleas check with G = spk_demoGUI
MLspike_params = load('MLpest2018rr', 'MLpest');
MLspike_params = MLspike_params.MLpest;
MLspike_params.dt = 1/freq;

nrois             = size(sig,2);
t = (1:length(sig)+1)/freq;

drifting_baseline = true; % default
if ~drifting_baseline
    % prevent drifting baseline estimation during MLspike deconvolution
    MLspike_params = rmfield(MLspike_params,'drift');
    dbsl_process = 'fixed baseline';
else
    dbsl_process = 'drifting baseline';
end

MLpar = repmat(MLspike_params, nrois, 1);

% data matrix preallocation
decon  = zeros(size(sig));
den    = zeros(size(sig));
spikes = cell(nrois,1);
% MLspike deconvolution
tic
if parallel_process
    if ~verLessThan('matlab','9.3')
        % define waitbar
        hWb = waitbar(0,sprintf('Running MLspike (Runtime %1.2f min)',0.00),...
                                'Name','Deconvolution');
        global h N p
        D = parallel.pool.DataQueue;
        afterEach(D, @retrievesignals_parWB);
        h = hWb;
        N = nrois;
        p = 1;
        parfor roi = 1:N % parallel processing should be possible, only the waitbar won't make much sense.
            [spikes{roi},den(:,roi),~] = spk_est(sig(:,roi), MLpar(roi));
            send(D,roi)
        end
        close(hWb)
        
    else % older versions of matlab (<2017b)
        
        N = nrois;
        parfor roi = 1:N % parallel processing should be possible, only the waitbar won't make much sense.           
            [spikes{roi},den(:,roi),~] = spk_est(sig(:,roi),MLpar(roi));
        end
    end
end
fprintf('MLspike took %1.2f minutes to process %d rois (%s).\n',toc/60,...
        nrois,dbsl_process)


parfor roi = 1:nrois % parallel processing should be possible, yields little.
   % decon(:,roi) = sum(spikes{roi}'>=t(1:end-1)&spikes{roi}'<t(2:end))';
   if ~isempty(spikes{roi})
    out = arrayfun(@(x) x>=t(1:end-1)&x<t(2:end), spikes{roi}', 'UniformOutput', false);
    decon(:,roi) = sum(cell2mat(out))';
   end
end

save(Sigfn, 'den', 'decon', '-append')

%% Deconvolving  [c, s, options] = deconvolveCa(y, varargin) (from ca source extraction also download cvx and setup in matlab)

%[fn, fp] = uigetfile('*_Sig.mat');
%load([fp fn])
%
% gp = gitpath();
% addpath( [gp '\ca_source_extraction'])
% addpath( [gp '\ca_source_extraction\utilities'])
% addpath( [gp '\ca_source_extraction\deconvolution'])
% addpath( [gp '\ca_source_extraction\deconvolution\functions'])

% parameter settings for OASIS deconvolution
% fit_type = 'ar2';
% pars = [];
% sn = [];
% b = 0;
% lambda = 0;
% optimize_b = true;
% optimize_pars = true;
% optimize_smin = true;
% method = 'thresholded';
% window = 200;
% shift = 100;
% smin = 0.0015; % smin seems to work: anecdotal evidence
% maxIter = 20;
% thresh_factor = 1.0;
% 
% den = zeros(size(sig));
% decon = zeros(size(sig));
% bsl = zeros(size(sig,2),1);
% siglength = size(sig,1);
% 
% for i = 1:size(sig,2)
%     try
%         [den(:,i), decon(:,i), optns] = deconvolveCa(sig(:,i)-1,fit_type,method,...
%         'pars',pars,'sn',sn,'b',b,'lambda',lambda,'optimize_b',optimize_b,...
%         'optimize_pars',optimize_pars,'optimize_smin',optimize_smin,...
%         'window',window,'shift',shift,'smin',smin,'maxIter',maxIter,...
%         'thresh_factor',thresh_factor);
%     catch ME
%         warning(['Error occured for ROI: ' num2str(i) '; ' ME.identifier])
%     end
% end
% 
% save(Sigfn,'den', 'decon', '-append')

%%
% Yp = sig(:,1) - baseline.b;
% %dario method
% shat = [];
% Est = [];
% a = lpc(Yp,2); % order=10
% for i = 1:size(Yp,2)
%     Est(:,i) = filter([0 -a(i,2:end)],1,Yp(:,i));
%     uhat = Yp(:,i)-Est(:,i); % estimated spike train
%     shat(:,i) = uhat > graythresh(uhat); 
% end
% 
% figure, plot(Yp, 'g')
% hold on, plot(Est, 'b')
% hold on, plot(shat*100, 'k')



% %% you can also load the rois to imagej roimanager, and use the original image sequence to get your data
% % this is a nice check to see that the previous import was correct 
% % read the actual data from the normcorre aligned sbx file and make virtual stack in ImageJ
%      sbx2ij
%      
%       %load these pnts to imagej roi manager
%      rois = PP.Con(:);
%      rois2ij('POLYGON', rois)
%  
%      %multimeasure rois in Imagej and get results table
%      %There are so much rois that normal MIJ.getResultsTable() does not work
%      %properly
%      RT = ij.measure.ResultsTable.getResultsTable();
%      lc = RT.getLastColumn();
%      %get all the columns separately in matlab arrays
%      sig = RT.getColumn(0);
%      for i = 1:lc
%          sig(:,i+1) = RT.getColumn(i);
%      end
%   
%      
%      save([pn fn(1:end-4) '_Sig.mat'], 'sig');
%      
%      %% Some other unneccessay stuff
%      %get the noise level of the roi traces
%      DF = median(abs(diff(R)));
%      %smooth and get maxima of ROI traces
%      Rs = smoothG(R,2);
%      MF = max(Rs);
%           
%      %find responses 25 x larger than noise level and the noise level should
%      %be at least 5 to avoid selecting dead pixels
%      Sel = find(MF./DF > 15 & DF > 5);
%      
% %      for i = 1:10:length(Sel)
% %          figure(2)
% %          for j = 1:10
% %              if (i-1+j) >length(Sel)
% %                  break;
% %              end
% %              subplot(10,1,j)
% %              plot(Rs(:,Sel(i-1+j)))
% %              hold on,
% %              line([1 length(R)],[DF(Sel(i-1+j))*3 DF(Sel(i-1+j))*3], 'color', 'red')
% %              hold off
% %              title(num2str(Sel(i-1+j)))
% %          end
% %          pause
% %      end
% %selected smoothed responses, and rois
%     R = Rs(:,Sel);
%     pnts = pnts(Sel);
% 
%  %% Lets have a look at the selected responses, tagged with their main spectral range
%      xas = (1:size(R,1))./freq;
% 
%      for l = 1:10:size(R,2)
%          figure(2)
%          ix = find(lns > Sel(l), 1, 'first');
%          str = {['Spectral component :' strfreq{ix}]};
%          ohndl = annotation('textbox', [0.15 0.95 0.3 0.03], 'String', str, 'FitBoxToText','off');
%          for j = 1:10
%              if (l-1+j) > size(R,2)
%                  break;
%              end
%              subplot(10,1,j)
%              plot(xas, R(:,l-1+j))
%              title(num2str(l-1+j))
% 
%          end
%          pause
%          delete(ohndl)
%      end
%     
% %% Now order to power in the high frequency range
%      xas = (1:size(R,1))./freq;
% 
%      for l = 1:5:size(R,2)
%          figure(2)
%          Si = ISx(l:l+4);
%          for j = 1:5
%             subplot(5,1,j)            
%             plot(xas, R(:,Si(j)))
%             title([num2str(Si(j)) ':' num2str(MS(Si(j))) ] )
%          end
%          pause
%      end
% 
%  
% %% which boutons are correlated
%  
% %simple correlations between pairs of responses
%  Nr = size(R,2);
%  C = cell(Nr*(Nr-1)/2+1,3);
%  cnt = 1;
%  for i = 1:size(R,2) %number of channels
%      x = R(:,i);
%      X = gpuArray(x);
%      for j = i+1:size(R,2);
%          y = R(:,j);
%          Y = gpuArray(y);
%          cnt = cnt + 1;
%            % tmp = corrcoef(x,y);
%            %C(cnt,1) = tmp(2,1);
%            [xc,lags] = xcov(X,Y,300, 'coeff');
%             xc = gather(xc);
%             C{cnt,1} = xc;
%             C{cnt,2} = i;
%             C{cnt,3} = j;
%           %  figure(2), plot(lags, xc)
%           %  pause
%      end
%  end
%  
%  %save C:\Users\togt\Desktop\prem_Hulk\transpose\Corr.mat C lags freq
%  
%  %significant correlations > 0.8
%  Mxdist = cellfun(@max,C(:,1));
%  Sc = find(Mxdist>0.6);
%  
%  %significant correlations > 0.8
%  Cs = C(C(:,1)>0.8,:);
%  pts = unique(Cs(:,2));
%  
%  figure(2)
%  for i = 1:length(Sc)
%      ix = Sc(i);
%      plot(lags./freq, C{ix,1})
%      title([num2str(C{ix,2}) ' - ' num2str(C{ix,3}) ])
%      pause
%  end
%  
%  %group rois with high correlations
%  grps = {};
%  cnt = 1;
%  bfound = false;
%  g = pts(1);
%  cp = find(Cs(:,2)== pts(1));
%  grps{1} = [g; Cs(cp,3)];
%  
%  for i = 2:length(pts)
%      g = pts(i);
%      cp = find(Cs(:,2)== g);
%      bfound = false;
%      for j = 1:cnt
%          pp = grps{j};
%          if ~isempty(find(pp == g, 1))
%              pp = [pp ; Cs(cp,3)];
%              grps{j} = pp;
%              bfound = true;
%              break;
%          end
%      end
%      if ~bfound
%         cnt = cnt +1;
%         grps{cnt} = [g; Cs(cp,3)];
%      end
%  end
%  
%  
%   %get image to plot correlated boutons
%  IMG = mean(Imgav(:,:,:),3);
%  
%  %show groups sharing high correlations
%  cols = {'r', 'g', 'b', 'y', 'c', 'm'};
%  for i =1:length(grps)
%      figure(1), imagesc(IMG), colormap gray , hold on   
%      pp = unique(grps{i});
%      %c = cols{mod(i,6)+1};
%      p1 = pnts(pp);
%      for j = 1:length(p1)
%         fill(p1(j).x, p1(j).y,  cols{4})
%      end
% 
%      hold off
%      pause
%  end
%              
% 
%  