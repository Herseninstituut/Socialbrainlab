% Created by Rajeev Rajendran
% 2019 July 03
% split quad data for 3 z-stacks
% makes the analysis and split folders and saves the splits files in
% appropriate folders

function quad_split3

% select the correct quadrature file from the correct raw image folder
% containing the original sbx and quadrature files
[fn , pn] = uigetfile('*_quadrature.mat');
% change current directory to the folder
cd(pn)
% load the correct quadrature file
quad = load(uigetfile('*_quadrature.mat')); 
quad = quad.quad_data;
% make blank files for the splits
quad_split1=NaN(1,floor((size(quad,2))/3));
quad_split2=NaN(1,floor((size(quad,2))/3));
quad_split3=NaN(1,floor((size(quad,2))/3));
% transfer the correct run file entry to the correct split
for i=1:(floor(size(quad,2)))
    if floor((i-1)/3)==(i-1)/3
    quad_split1(1,floor(i/3)+1)= quad(1,i);
    elseif floor((i-2)/3)==(i-2)/3
    quad_split2(1,floor(i/3)+1)= quad(1,i);
    elseif floor((i-3)/3)==(i-3)/3
    quad_split3(1,floor(i/3)+1)= quad(1,i);
    end 
end
% check and/or make split directories
if exist([pn 'analysis\split1\'])  
else
    mkdir([pn 'analysis\split1\'])
end
if exist([pn 'analysis\split2\'])  
else
    mkdir([pn 'analysis\split2\'])
end
if exist([pn 'analysis\split3\'])  
else
    mkdir([pn 'analysis\split3\'])
end
% save split files to the correct split directories
save([pn 'analysis\split1\' fn(1:end-15) '_quad_split1.mat'],'quad_split1')
save([pn 'analysis\split2\' fn(1:end-15) '_quad_split2.mat'],'quad_split2')
save([pn 'analysis\split3\' fn(1:end-15) '_quad_split3.mat'],'quad_split3')

