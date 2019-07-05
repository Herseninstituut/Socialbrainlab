% Created by Rajeev Rajendran
% 2019 July 03
% split quad data for 3 z-stacks
function quad_split3
[fn , pn] = uigetfile('*_quadrature.mat'); % select the correct quadrature file
cd(pn)
quad = load(uigetfile('*_quadrature.mat')); % select the correct quadrature file again
quad = quad.quad_data;
quad_split1=NaN(1,floor((size(quad,2))/3));
quad_split2=NaN(1,floor((size(quad,2))/3));
quad_split3=NaN(1,floor((size(quad,2))/3));

for i=1:(floor(size(quad,2)))
    if floor((i-1)/3)==(i-1)/3
    quad_split1(1,floor(i/3)+1)= quad(1,i);
    elseif floor((i-2)/3)==(i-2)/3
    quad_split2(1,floor(i/3)+1)= quad(1,i);
    elseif floor((i-3)/3)==(i-3)/3
    quad_split3(1,floor(i/3)+1)= quad(1,i);
    end 
end
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
save([pn 'analysis\split1\' fn(1:end-15) '_quad_split1.mat'],'quad_split1')
save([pn 'analysis\split2\' fn(1:end-15) '_quad_split2.mat'],'quad_split2')
save([pn 'analysis\split3\' fn(1:end-15) '_quad_split3.mat'],'quad_split3')

