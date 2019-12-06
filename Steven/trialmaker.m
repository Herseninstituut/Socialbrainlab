function triallist = trialmaker(blocks,conditions) %% creates triallist 

  trials = zeros(blocks,conditions)'
  triallist = zeros(blocks*conditions, 1)
for j = 1:blocks
   
    trials(:,j)=randperm(conditions,conditions)
    
end

for k = 1:(blocks*conditions)
    
  triallist(k,1)= trials(k)
end