function [mhat] = find_ML_path(node,k)
%Function:  Computes the ML estimate by traversing the Trellis Map, looking
%           for the survivors.

%node - nodes corresponding to the trellis map, contains survivors and
%       state transition
%k    - length of message that we are seeking.

%Initialize survivor and cost list.
p_survivor = zeros(1,length(node)+1);%length(node)=k+m+1，为时间点的个数
%存储每个时间点幸存路径上的节点状态序号
cost       = zeros(1,length(node)+1);%路径度量

%Initialize branch output, 1 = top branch taken, 0 = lower branch
branch     = ones(1,length(node));%说明输入的是0还是1

%Initialize Survivor Trackback
branch(end)       = 0;
p_survivor(end)   = 1; %存储最后一个时间点（上一个）的第一个节点在幸存路径上
cost(end)         = node{length(node)}{1}.cost;%总路径度量
p_survivor(end-1) = node{length(node)}{1}.surv;%右式存储的是上一时刻幸存路径上的节点状态序号
cost(end-1)       = node{length(node)}{1}.cost;%

%Traverse Backwards -- look for surviving branches回溯
for n=length(node)-1:-1:1
    p_survivor(n) = node{n}{p_survivor(n+1)}.surv;
    cost(n) = node{n}{p_survivor(n+1)}.cost;
    
    %If we take the lower branch, assign a 0.  Otherwise top branch is 
    %assigned a 0.
    if(node{n}{p_survivor(n+1)}.f{1} == p_survivor(n+2))
        branch(n) = 0;
    end
end

%The code word is only the first kth bits.  The last length(branch)-k bits
%are by defination, 0.
mhat = branch(1:k);
    