function [mhat,node] = decode_1_3(r,n,mem,k,flag)
%Function: This is the decoder for a generalized convolution encoder.
%This file is independent of the desired rate and memory elements.  It
%first produces a trellis map where we have assigned node states (previous
%and forward) as well as the cost functions associated with received vector
%and acceptable code words.  Note:  Only implements for our binary case.

%r    - codeword接收到的码元
%n    - output (1/n) convolution coder
%mem  - number of memory elements卷积码存储
%k    - number of original message bits信息码元个数
%flag - Soft/Hard Decoding ==> 0/1.

%%%%%%%%%%%%%%%%%%%%%%%%%%%PSEUDO CODE%%%=%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%1) Initialize next set of nodes, states, stages, etc.

%2) Traverse through nodes that have only been visited (i.e., no need to
%   check on nodes (2,3,4) of a 4 state Trellis Map at time = 1 since we
%   know that we should only visit node 1 given that we begin here.初始化第一个
%节点，所以只遍历第一个节点（访问过的）

%3) Given q = 2, we have two inputs = {0,1}.  Input each value and update
%   nodes accordingly 
%   
%4) Compute acceptable output for each branch of the map, which is
%   formed from the function circuit_logic.  This changes
%   with each encoder.

%5) Compute distance between received vector and acceptable vector for each
%   of the branch at the ith stage.

%6) Update Nodes - a) Next State of the node (Could have 2 Possibilities) 
%                    b) Previous State of Node (Could have 2 Possibilities)
%                    c )Node has been visited?  
%                    d)Total Cost assigned to Node
%                    e) List/Determine Surviving Branch
%
%7) Form the decoded message by traversing backwards and finding the 
%   surviving Branches and Maximum Likelihood (ML) estimate.  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%(STEP 1)%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Develop Trellis Map for decoder.  Find # of stages and # of states.
stages = k+mem;%状态转移的次数（0-stages,时间节点），k为信息位的个数
states = 2^mem; %每个节点出发有2^k(k=1)条分支，每个分支有n个数据，最多可能有2^(km)种不同的状态
block_st= 1;%分支对应的接收序列的序号（1bit变成n=3bit,每次转移计算分支度量都要加n）
            %即，块状态
%------

%If flag = 1 -> Hard Decoding.  We must first make "hard" decisions of the
%input received vector.
if(flag)
    ind_1 = r>0;
    ind_0 = r<=0;
    r(ind_1) = 1;  %接收序列中大于0的改为1
    r(ind_0) = -1; %接受序列中小于等于0的改为-1
end
%-----

%Initialize State or Memory/Input Elements初始化
for i = 1:mem; c_S.m{i} = 0; end   %表示该时刻节点的状态序号初始化为000（0-2^m-1）
c_S.st    = 1;%表示该时刻节点的状态序号（可由.m求得）
c_S.in    = 0; %使状态转移的输入
%-----

%Initialize Trellis Map, which state we start with etc.
for l = 1:states
    %初始化第一个时间节点的所有状态
    node{1}{l}.p{1}   = NaN;%与之相连的上一时刻节点的状态序号（从1开始）
    node{1}{l}.p{2}   = NaN;%与之相连的上一时刻节点的状态序号（从1开始），有两条路径同时到达i+1时刻
    node{1}{l}.f{1}   = NaN;%存储输入为0的下一时刻节点的状态序号
    node{1}{l}.f{2}   = NaN;%存储输入为1的下一时刻节点的状态序号
    node{1}{l}.cost   = -100000;%存储至当前时刻的路径度量
    node{1}{l}.visit  = 0;%表示是否访问过，0表示没有遍历（不是路径上的节点）
    node{1}{l}.surv   = NaN;%存储幸存路径上的上一时刻节点的状态序号
end
node{1}{1}.visit = 1; %初始化从第一个时刻的第一个状态开始状态转移
node{1}{1}.surv  = 1;  
node{1}{1}.cost  = 0;  %路径度量初始化为0
%-----

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%(STEP 2)%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i = 1:stages
    %Update block status（块状态初始化为1）, and initialize next set of nodes
    if(i~=1);block_st = block_st+n; end
    for l = 1:states                    %从第二个转移时间节点开始初始化
        node{i+1}{l}.p{1}   = NaN; %2-stages 
        node{i+1}{l}.p{2}   = NaN;
        node{i+1}{l}.f{1}   = NaN; 
        node{i+1}{l}.f{2}   = NaN;
        node{i+1}{l}.cost   = -100000;
        node{i+1}{l}.visit  = 0; 
        node{i+1}{l}.surv   = NaN;
    end
    
    %For each state or node, check if we need to do processing on.
    for l = 1:states
        if(node{i}{l}.visit)  %如果被遍历
            %If we do process this node, determine its current numerical
            %state
            c_S.st = l; %当前的状态序号
            %将十进制转化成二进制（二进制从000开始表示，所以十进制先减1）
            val = l-1;  
            for j = mem-1:-1:0
                if((val - 2^j)>=0)
                    c_S.m{j+1} = 1;
                    val = val-2^j;
                else
                    c_S.m{j+1} = 0;
                end
            end
            
            %%%%%%%%%%%%%%%%%%%%(STEP 3-5)%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %State Input = 0; (Binary, q=2)
            c_S.in   = 0;%对应篱笆图，输入0走上面分支，输入为1走下面分支，恩恩
            
            %Determine acceptable output from circuit logic
            [o,n_S] = circuit_logic(c_S,n,mem);%输出相对应的数据和下一时刻到达的节点状态
            
            %Soft or Hard Decoding (e.g., Hard => use Hamming Distance)
            %计算分支度量，转移一个时间节点
            if(flag)
                dist    = compute_Hamm(o,r,block_st,n);
            else
                dist    = compute_Lp(o,r,block_st,n);
            end
            
            %%%%%%%%%%%%%%%%%%%%(STEP 6)%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %Update node's status (e.g., node's status, total cost, is it a
            %possible survivor?)
            node{i}{c_S.st}.f{1}  = n_S.st;%存储输入为0时下一时刻节点的状态序号
            node{i}{c_S.st}.visit = 1;
           
            if(isnan(node{i+1}{n_S.st}.p{1}))  %如果p{1}存储的是NAN（初始化为NAN）
                node{i+1}{n_S.st}.p{1}  = c_S.st;%存储上一时刻节点的状态序号
                node{i+1}{n_S.st}.visit = 1;
                node{i+1}{n_S.st}.cost  = node{i}{c_S.st}.cost+dist; %计算路径度量
                node{i+1}{n_S.st}.surv  = c_S.st;%上一时刻节点的状态序号为幸存路径的最后一个节点
            else      %如果i+1时刻到达该节点的已经有上一节点信息（i时刻2个节点状态转移到i+1时刻节点）
                node{i+1}{n_S.st}.p{2} = c_S.st;
                node{i+1}{n_S.st}.visit = 1;
                
                %Two Possible Survivors: Determine surviving branch
                if(node{i+1}{n_S.st}.cost<=node{i}{c_S.st}.cost+dist)%已有的节点到i+1时刻节点的路径度量小
                    node{i+1}{n_S.st}.surv  =node{i+1}{n_S.st}.p{1};
                    %幸存路径存储已有的节点
                else
                    node{i+1}{n_S.st}.surv  = c_S.st;
                    node{i+1}{n_S.st}.cost  = node{i}{c_S.st}.cost+dist; 
                end       
            end
            
            %%%%%%%%%%%%%%%%%%%%(STEP 3-5)%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %State Input = 1; (Binary, q=2，软判决中，二电平量化)
            if(i<=k)
                %Only process input 1 for first k stages
                c_S.in   = 1; 
                
                 %Determine acceptable output from circuit logic
                [o,n_S] = circuit_logic(c_S,n,mem);%计算出下一个时刻输出数据（n位）
                                                   %以及节点状态
                
                %Update node's status (e.g., node's status, total cost, is it a
                %possible survivor?)
                if(flag)
                    dist    = compute_Hamm(o,r,block_st,n);
                else
                    dist    = compute_Lp(o,r,block_st,n);  %计算分支度量
                end;
                 
                %%%%%%%%%%%%%%%%%%%%(STEP 6)%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                node{i}{c_S.st}.f{2}  = n_S.st;%存储输入为1时下一节点的状态序号
                node{i}{c_S.st}.visit = 1;
                if(isnan(node{i+1}{n_S.st}.p{1}))   
                    node{i+1}{n_S.st}.p{1} = c_S.st;
                    node{i+1}{n_S.st}.visit = 1;
                    node{i+1}{n_S.st}.cost  = node{i}{c_S.st}.cost+dist;
                    node{i+1}{n_S.st}.surv  = c_S.st;
                else                                 
                    node{i+1}{n_S.st}.p{2} = c_S.st;
                    node{i+1}{n_S.st}.visit = 1;
                    
                    %Two Possible Survivors: Determine surviving branch
                    if(node{i+1}{n_S.st}.cost<=node{i}{c_S.st}.cost+dist)
                        node{i+1}{n_S.st}.surv  =node{i+1}{n_S.st}.p{1};
                    else
                        node{i+1}{n_S.st}.surv  = c_S.st;
                        node{i+1}{n_S.st}.cost  = node{i}{c_S.st}.cost+dist; 
                    end
                end
            end
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%(STEP 7)%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mhat = find_ML_path(node,k);
end