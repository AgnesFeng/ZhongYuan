function [mhat,node] = decode_1_3(r,n,mem,k,flag)
%Function: This is the decoder for a generalized convolution encoder.
%This file is independent of the desired rate and memory elements.  It
%first produces a trellis map where we have assigned node states (previous
%and forward) as well as the cost functions associated with received vector
%and acceptable code words.  Note:  Only implements for our binary case.

%r    - codeword���յ�����Ԫ
%n    - output (1/n) convolution coder
%mem  - number of memory elements�����洢
%k    - number of original message bits��Ϣ��Ԫ����
%flag - Soft/Hard Decoding ==> 0/1.

%%%%%%%%%%%%%%%%%%%%%%%%%%%PSEUDO CODE%%%=%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%1) Initialize next set of nodes, states, stages, etc.

%2) Traverse through nodes that have only been visited (i.e., no need to
%   check on nodes (2,3,4) of a 4 state Trellis Map at time = 1 since we
%   know that we should only visit node 1 given that we begin here.��ʼ����һ��
%�ڵ㣬����ֻ������һ���ڵ㣨���ʹ��ģ�

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
stages = k+mem;%״̬ת�ƵĴ�����0-stages,ʱ��ڵ㣩��kΪ��Ϣλ�ĸ���
states = 2^mem; %ÿ���ڵ������2^k(k=1)����֧��ÿ����֧��n�����ݣ���������2^(km)�ֲ�ͬ��״̬
block_st= 1;%��֧��Ӧ�Ľ������е���ţ�1bit���n=3bit,ÿ��ת�Ƽ����֧������Ҫ��n��
            %������״̬
%------

%If flag = 1 -> Hard Decoding.  We must first make "hard" decisions of the
%input received vector.
if(flag)
    ind_1 = r>0;
    ind_0 = r<=0;
    r(ind_1) = 1;  %���������д���0�ĸ�Ϊ1
    r(ind_0) = -1; %����������С�ڵ���0�ĸ�Ϊ-1
end
%-----

%Initialize State or Memory/Input Elements��ʼ��
for i = 1:mem; c_S.m{i} = 0; end   %��ʾ��ʱ�̽ڵ��״̬��ų�ʼ��Ϊ000��0-2^m-1��
c_S.st    = 1;%��ʾ��ʱ�̽ڵ��״̬��ţ�����.m��ã�
c_S.in    = 0; %ʹ״̬ת�Ƶ�����
%-----

%Initialize Trellis Map, which state we start with etc.
for l = 1:states
    %��ʼ����һ��ʱ��ڵ������״̬
    node{1}{l}.p{1}   = NaN;%��֮��������һʱ�̽ڵ��״̬��ţ���1��ʼ��
    node{1}{l}.p{2}   = NaN;%��֮��������һʱ�̽ڵ��״̬��ţ���1��ʼ����������·��ͬʱ����i+1ʱ��
    node{1}{l}.f{1}   = NaN;%�洢����Ϊ0����һʱ�̽ڵ��״̬���
    node{1}{l}.f{2}   = NaN;%�洢����Ϊ1����һʱ�̽ڵ��״̬���
    node{1}{l}.cost   = -100000;%�洢����ǰʱ�̵�·������
    node{1}{l}.visit  = 0;%��ʾ�Ƿ���ʹ���0��ʾû�б���������·���ϵĽڵ㣩
    node{1}{l}.surv   = NaN;%�洢�Ҵ�·���ϵ���һʱ�̽ڵ��״̬���
end
node{1}{1}.visit = 1; %��ʼ���ӵ�һ��ʱ�̵ĵ�һ��״̬��ʼ״̬ת��
node{1}{1}.surv  = 1;  
node{1}{1}.cost  = 0;  %·��������ʼ��Ϊ0
%-----

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%(STEP 2)%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i = 1:stages
    %Update block status����״̬��ʼ��Ϊ1��, and initialize next set of nodes
    if(i~=1);block_st = block_st+n; end
    for l = 1:states                    %�ӵڶ���ת��ʱ��ڵ㿪ʼ��ʼ��
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
        if(node{i}{l}.visit)  %���������
            %If we do process this node, determine its current numerical
            %state
            c_S.st = l; %��ǰ��״̬���
            %��ʮ����ת���ɶ����ƣ������ƴ�000��ʼ��ʾ������ʮ�����ȼ�1��
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
            c_S.in   = 0;%��Ӧ���ͼ������0�������֧������Ϊ1�������֧������
            
            %Determine acceptable output from circuit logic
            [o,n_S] = circuit_logic(c_S,n,mem);%������Ӧ�����ݺ���һʱ�̵���Ľڵ�״̬
            
            %Soft or Hard Decoding (e.g., Hard => use Hamming Distance)
            %�����֧������ת��һ��ʱ��ڵ�
            if(flag)
                dist    = compute_Hamm(o,r,block_st,n);
            else
                dist    = compute_Lp(o,r,block_st,n);
            end
            
            %%%%%%%%%%%%%%%%%%%%(STEP 6)%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %Update node's status (e.g., node's status, total cost, is it a
            %possible survivor?)
            node{i}{c_S.st}.f{1}  = n_S.st;%�洢����Ϊ0ʱ��һʱ�̽ڵ��״̬���
            node{i}{c_S.st}.visit = 1;
           
            if(isnan(node{i+1}{n_S.st}.p{1}))  %���p{1}�洢����NAN����ʼ��ΪNAN��
                node{i+1}{n_S.st}.p{1}  = c_S.st;%�洢��һʱ�̽ڵ��״̬���
                node{i+1}{n_S.st}.visit = 1;
                node{i+1}{n_S.st}.cost  = node{i}{c_S.st}.cost+dist; %����·������
                node{i+1}{n_S.st}.surv  = c_S.st;%��һʱ�̽ڵ��״̬���Ϊ�Ҵ�·�������һ���ڵ�
            else      %���i+1ʱ�̵���ýڵ���Ѿ�����һ�ڵ���Ϣ��iʱ��2���ڵ�״̬ת�Ƶ�i+1ʱ�̽ڵ㣩
                node{i+1}{n_S.st}.p{2} = c_S.st;
                node{i+1}{n_S.st}.visit = 1;
                
                %Two Possible Survivors: Determine surviving branch
                if(node{i+1}{n_S.st}.cost<=node{i}{c_S.st}.cost+dist)%���еĽڵ㵽i+1ʱ�̽ڵ��·������С
                    node{i+1}{n_S.st}.surv  =node{i+1}{n_S.st}.p{1};
                    %�Ҵ�·���洢���еĽڵ�
                else
                    node{i+1}{n_S.st}.surv  = c_S.st;
                    node{i+1}{n_S.st}.cost  = node{i}{c_S.st}.cost+dist; 
                end       
            end
            
            %%%%%%%%%%%%%%%%%%%%(STEP 3-5)%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %State Input = 1; (Binary, q=2�����о��У�����ƽ����)
            if(i<=k)
                %Only process input 1 for first k stages
                c_S.in   = 1; 
                
                 %Determine acceptable output from circuit logic
                [o,n_S] = circuit_logic(c_S,n,mem);%�������һ��ʱ��������ݣ�nλ��
                                                   %�Լ��ڵ�״̬
                
                %Update node's status (e.g., node's status, total cost, is it a
                %possible survivor?)
                if(flag)
                    dist    = compute_Hamm(o,r,block_st,n);
                else
                    dist    = compute_Lp(o,r,block_st,n);  %�����֧����
                end;
                 
                %%%%%%%%%%%%%%%%%%%%(STEP 6)%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                node{i}{c_S.st}.f{2}  = n_S.st;%�洢����Ϊ1ʱ��һ�ڵ��״̬���
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