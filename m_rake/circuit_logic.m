function [o,next_State] = circuit_logic(cur_State,n,m)
%Function:  This defines the circuit logic for a specific convolution
%           encoder.  For now, I hand code the n outputs, but this can be
%           easily done automatically via the generator impulse responses.

%cur_State - The current state of the filter
%n         - number of output words [y(1),y(2),...,y(n)]
%m         - number of memory elements

%Define output (n bit) in terms of states？？？？？？
y{1} = mod(cur_State.in + cur_State.m{2} +cur_State.m{3},2);
%g*[cur_State.in;cur_State[m](end:-1:1)];
%y{1} = mod(cur_State.in +cur_State.m{3}+ cur_State.m{4}+ cur_State.m{5},2);
%y{1} = mod(cur_State.in +cur_State.m{3}+ cur_State.m{5}+ cur_State.m{7},2);
y{2} = mod(cur_State.in + cur_State.m{1} + cur_State.m{3},2);
%y{2} = mod(cur_State.in + cur_State.m{2} + cur_State.m{4}+ cur_State.m{5},2);
%y{2} = mod(cur_State.in + cur_State.m{1}+cur_State.m{3}+ cur_State.m{4}+ cur_State.m{7},2);
y{3} = mod(cur_State.in + cur_State.m{1} + cur_State.m{2} + cur_State.m{3},2);
%y{3} = mod(cur_State.in + cur_State.m{1} + cur_State.m{2} + cur_State.m{3}...
% + cur_State.m{5},2);
%y{3} = mod(cur_State.in +cur_State.m{1}+ cur_State.m{2}+ cur_State.m{3}...
  %  + cur_State.m{5}+ cur_State.m{6}+ cur_State.m{7},2);
   
%%从前一个时刻的节点状态得到输出数据（n位）

%Initialize Output Word
o = zeros(1,n*length(y{1}));
for i = 0:n-1; o(1+i:n:end) = y{i+1}; end %并串转换
o(o==0) = -1;%极性变换，r也经过了极性变换

next_State.st = 0;
%Convert binary vec to state value, and Update State.
for i = 0:m-1
    if(i+1==1); next_State.m{i+1} = cur_State.in;
    else       next_State.m{i+1} = cur_State.m{i};%得到下一节点的状态序号
    end
    next_State.st = next_State.st+ (2^i)*next_State.m{i+1};%.m的数据二进制变成十进制（节点状态的顺序）
end
next_State.st = next_State.st+1;%节点状态序号从1开始编号
end