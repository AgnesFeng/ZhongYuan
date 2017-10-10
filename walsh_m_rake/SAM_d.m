function [selq] = SAM_d(mt,sample,NUM)
%方波采样函数。 原序列，单个脉冲采样点数，码元（脉冲）数,码元周期
N = NUM; % 码元数
gt1 = ones(1, sample); % 不归零波形
% 生成随机序列
RAN = round(mt); % 随机0 1序列
selq = [];
for i = 1 : N % 生成序列
   if RAN(i)==1
       selq = [selq gt1]; 
   else
       selq = [selq -1*gt1];
   end
end

end