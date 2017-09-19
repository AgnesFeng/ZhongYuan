function [selq] = SAM_0(mt,sample,NUM,Ts)
%方波采样函数。 原序列，单个脉冲采样点数，码元（脉冲）数,码元周期

N_sample = sample; % 单个码元抽样点数
dt = Ts / N_sample; % 抽样时间间隔
N = NUM; % 码元数
t = 0 : dt : (N * N_sample - 1) * dt; % 序列传输时间
gt1 = ones(1, N_sample); % 不归零波形
% 生成随机序列
RAN = round(mt); % 随机0 1序列
selq = [];
for i = 1 : N % 生成序列
   if RAN(i)==1
       selq = [selq gt1]; 
   else
       selq = [selq 0*gt1];
   end
end

end