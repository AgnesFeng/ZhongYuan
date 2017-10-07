function p_sl = c_3(sample,N)
%sample = 1500e6;N = 100000;
dt = 1/sample; %时间采样间隔 
P=[0 -9 -12];%每一阶的功率衰减
K=[6 0 0];%直线标度中莱斯分布的K因子   莱斯因子
tau=[0 1 2];%每个抽头的时延(微秒)
delay = [0,0.000001,0.000002];
fd=25;%最大多普勒频移参数

%%%%%%%%%%%计算每一径莱斯分布中常数部分和随机部分的功率%%%%%%%%%%%%%%%%%%%%%%%
P = 10.^(P/10);%计算线性额定功率
s2 = P./(K+1); % 计算方差
m2 = P.*(K./(K+1)); % 计算常数功率
m = sqrt(m2); % 计算常数部分
%%%%%%%%%%%均方根(Root Mean Square) 延时%%%%%%%%%%%%%%%%
rmsdel = sqrt( sum(P.*(tau.^2))/sum(P) - (sum(P.*tau)/sum(P))^2 );
fprintf('均方根时延rms ： %6.3f μs\n', rmsdel);

%%%%%%%%%在特定功率生成莱斯信道%%%%%%%%%%%%%%%%
L = length(P); % 阶数， L=3
paths_r = sqrt(1/2)*(randn(L,N) + 1i*randn(L,N)).*((sqrt(s2))' * ones(1,N)); %L*N矩阵每阶的数据噪声，与K有关
paths_c = m' * ones(1,N);%常数部分
%位移
t_shift=floor(delay/dt);%归一化各径延时  大于1的
%为各径的输入信号做延迟处理
rr1 = paths_r(1,:); 
rr2 = paths_r(2,:);
rr2(1,t_shift(2)+1:end) = rr2(1,1:end-t_shift(2));
rr3 = paths_r(3,:);
rr3(1,t_shift(3)+1:end) = rr3(1,1:end-t_shift(3));
%叠加
paths_sum = rr1 + rr2 + rr3;
doppler = ray_doppler(fd,dt,N); %25Hz,确实比0.4Hz的影响大
paths_sum1 = paths_sum.*doppler; 
paths_OR = paths_sum1 + paths_c(1,:); %注意，若是分开的单径的，常数部分只加到第一径上面
pp1 = abs(paths_OR); %倍数
pp2 = 10*log10(pp1.^2); %dB

%加大尺度衰落
f = 1900e6;
L = 3*10^8/f;             %波长  0.1579m
D = 1;                    %天线的最大尺寸 m
df = 2*D^2/L;             %天线远场参考距离  12m
d0 = 100;                 %接受功率的参考点
PL_d0 = 10*log10(L^2/(16*pi^2)*d0^2);  %自由路径损耗
%对于时变信道，路径损耗需考虑随机因素，采用对数正态阴影
X = normrnd(0, 1);        %高斯随机变量
pn = 4;                   %路径损耗指数
d = input('请输入发送距离：'); 
PL = PL_d0 - (10*pn*log(d/d0) + X);     %特定距离d下的路径损耗
PSL = PL + pp2; %dB
p_sl = 10.^(PSL./10);
end