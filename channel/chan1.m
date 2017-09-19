% close all;
clear all;
sample = 1500e6;        
N = 3000;

dt = 1/sample; %时间采样间隔 
T = N*dt;
t = 0:dt:T-dt;                            
OR=sample;%观察频率（Hz）

P=0;
K=12;%直线标度中莱斯分布的K因子   莱斯因子
fd=25;%最大多普勒频移参数

%%%%%%%%%%%计算莱斯分布中常数部分和随机部分的功率%%%%%%%%%%%%%%%%%%%%%%%
P = 10.^(P/10);%计算线性额定功率
s2 = P./(K+1); % 计算方差
m2 = P.*(K./(K+1)); % 计算常数功率
m = sqrt(m2); % 计算常数部分

%%%%%%%%%在特定功率生成莱斯信道%%%%%%%%%%%%%%%%
L = length(P); % 阶数， L=3
paths_r = sqrt(1/2)*(randn(L,N) + 1i*randn(L,N)).*((sqrt(s2))' * ones(1,N)); %L*N矩阵每阶的数据噪声，与K有关
paths_c = m' * ones(1,N);%常数部分
doppler = ones(1,N);
doppler = ray_doppler(fd,dt,N);
paths_OR = paths_r.*doppler + paths_c(1,:); 
pp1 = abs(paths_OR);
pp2 = 10*log10(pp1.^2);
x = 1:N;
x1 = x/sample;
figure
plot(x1,pp2);axis([0,0.000002,-20,15])
xlabel('时间（s）观察频率1500MHz'),ylabel('小尺度衰落幅度dB');

%加大尺度衰落
f = 1900e6;
L = 3*10^8/f;             %波长  0.1579m
D = 1;                    %天线的最大尺寸 m
df = 2*D^2/L;             %天线远场参考距离  12m
d0 = 100;                 %接受功率的参考点
PL_d0 = 10*log10(L^2/(16*pi^2)*d0^2);  %自由路径损耗
%对于时变信道，路径损耗需考虑随机因素，采用对数正态阴影
X = normrnd(0, 1);        %高斯随机变量
pn = 2.7;                   %路径损耗指数
dt1 = (1000-100)/N;
d = 100:dt1:1000-dt1;
d1 = log(d);              %log默认是log2（）
PL = PL_d0 - (10*pn*log(d/d0) + X);     %特定距离d下的路径损耗
PSL = PL + pp2;
figure
plot(d1,PL,'b');hold on;
plot(d1,PSL,'r');hold on;
%axis([4.5,6,-80,-30]);
xlabel('发送距离log（d）');
ylabel('纯路径损耗 dB');
