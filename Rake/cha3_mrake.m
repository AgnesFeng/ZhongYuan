function [P_sl] = cha3_mrake(sample,n)
N=n;%独立随机数实现的数目， 因为计算的是幅度，采样点数社称已向方便逐点相乘

OR=sample;%后面画图横坐标的间隔是1/OR
M=256;%多普勒滤波器的阶数

P=[0 -9 -12];%每一阶的功率衰减
K=[6 0 0];%直线标度中莱斯分布的K因子   莱斯因子
tau=[0 1 2];%每个抽头的时延(微秒)
Dop=[25 25 25];%最大多普勒频移参数

%%%%%%%%%%%计算每一径莱斯分布中常数部分和随机部分的功率%%%%%%%%%%%%%%%%%%%%%%%
P = 10.^(P/10);%计算线性额定功率
s2 = P./(K+1); % 计算方差
m2 = P.*(K./(K+1)); % 计算常数功率
m = sqrt(m2); % 计算常数部分
%%%%%%%%%%%均方根(Root Mean Square) 延时%%%%%%%%%%%%%%%%
rmsdel = sqrt( sum(P.*(tau.^2))/sum(P) - (sum(P.*tau)/sum(P))^2 );
fprintf('rms delay spread %6.3f μs\n', rmsdel);
%%%%%%%%%在特定功率生成莱斯信道%%%%%%%%%%%%%%%%
L = length(P); % 阶数， L=3
paths_r = sqrt(1/2)*(randn(L,N) + 1i*randn(L,N)).*((sqrt(s2))' * ones(1,N)); %L*N矩阵每阶的数据噪声，与K有关
paths_c = m' * ones(1,N);%常数部分

for p = 1:L %矩阵中的三项
    D = Dop(p) / max(Dop) / 2; % 归一化最大多普勒频移 
    f0 = (0:M*D)/(M*D); % 频率因子 
    PSD = 0.785*f0.^4 - 1.72*f0.^2 + 1.0; % PSD近似法
    filt = [ PSD(1:end-1) zeros(1,floor(M-2*M*D)) PSD(end:-1:2) ]; % S(f) 
    
    filt = sqrt(filt); %从S(f)到|H(f)| 
    filt = ifftshift(ifft(filt)); % 获得脉冲响应 
    filt = real(filt); % 寻找实数滤波器
    filt = filt / sqrt(sum(filt.^2)); %归一化滤波器
    path = fftfilt(filt, [ paths_r(p,:) zeros(1,M) ]); 
    paths_r(p,:) = path(1+M/2:end-M/2); %随机部分
end; 
paths = paths_r + paths_c;  %增益衰减因数(相比SUI，少了增益归一化因子)    归一化因子的意义？？？？？？？？？？
Pest = mean(abs(paths).^2, 2);%平均抽头总功率mean（A，2）：求矩阵A各行的均值，mean（A,1）=mean(A)
fprintf('tap mean power level: %0.2f dB\n', 10*log10(Pest));%打印每行的平均值

paths_OR = paths;
paths_OR1=10*log10(paths_OR); %方便观察数据用
%NN=length(paths_OR(1,:)); %矩阵的列数
y1=abs(paths_OR).^2; %功率   倍数     直接取模得到幅度，平方后是功率
y2=10*log10(y1); %转换为dB
%t=60; %时间长度
x=1:n; %(换成n就不行)
x1=x/OR;

%%%%%%%%%%%%%%三径小尺度衰落倍数%%%%%%%%%%%%%
pp = paths_OR(1,:) + paths_OR(2,:) + paths_OR(3,:); %接收端同一时刻到达的幅度，因此不加延时直接叠加 
pp1 = abs(pp);  %倍数
pp2 = 10*log10(pp1.^2); %功率dB

%----------------加大尺度衰落--------------------
f = 1900e6;
L = 3*10^8/f;             %波长  0.1579m
D = 1;                    %天线的最大尺寸 m
df = 2*D^2/L;             %天线远场参考距离  12m
d0 = 100;                 %接受功率的参考点
PL_d0 = 10*log10(L^2/(16*pi^2)*d0^2);  %自由路径损耗
%对于时变信道，路径损耗需考虑随机因素，采用对数正态阴影
X = normrnd(0, 1);        %高斯随机变量
pn = 2.8;                   %路径损耗指数
d = input('请输入发送距离：'); 
d1 = log(d);              %log默认是log2（）
PL = PL_d0 - (10*pn*log(d/d0) + X);     %特定距离d下的路径损耗
%PSL = PL + pp2(1:901);
pl = 10^(PL./10);         %将log换成倍数
P_sl = pl*ones(3,n);
P_sl(1,:) = abs(paths_OR(1,:))+ P_sl(1,:);
P_sl(2,:) = abs(paths_OR(2,:))+ P_sl(2,:);
P_sl(3,:) = abs(paths_OR(3,:))+ P_sl(3,:);
end