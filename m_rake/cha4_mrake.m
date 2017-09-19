function [p_sl] = cha4_mrake(OR,N)
%独立随机数实现的数目， 因为计算的是幅度，采样点数社称已向方便逐点相乘    ;%后面画图横坐标的间隔是1/OR

dt = 1/OR;
P=[0 -3 -6 -9];%每一阶的功率衰减
K=[0 0 0 0];%直线标度中莱斯分布的K因子   莱斯因子
tau=[0 1 2 3];%每个抽头的时延(微秒)

fd = 10;
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
doppler = ray_doppler(fd,dt,N);
for p = 1:L %矩阵中的三项
     paths_r(p,:) = paths_r(p,:).*doppler;
end; 

paths = paths_r + paths_c;  %增益衰减因数(相比SUI，少了增益归一化因子)    归一化因子的意义？？？？？？？？？？
Pest = mean(abs(paths).^2, 2);%平均抽头总功率mean（A，2）：求矩阵A各行的均值，mean（A,1）=mean(A)
fprintf('tap mean power level: %0.2f dB\n', 10*log10(Pest));%打印每行的平均值

paths_OR = paths;
%y1=abs(paths_OR).^2; %功率   倍数     直接取模得到幅度，平方后是功率 

%----------------加大尺度衰落--------------------
f = 1900e6;
L = 3*10^8/f;             %波长  0.1579m
D = 1;                    %天线的最大尺寸 m
df = 2*D^2/L;             %天线远场参考距离  12m
d0 = 100;                 %接受功率的参考点
PL_d0 = 10*log10(L^2/(16*pi^2)*d0^2);  %自由路径损耗
%对于时变信道，路径损耗需考虑随机因素，采用对数正态阴影
X = normrnd(0, 1);        %高斯随机变量
pn = 5;                   %路径损耗指数
d = input('请输入发送距离：'); 
d1 = log(d);              %log默认是log2（）
PL = PL_d0 - (10*pn*log(d/d0) + X);     %特定距离d下的路径损耗
%PSL = PL + pp2(1:901);
pl = 10^(PL./10);         %将log换成倍数
p_sl = pl*ones(4,N);
p_sl(1,:) = abs(paths_OR(1,:))+ p_sl(1,:);
p_sl(2,:) = abs(paths_OR(2,:))+ p_sl(2,:);
p_sl(3,:) = abs(paths_OR(3,:))+ p_sl(3,:);
p_sl(4,:) = abs(paths_OR(4,:))+ p_sl(4,:);
end