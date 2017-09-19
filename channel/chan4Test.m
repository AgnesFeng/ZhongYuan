
clear all;
close all;

f=300;                %信源频率      
sample = 1500;
dt = 1/sample;
n = 40000;          %采样点数
T = n*dt;
t = 0:dt:T-dt; 
mt = sqrt(2)*cos(2*pi*f*t);  %信源

ts = dt;
fd = 10;           %最大多普勒频移
tau = [0 0.000001 0.000002 0.000003];
tau1=[0 1 2 3];
pdb = [0 -3 -6 -9];
chan = rayleighchan(ts,fd,tau,pdb); %多径，要求每一径都是瑞丽衰落。ts：输入信号的采样时间（s）；fd：最大多普勒频移；pdb：平均路径增益（dB）
y = filter(chan,mt);       %为复数形式
y1 = abs(y); %衰落后的信道取模
fad = y1./abs(mt); %衰落倍数
fad1 = 10*log10(fad);  %dB
average = mean(fad1); 
rmsdel = sqrt( sum(pdb.*(tau1.^2))/sum(pdb) - (sum(pdb.*tau1)/sum(pdb))^2 );
fprintf('均方根时延rms ： %6.3f μs\n', rmsdel);
figure
subplot(211);
plot(t,mt);
xlabel('时间');ylabel('原始信号幅度')  
subplot(212);
plot(t,y);
xlabel('时间');ylabel('衰落后的信道幅度 ')  
figure
plot(t,fad1);
xlabel('时间');ylabel('小尺度衰落幅度 dB')  
%-----仅小尺度衰落-----%


%----------------加大尺度衰落--------------------
f = 1900e6
L = 3*10^8/f;             %波长  0.1579m
D = 1;                    %天线的最大尺寸 m
df = 2*D^2/L;             %天线远场参考距离  12m
d0 = 100;                 %接受功率的参考点
PL_d0 = 10*log10(L^2/(16*pi^2)*d0^2);  %自由路径损耗
%对于时变信道，路径损耗需考虑随机因素，采用对数正态阴影
X = normrnd(0, 1);        %高斯随机变量
pn = 5;                   %路径损耗指数
dt1 = (1000-100)/n;
d = 100:dt1:1000-dt1;
d1 = log(d);              %log默认是log2（）
PL = PL_d0 - (10*pn*log(d/d0) + X);     %特定距离d下的路径损耗
PSL = PL + fad1;
figure
plot(d1,PL,'b');hold on;
%axis([4.5,6,-80,-30]);
xlabel('发送距离log（d）')
ylabel('纯路径损耗 dB')

plot(d1,PSL,'r');
%axis([4.5,6,-80,-30]);
xlabel('发送距离log（d）')
ylabel('路径损耗加阴影衰落和多径效应 dB')
aa=max(PSL);
bb=min(PSL);
fprintf('四径衰落上限： %6.3f dB\n', bb);
fprintf('四径衰落下限： %6.3f dB\n', aa);
