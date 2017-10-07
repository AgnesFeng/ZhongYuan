clear all;        %仅生成一个周期的m序列画图效果不好
close all;
clc;
primpoly(7,'all')%得到所有7阶(n为7)的所有本原多项式
m = m_sequence([0,0,0,1,0,0,1]); %得到一个周期的m序列 (各寄存器的初始状态，从左到右为C1~Cn，c0默认为1) 
N = length(m);%m长:127
v = 225000;
Tb = 1/v; %码元持续时间  
Tc = Tb/N;
Tm_sample = 4;%采样率越大，结果图越精准（注意，m序列重复后长度可能不够）
T_sample = Tm_sample * N;  %m序列的采样点数为6
m1 =  SAM_0(m,Tm_sample,length(m),Tc); %给持续时间
m2 = 1-2*m1;  %变双极性
L = length(m2);

NP = 300;
mm = repmat(m,1,NP); %多周期
mm1 = SAM_0(mm,Tm_sample,length(mm),Tc);%多周期，有持续时间
mm2 = 1-2*mm1;       %变双极性
dt = Tc/Tm_sample;
t=0:dt:Tc-dt;
rt1=conv(mm2(1:3*L),m2(end:-1:1))/(N*Tm_sample);%取中间的一部分
figure
index = L+Tm_sample;
plot(rt1);title('单周期m序列的自相关函数(有持续时间)');
% x1 = m2; %函数的实验
% x2 = [1 1 1];
% [a,b] = xcorr(x1,x2);
%---------------------------------原始信号---------------------------------
Tlen = 100; %数据长度  编码后是三倍加9的关系
s_initial = rands( 1, Tlen ); 
figure,plot(s_initial);
s_initial1 =  SAM2(s_initial,T_sample,length(s_initial),Tb);
mseq = mm2(1:length(s_initial1));

spreadData = mseq.*s_initial1;
first = spreadData(1:L);%一个周期的m序列的长度，第一个码元
test = conv(mm2(1:3*L),first(end:-1:1))/L;%移位后最大值会变化，需要设定一个阈值，捕获时滑动

figure
plot(test);
%---------------------------------位移---------------------------------
rr2 = spreadData;   %1us：1/4码元
shift = 3;
rr2(1,round(T_sample/shift+1):end) = rr2(1,1:end-round(T_sample/shift));    %rr2(1,1:T_sample/4)=rr2(1,end-T_sample/4+1:end);
%说明：码元不对齐时有峰值，但不为1，可以根据最大值的位置确定位移，也可以根据滑动检索的距离确定位移，但都需要提前设定阈值
%打算用滑动距离确定，不用最大值的索引
RR2 = rr2(end-L:end);
test1 = conv(mm2(1:3*L),RR2(end:-1:1))/L;
figure
plot(test1);
%观察m序列移位后再相关会不会到最大值
mafter = mm2(1:3*L);
mafter(1,round(T_sample/shift)+1:end) = mafter(1,1:end-round(T_sample/shift));   
mafter(1,1:round(T_sample/shift)) = mafter(1,end-round(T_sample/shift)+1:end);
test2 = conv(mafter,RR2(end:-1:1))/L;
figure
plot(test2);
%---------------------------------捕获---------------------------------

