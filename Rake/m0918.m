clear all;        %仅生成一个周期的m序列画图效果不好
close all;
clc;
primpoly(7,'all')%得到所有7阶(n为7)的所有本原多项式
m = m_sequence([0,0,0,1,1,1,0,1]); %得到一个周期的m序列 (各寄存器的初始状态，从左到右为C1~Cn，c0默认为1) 
N = length(m);
Tc = 1;
Tm_sample = 6;%采样率越大，结果图越精准（注意，m序列重复后长度可能不够）
T_sample = 6 * N;  %m序列的采样点数为6
m1 =  SAM_0(m,Tm_sample,length(m),Tc); %给持续时间
m2 = 1-2*m1;  %变双极性
L = length(m2);

NP = 3;
mm = repmat(m,1,NP); %多周期
mm1 = SAM_0(mm,Tm_sample,length(mm),Tc);%多周期，有持续时间
mm2 = 1-2*mm1;       %变双极性
dt = Tc/Tm_sample;
t=0:dt:Tc-dt;
rt1=conv(mm2,m2(end:-1:1))/(N*Tm_sample);%取中间的一部分
figure(2)
index = L+Tm_sample;
plot(rt1(index:end-index));title('单周期m序列的自相关函数(有持续时间)');
% x1 = m2; %函数的实验
% x2 = [1 1 1];
% [a,b] = xcorr(x1,x2);
%--------------------------3、设置阈值------(应该放在扩频后的序列)-----------------

