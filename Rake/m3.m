clear all;
close all;
g=19; % G=10011  n=4  m序列长15
state=8; % state=1000  移位寄存器初始状态
L=1000;
% m序列产生
N=15;
mq=mgen(g,state,L); %长为1000，周期为15的m序列
mq1 = [1 1 1 0 0 1 0];
% 求序列自相关
ms=conv(1-2*mq,1-2*mq(15:-1:1))/N;%前后顺序改变无所谓，其中一个倒着与另一个做卷积，必须是“一个周期的倒着，与整个周期卷积”
[ms1,b] = xcorr(1-2*mq,'biased');
figure(1)
subplot(222)
stem(ms(15:end));% 画针状图 (头尾多了15个)
%axis([0 63 -0.3 1.2]); title('m序列的自相关序列')

% m序列构成的信号（矩形脉冲）
N_sample=8;
Tc=1;
dt=Tc/N_sample;
t=0:dt:Tc*L-dt;  %%%%%%5
gt=ones(1,N_sample);
mt=sigexpand(1-2*mq,N_sample);%采样，给持续时间，mq中的1变-1，
mt=conv(mt,gt); %矩形成形 是双极性的，全长
figure(1)
subplot(221);
plot(t,mt(1:length(t)));
axis([0 63 -1.2 1.2]); title('m序列矩形成形信号');
% 
st=sigexpand(1-2*mq(1:15),N_sample); %双极性m序列，一个周期
s=conv(st,gt);     %给持续时间
st=s(1:length(st));%卷积后长度变了，取第一个周期
rt1=conv(mt,st(end:-1:1))/(N*N_sample);%最后一行到第一行换反过来,做相关
subplot(223)
plot(t,rt1(length(st):length(st)+length(t)-1));
%plot(t,rt1(1:length(t)));
%plot(rt1);
axis([0 63 -0.2 1.2]); title('m序列矩形成形信号的自相关');xlabel('t');















% sinc脉冲
% Tc=1;
% dt=Tc/N_sample;
% t=-20:dt:20;
% gt=sinc(t/Tc);
% mt=sigexpand(1-2*mq,N_sample);
% mt=conv(mt,gt);
% st2=sigexpand(1-2*mq(1:15),N_sample);
% s2=conv(st2,gt);
% st2=s2;
% rt2=conv(mt,st2(end:-1:1))/(N*N_sample);
% subplot(224)
% t1=-55+dt:dt:Tc*L-dt;
% plot(t1,rt2(1:length(t1)));
% axis([0 63 -0.5 1.2]); title('m序列sinc成形信号的自相关');xlabel('t');