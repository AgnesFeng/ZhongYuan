clear all;
close all;
fm=25;                %信源频率
T=50;                %信号时长
dt = 0.05;           %时间采样间隔 
sample = 1/dt;       %采样频率
n = T/dt+1;          %采样点数

t = 0:dt:T; 
mt = sqrt(2)*cos(2*pi*fm*t);  %信源
fd = 10;
%function [magn]=rayleigh(fd,t)  
%改进的jakes模型来产生单径的平坦型瑞利衰落信道  
%Yahong R.Zheng and Chengshan Xiao "Improved Models for   
%the Generation of Multiple Uncorrelated Rayleigh Fading Waveforms"   
%IEEE Commu letters, Vol.6, NO.6, JUNE 2002  
%输入变量说明：  
%  fd：信道的最大多普勒频移 单位Hz       
%  t :信号的抽样时间序列，抽样间隔单位s    
%  h为输出的瑞利信道函数，是一个时间函数复序列   
  
    %假设的入射波数目  
    N=30;   
    wm=2*pi*fd;  %最大多普勒频移，w值
    %每象限的入射波数目即振荡器数目  
    N0=N/4;  
    %信道函数的实部  
    Tc=zeros(1,length(mt));  
    %信道函数的虚部  
    Ts=zeros(1,length(mt));  
    %归一化功率系数  
    P_nor=sqrt(1/N0);  
    %区别个条路径的均匀分布随机相位  
    theta=2*pi*rand(1,1)-pi;  
    for ii=1:N0  
          %第i条入射波的入射角   
            alfa(ii)=(2*pi*ii-pi+theta)/N;  
            %对每个子载波而言在(-pi,pi)之间均匀分布的随机相位  
            fi_tc=2*pi*rand(1,1)-pi;  
            fi_ts=2*pi*rand(1,1)-pi;  
            %计算冲激响应函数  
            Tc=Tc+cos(cos(alfa(ii))*wm*mt+fi_tc);  
            Ts=Ts+cos(sin(alfa(ii))*wm*mt+fi_ts);  
    end;  
    %乘归一化功率系数得到传输函数  
   h=P_nor*(Tc+1j*Ts ); %衰落后的信道 
   y1 = abs(h); %衰落后的信道
   fad = y1./mt; %衰落倍数
   magn = 10*log10(fad);  %dB
   figure
   plot(t(1:500),mt(1:500));
   figure
   plot(t(1:500),h(1:500),'r');
   figure
   plot(t,magn);