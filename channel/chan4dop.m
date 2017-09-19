clear all;
close all;
fm=25;                %��ԴƵ��
T=50;                %�ź�ʱ��
dt = 0.05;           %ʱ�������� 
sample = 1/dt;       %����Ƶ��
n = T/dt+1;          %��������

t = 0:dt:T; 
mt = sqrt(2)*cos(2*pi*fm*t);  %��Դ
fd = 10;
%function [magn]=rayleigh(fd,t)  
%�Ľ���jakesģ��������������ƽ̹������˥���ŵ�  
%Yahong R.Zheng and Chengshan Xiao "Improved Models for   
%the Generation of Multiple Uncorrelated Rayleigh Fading Waveforms"   
%IEEE Commu letters, Vol.6, NO.6, JUNE 2002  
%�������˵����  
%  fd���ŵ�����������Ƶ�� ��λHz       
%  t :�źŵĳ���ʱ�����У����������λs    
%  hΪ����������ŵ���������һ��ʱ�亯��������   
  
    %��������䲨��Ŀ  
    N=30;   
    wm=2*pi*fd;  %��������Ƶ�ƣ�wֵ
    %ÿ���޵����䲨��Ŀ��������Ŀ  
    N0=N/4;  
    %�ŵ�������ʵ��  
    Tc=zeros(1,length(mt));  
    %�ŵ��������鲿  
    Ts=zeros(1,length(mt));  
    %��һ������ϵ��  
    P_nor=sqrt(1/N0);  
    %�������·���ľ��ȷֲ������λ  
    theta=2*pi*rand(1,1)-pi;  
    for ii=1:N0  
          %��i�����䲨�������   
            alfa(ii)=(2*pi*ii-pi+theta)/N;  
            %��ÿ�����ز�������(-pi,pi)֮����ȷֲ��������λ  
            fi_tc=2*pi*rand(1,1)-pi;  
            fi_ts=2*pi*rand(1,1)-pi;  
            %����弤��Ӧ����  
            Tc=Tc+cos(cos(alfa(ii))*wm*mt+fi_tc);  
            Ts=Ts+cos(sin(alfa(ii))*wm*mt+fi_ts);  
    end;  
    %�˹�һ������ϵ���õ����亯��  
   h=P_nor*(Tc+1j*Ts ); %˥�����ŵ� 
   y1 = abs(h); %˥�����ŵ�
   fad = y1./mt; %˥�䱶��
   magn = 10*log10(fad);  %dB
   figure
   plot(t(1:500),mt(1:500));
   figure
   plot(t(1:500),h(1:500),'r');
   figure
   plot(t,magn);