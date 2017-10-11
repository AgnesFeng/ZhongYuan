clear all;        %仅生成一个周期的m序列画图效果不好
close all;
clc;
%%
%----------------------2、walsh序列---------------------
wal2 = [ 1 1; 1 -1 ];
wal4 = [wal2 wal2; wal2 wal2*(-1)];  
wal8 = [wal4 wal4; wal4 wal4*(-1)];  %8*8
wal16 = [wal8 wal8; wal8 wal8*(-1)]; %16*16
p = wal16(5,:);
KP = length(p); %扩频因子
Tm_sample = 4;%采样率越大，结果图越精准
pp = SAM_d(p,Tm_sample,length(p));%多周期，有持续时间
%%
%%--------------------------1、m序列----------------
m = m_sequence([0,0,0,1,1,1,0,1]); %得到一个周期的m序列 (各寄存器的初始状态，从左到右为C1~Cn，c0默认为1) 
N = length(m);%m长:225
m1 =  SAM_m(m,Tm_sample,length(m)); %给持续时间
m2 = 1-2*m1;  %变双极性
L = length(m2);
peri = 3;
mm = repmat(m,1,peri); %3个周期的m序列，方便做相关
mm1 = SAM_m(mm,Tm_sample,length(mm));%多周期，有持续时间
mm2 = 1-2*mm1;       %变双极性

rt1=conv(mm2(1:3*L),m2(end:-1:1))/(N*Tm_sample);%取中间的一部分
figure
index = L+Tm_sample;
figure(1);subplot(211);plot(rt1(index:end-index));title('单周期m序列的自相关函数conv,最大值的位置差三个');
xcor = dsp.Crosscorrelator;
y = step(xcor,m2',m2'); 
figure(1);subplot(212);plot(y);title('单周期m序列的自相关函数，相关器，首尾不相接');
%%
%---------------------3、信号序列，扩频,加前导序列---------------------
Tlen = 6000; 
v = 225000;
Tb = 1/v; %码元持续时间  
Tc = Tb/KP;
dt = Tc/Tm_sample;
s_initial = randsrc( 1, Tlen );
T_sample = Tm_sample * length(p);  %m序列的采样点数为4
s_initial1 =  SAM_d(s_initial,T_sample,length(s_initial));
pt = repmat(pp,1,Tlen);
spreadDate = pt.*s_initial1;%扩频
spreadDate_m = ones(1,length(m2)+length(s_initial1));%加前导序列
spreadDate_m(1:length(m2)) = m2;
spreadDate_m(length(m2)+1:end) = s_initial1;
len_spread = length(spreadDate_m);
%%
%找到相位对齐时的下标
x1 = spreadDate_m(1:length(m2));
x2 = m2;
xcor = dsp.Crosscorrelator;
%  delay = dsp.Delay(30);
%  x2 = step(delay,m2');
y = step(xcor,x2',x1'); %computes cross-correlation of x1 and x2
figure(2), subplot(211);plot(y); title('Correlated output')

%%
%--------------------调制-------------------
%spreadDate_m = rand(1,128);len_spread = length(spreadDate_m);
%fs = 1/dt;
tx = ones(1,len_spread/2);%串并转换后的数组
tx_re = tx; 
tx_im = tx;
k2 = 1;
for k1=1:Tm_sample*2:len_spread-2*Tm_sample+1 
    tx_re(k2:k2+Tm_sample-1) = spreadDate_m(k1:k1+Tm_sample-1);
    tx_im(k2:k2+Tm_sample-1) = spreadDate_m(k1+Tm_sample:k1+2*Tm_sample-1);
    k2 = k2+Tm_sample;
end

txx = tx_re + 1i*tx_im;
%--------------------过信道-----------------
 ts =dt;
fd = 25;           %最大多普勒频移 
k = 6;
tau = [0 0.000001 0.000002];
pdb = [0,-2,-2];
chan = ricianchan(ts,fd,k,tau,pdb); %多径，要求每一径都是瑞丽衰落。ts：输入信号的采样时间（s）；fd：最大多普勒频移；pdb：平均路径增益（dB）

% fd = 10;           %最大多普勒频移 
% k = 6;
% tau = [0 0.000001 0.000002 0.000003];
% pdb = [0,-3,-6,-9];
% chan = rayleighchan(ts,fd,tau,pdb);

y = filter(chan,txx); 
EbNo = 0:1:15;
for snr = 1:length(EbNo)
   y_snr = awgn(y,snr,'measured');

    %--------------------解调-------------------
    re_2=real(y_snr);
    im_2=imag(y_snr);
    rx = ones(1,len_spread);
    k2 = 1; 
    for k1=1:Tm_sample*2:len_spread-2*Tm_sample+1 
        rx(k1:k1+Tm_sample-1) = re_2(k2:k2+Tm_sample-1);
        rx(k1+Tm_sample:k1+2*Tm_sample-1) = im_2(k2:k2+Tm_sample-1);
        k2 = k2+Tm_sample;
    end
   %判决
   receive = sign(rx);%判决后其他径的峰值大一些
   %receive = rx; 
Bit_error = length(find(receive ~= spreadDate_m)); 
error_rate(snr) = Bit_error/len_spread; 
end

%%
%------------相关峰捕获
x2 = rx(1:length(m2));
xcor = dsp.Crosscorrelator;
y = step(xcor,x2',m2'); %computes cross-correlation of x1 and x2
figure(2), subplot(212);plot(abs(y)); title('Correlated output')


