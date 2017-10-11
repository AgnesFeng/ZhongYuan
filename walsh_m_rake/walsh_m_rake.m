clear all;        %仅生成一个周期的m序列画图效果不好
close all;
clc;
%%
%----------------------1、walsh序列---------------------
wal2 = [ 1 1; 1 -1 ];
wal4 = [wal2 wal2; wal2 wal2*(-1)];  
wal8 = [wal4 wal4; wal4 wal4*(-1)];  %8*8
wal16 = [wal8 wal8; wal8 wal8*(-1)]; %16*16
p = wal16(5,:);
KP = length(p); %扩频因子
Tm_sample = 4;%采样率越大，结果图越精准
pp = SAM_d(p,Tm_sample,length(p));%多周期，有持续时间
%%
%%--------------------------2、m序列----------------
m = m_sequence([0,0,0,1,1,1,0,1]); %得到一个周期的m序列 (各寄存器的初始状态，从左到右为C1~Cn，c0默认为1) 
N_m = length(m);%m长:225
m1 =  SAM_m(m,Tm_sample,length(m)); %给持续时间
m2 = 1-2*m1;  %变双极性
L_m = length(m2);
peri = 3;
mm = repmat(m,1,peri); %3个周期的m序列，方便做相关
mm1 = SAM_m(mm,Tm_sample,length(mm));%多周期，有持续时间
mm2 = 1-2*mm1;       %变双极性

rt1=conv(mm2(1:3*L_m),m2(end:-1:1))/(N_m*Tm_sample);%取中间的一部分
index = L_m+Tm_sample;
% figure(1);subplot(211);plot(rt1(index:end-index));title('单周期m序列的自相关函数conv,最大值的位置差三个');
% xcor = dsp.Crosscorrelator;
% y = step(xcor,m2',m2'); 
% figure(1);subplot(212);plot(y);title('单周期m序列的自相关函数，相关器，首尾不相接');
%%
%---------------------3、信号序列，扩频,加前导序列---------------------
Tlen = 6000; 
v = 225000;
Tb = 1/v; %码元持续时间  
Tc = Tb/KP;
dt = Tc/Tm_sample;
s_initial = randsrc( 1, Tlen );
T_sample = Tm_sample * KP;  %m序列的采样点数为4
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
y = step(xcor,x2',x1'); %computes cross-correlation of x1 and x2
figure(2), subplot(211);plot(y); title('Correlated output');
aim = find(max(y) == y);

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
%%
%--------------------过信道-----------------
 ts =dt;
% fd = 25;           %最大多普勒频移 
% k = 6;
% tau = [0 0.000001 0.000002];
% pdb = [0,-9,-12];
% chan = ricianchan(ts,fd,k,tau,pdb); %多径，要求每一径都是瑞丽衰落。ts：输入信号的采样时间（s）；fd：最大多普勒频移；pdb：平均路径增益（dB）

fd = 10;           %最大多普勒频移 
k = 6;
tau = [0 0.000001 0.000002 0.000003];
pdb = [0,-3,-6,-9];
chan = rayleighchan(ts,fd,tau,pdb);

c_out = filter(chan,txx); 
EbNo = 29:1:30;
for snr = 1:length(EbNo)
   y_snr = awgn(c_out,snr,'measured');

% psl = c_3(1/dt,500); %（采样频率，序列长度，序列）变三径
% rxx = filter(psl,1,txx);
% EbN0db  = 0:1:8;
% for snr = 1:length( EbN0db ) 
%     y_snr = awgn(rxx,snr,'measured'); 
%%   
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
   receive1 = sign(rx);%判决后其他径的峰值大一些
   Bit_error1 = length(find(receive1 ~= spreadDate_m)); 
   error_rate1(snr) = Bit_error1/len_spread; 

%%
    %------------相关峰捕获
    x2 = rx(1:length(m2));
    xcor = dsp.Crosscorrelator;
    y = step(xcor,x2',m2'); %computes cross-correlation of x1 and x2
    y1 = abs(y);
    figure(2), subplot(212);plot(y1); title('Correlated output')

    y_max = max(y1);
    id = find(y_max == y1); %找到最大值处，左右开窗,各1码元（延时设在4微秒以内）
    L_catch = Tm_sample * KP *2;
    tempx = ones(1,10);
    a = 1;
    b1 = id - L_catch/2;
    b2 = id + L_catch/2;
    for x =b1 : 1 : b2
        if y1(x) >= y_max/5 && y1(x-1) < y1(x) && y1(x) > y1(x+1)
            tempx(a) = x;
            a = a+1;
        end
    end
    tempxx = tempx;
    for j = 1:a-1   %把最大值的横坐标排序，找出最大的三个值的横坐标
        for i = 1:a-1-j
            if y1(tempxx(i+1)) > y1(tempxx(i))
              temp = tempxx(i);
              tempxx(i) = tempxx(i+1);
              tempxx(i+1) = temp;
            end
        end
    end
%%
    %----------------相位对齐
    delay1 = dsp.Delay(tempxx(1)-aim); %delay好像不可以是数组
    delay2 = dsp.Delay(tempxx(2)-aim);
    delay3 = dsp.Delay(tempxx(3)-aim); %有时候没有捕捉到第三个峰，要重跑一次
    path(1,:) = step(delay1,rx');
    path(2,:) = step(delay2,rx');
    path(3,:) = step(delay3,rx');
    A = 10.^(pdb./10);
    merge = path(1,:) *A(1) + path(2,:) *A(2) + path(3,:) *A(3); %将增益化成倍数
    receive2 = sign(merge);
    Bit_error2 = length(find(receive2 ~= spreadDate_m)); 
    error_rate2(snr) = Bit_error2/len_spread; 
end
%相位对齐要再思考
%叠加后观察一下相关峰
bef = min(tempxx(1:3))-aim;
aa = merge(bef+1:bef+L_m);
y = step(xcor,aa',m2'); %computes cross-correlation of x1 and x2
y1 = abs(y);
figure(4), subplot(212);plot(y1); title('Correlated output')

    