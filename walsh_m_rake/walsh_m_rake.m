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
figure(1);subplot(211);plot(rt1(index:end-index));title('单周期m序列的自相关函数conv,最大值的位置差三个');
xcor = dsp.Crosscorrelator;
y = step(xcor,m2',m2'); 
figure(1);subplot(212);plot(y);title('单周期m序列的自相关函数，相关器，首尾不相接');
%%
%---------------------3、信号序列，扩频,加前导序列---------------------
Tlen = 3000; 
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
%--------------------qpsk调制，仅串并转换，无键控变换-------------------
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
%-------------------信道类型选择-----------------
%channel_type = input('请输入信道类型：'); 
channel_type=2;

ts =dt;
if(channel_type==1)
    fd = 25;
    k = 10^(12/10);
    chan = ricianchan(ts,fd,k);  
elseif(channel_type==2)
    fd = 25;           %最大多普勒频移 
    k = 10^(6/10);
    tau = [0 0.000001 0.000002];
    pdb = [0,-9,-12];
    chan = ricianchan(ts,fd,k,tau,pdb); %多径，要求每一径都是瑞丽衰落。ts：输入信号的采样时间（s）；fd：最大多普勒频移；pdb：平均路径增益（dB）
elseif(channel_type==3)
    fd = 10;     %最大多普勒频移 
    tau = [0 0.000001 0.000002 0.000003];
    pdb = [0,-3,-6,-9];
    chan = rayleighchan(ts,fd,tau,pdb);
else                                    % (channel_type~=3 &&channel_type~=1 &&channel_type~=2)
    error('Error! channel type should be one of 1 2 3!');
end 
%%
%--------------------过信道-----------------
c_out = filter(chan,txx);
EbNo = 1:1:3;
for snr = 1:length(EbNo)
   y_snr = awgn(c_out,snr,'measured');
   %y_snr = c_out;
% psl = c_3(1/dt,500); %（采样频率，序列长度，序列）变三径
% rxx = filter(psl,1,txx);
% EbN0db  = 0:1:8;
% for snr = 1:length( EbN0db ) 
%     y_snr = awgn(rxx,snr,'measured'); 
%%   
    %--------------------解调，并串转换-------------------
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
   BitErrorBeforeRake = length(find(receive1 ~= spreadDate_m)); 
   BitErrorBeforeRake_rate(snr) = BitErrorBeforeRake/(len_spread); 

%%
    %------------相关峰捕获
    x2 = rx(1:length(m2));
    xcor = dsp.Crosscorrelator;
    buHuo = step(xcor,x2',m2'); %computes cross-correlation of x1 and x2
    y1 = abs(buHuo);
    figure(2), subplot(212);plot(y1); title('解调后的相关峰');

    y_max = max(y1);
    id = find(y_max == y1); %找到最大值处，左右开窗,各1码元（延时设在4微秒以内）
    L_catch = Tm_sample * KP*2 ;
    feng_x = ones(1,10);
    a = 1;
    b1 = id - L_catch/2;
    b2 = id + L_catch/2;
    for x =b1 : 1 : b2
        if y1(x) >= y_max/4 && y1(x-1) < y1(x) && y1(x) > y1(x+1) && y1(x-2) < y1(x) && y1(x) > y1(x+2)
            feng_x(a) = x;
            a = a+1;
        end
    end
    fengPaiXu_x = feng_x;
    for j = 1:a-1   %把最大值的横坐标排序，找出最大的三个值的横坐标
        for i = 1:a-1-j
            if y1(fengPaiXu_x(i+1)) > y1(fengPaiXu_x(i))
              temp = fengPaiXu_x(i);
              fengPaiXu_x(i) = fengPaiXu_x(i+1);
              fengPaiXu_x(i+1) = temp;
            end
        end
    end
%%
    %----------------相位对齐
    %%%%%%%%捕捉到3个或以上
    if(feng_x(3)~=1)   
        shift_x = ones(1,3);  
        feng = ones(1,3);
        Conj = ones(1,3);  
        for i = 1:3
          shift_x(i) = fengPaiXu_x(i)-aim; 
          feng(i) = buHuo(fengPaiXu_x(i));
          Conj(i) = conj(feng(i));
        end
        path_1 = rx;
        path_2 = rx;
        path_3 = rx;
        path_1(1,1:end-shift_x(1)) = path_1(1,shift_x(1)+1:end);%左移
        path_2(1,1:end-shift_x(2)) = path_2(1,shift_x(2)+1:end);
        path_3(1,1:end-shift_x(3)) = path_3(1,shift_x(3)+1:end);
        p = y1(fengPaiXu_x(1)) + y1(fengPaiXu_x(2)) + y1(fengPaiXu_x(3));
        %p = y1(tempxx(1)) + y1(tempxx(2)) + y1(tempxx(3))+y1(tempxx(4));   %计算每一径的加权系数
        u1 = y1(fengPaiXu_x(1))/p;             
        u2 = y1(fengPaiXu_x(2))/p;
        u3 = y1(fengPaiXu_x(3))/p; 
        merge = path_1 * u1 + path_2 * u2 + path_3 * u3;
        %merge = path_1 * u1 + path_2 * u2 + path_3 * u3 + path_4 * u4; %权重要对应
    end
    %%%%%%%%%%%%%%%%%%捕捉到2个
    if(feng_x(2)~=1 && feng_x(3)==1) 
      shift_x = ones(1,2);  
      feng = ones(1,2);
      Conj = ones(1,2);
      for i = 1:2
        shift_x(i) = fengPaiXu_x(i)-aim; 
        feng(i) = buHuo(fengPaiXu_x(i));
        Conj(i) = conj(feng(i));
      end    
      path_1 = rx;
      path_2 = rx;
      path_1(1,1:end-shift_x(1)) = path_1(1,shift_x(1)+1:end);%左移
      path_2(1,1:end-shift_x(2)) = path_2(1,shift_x(2)+1:end);
      merge = path_1 * Conj(1) + path_2 * Conj(2);
    end
    %%%%%%%%%%%%%%%%%%%捕捉到1个
    if(feng_x(1)~=1 && feng_x(3)==1 && feng_x(2)==1) 
      shift_x = fengPaiXu_x(1)-aim; 
      feng = buHuo(fengPaiXu_x(1));
      Conj = conj(feng(1));  
      path_1 = rx;
      path_1(1,1:end-shift_x(1)) = path_1(1,shift_x(1)+1:end);%左移
      merge = path_1 * Conj(1);
    end
    %%%%解扩前计算误码率
    receive2 = sign(merge);
    BitErrorRake = length(find(receive2 ~= spreadDate_m)); 
    BitErrorRake_rate(snr) = BitErrorRake/(len_spread); 
end
%叠加后观察一下相关峰
y = step(xcor,merge(1:L_m)',m2'); %computes cross-correlation of x1 and x2
y2 = abs(y);
figure(4), subplot(212);plot(y2); title('Correlated output')

    