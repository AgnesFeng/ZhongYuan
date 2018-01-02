clear all;        %仅生成一个周期的m序列画图效果不好
close all;
clc;
%%
%----------------------1、walsh序列(舍弃)---------------------
% wal2 = [ 1 1; 1 -1 ];
% wal4 = [wal2 wal2; wal2 wal2*(-1)];  
% wal8 = [wal4 wal4; wal4 wal4*(-1)];  %8*8
% wal16 = [wal8 wal8; wal8 wal8*(-1)]; %16*16
% p = wal16(5,:);
% KP = length(p); %扩频因子
Tm_sample = 4;%采样率越大，结果图越精准
% pp = SAM_d(p,Tm_sample,length(p));%多周期，有持续时间

%----------------------1、扩频的m序列---------------------
m_spread = m_sequence([0,1,0,0,1]);%能得到31位的m序列
N_ms = length(m_spread);%m长:31
ms1 =  SAM_m(m_spread,Tm_sample,length(m_spread)); %给持续时间
ms2 = 1-2*ms1;  %变双极性
L_ms = length(ms2);
KP = N_ms;
%%%%%%%画图
mms = repmat(ms2,1,3); %3个周期的m序列，方便做相关
rt1=conv(mms,ms2(end:-1:1))/(L_ms);%取中间的一部分
index = L_ms+Tm_sample;
figure(11),subplot(211);plot(rt1(index:end-index));title('周期内m序列的自相关');

xcor = dsp.Crosscorrelator;
y = step(xcor,ms2',ms2'); 
figure(11),subplot(212);plot(y);title('相关器中m序列的自相关');
%%
%%--------------------------2、前缀m序列----------------
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
figure(1);subplot(311);plot(rt1);title('单周期m序列的自相关函数conv，峰值处索引值在m序列长的2倍处');
figure(1);subplot(312);plot(rt1(index:end-index));title('单周期m序列的自相关函数conv,峰值处索引差3个');
xcor = dsp.Crosscorrelator;
y = step(xcor,m2',m2'); 
figure(1);subplot(313);plot(y);title('单周期m序列的自相关函数，相关器，首尾不相接');
%%
%---------------------3、信号序列，扩频,加前导序列---------------------
Tlen = 5000; 
v = 225000;
Tb = 1/v; %码元持续时间  
Tc = Tb/KP;
dt = Tc/Tm_sample;
s_initial = randsrc( 1, Tlen );
T_sample = Tm_sample * KP;  %m序列的采样点数为4
s_initial1 =  SAM_d(s_initial,T_sample,length(s_initial));
pt = repmat(ms2,1,Tlen);%扩频的序列
spreadData = pt.*s_initial1;%扩频
spreadData_m = ones(1,length(m2)+length(s_initial1));%加前导序列
spreadData_m(1:length(m2)) = m2;
spreadData_m(length(m2)+1:end) = spreadData;
len_spread = length(spreadData_m);
%%
%找到相位对齐时的下标
x1 = spreadData_m(1:length(m2));
x2 = m2;
xcor = dsp.Crosscorrelator;
y = step(xcor,x2',x1'); %computes cross-correlation of x1 and x2
figure(2), subplot(211);plot(y); title('Correlated output');
aim = find(max(y) == y);

%%
%--------------------调制整个序列-------------------
txx = mo_bpsk(spreadData_m,v,Tm_sample);
%plot(txx(1:100));
%%
%channel_type = input('请输入信道类型：'); 
channel_type=3;
%-------------------信道类型选择-----------------
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
c_out = filter(chan,txx);
% figure(11),
% subplot(311);plot(s_initial1);axis([0,500,-1.3,1.3]);
% subplot(312),plot(spreadData);axis([0,800,-1.3,1.3]);
% subplot(313),plot(sign(real(c_out)));axis([0,800,-1.3,1.3]);
%--------------------高斯噪声-----------------------
EbNo = 1:1:12;
for snr = 1:length(EbNo)
   y_snr = awgn(c_out,snr,'measured');

%%
    %--------------------省略了解调-------------------
    rx = y_snr;

%%
    %------------调制前导序列和用于扩频的m序列
    m_modu = mo_bpsk(m2,v,Tm_sample);    %调制前导序列
    m_despread = mo_bpsk(pt,v,Tm_sample);  %调制用于扩频的m序列，去软解扩
    %plot(m_modu);
    %------------相关峰捕获
    xcor = dsp.Crosscorrelator;
    buHuo = step(xcor,rx',m_modu'); %computes cross-correlation of x1 and x2
    buHuoQuMo = abs(buHuo);
    figure(2), subplot(212);plot(buHuoQuMo); title('软相关的相关峰');axis([0 2000 0 600]);

    y_max = max(buHuoQuMo);
    id = find(y_max == buHuoQuMo); %找到最大值处，左右开窗,各1码元（延时设在4微秒以内）
    L_catch_window = Tm_sample * KP*2 ;
    feng_x = ones(1,10);
    a = 1;
    b1 = id - L_catch_window/2;
    b2 = id + L_catch_window/2;
    for x =b1 : 1 : b2
        if buHuoQuMo(x) >= y_max/4 && buHuoQuMo(x-1) < buHuoQuMo(x) && buHuoQuMo(x) > buHuoQuMo(x+1) && buHuoQuMo(x-2) < buHuoQuMo(x) && buHuoQuMo(x) > buHuoQuMo(x+2)
            feng_x(a) = x;
            a = a+1;
        end
    end
    fengPaiXu_x = feng_x;
    for j = 1:a-1   %把最大值的横坐标排序，找出最大的三个值的横坐标
        for i = 1:a-1-j
            if buHuoQuMo(fengPaiXu_x(i+1)) > buHuoQuMo(fengPaiXu_x(i))
              temp = fengPaiXu_x(i);
              fengPaiXu_x(i) = fengPaiXu_x(i+1);
              fengPaiXu_x(i+1) = temp;
            end
        end
    end
%%
    %--------------------矫正相位
%     jiaoZheng = step(xcor,rx(L_m+1:end)',m_despread(1:T_sample)');%得要一个周期的m序列
%     figure(3),plot(abs(jiaoZheng)); title('数据部分的相关峰');%axis([0 2000 0 600]);
    
%%
    %----------------判断搜索到几径
    %%%%%%%%捕捉到3个或以上
    if(feng_x(3)~=1)              
      shift_x = ones(1,3);  
      feng = ones(1,3);
      Conj = ones(1,3);
      for i = 1:3
        shift_x(i) = fengPaiXu_x(i)-aim;      %记录位移
        if(shift_x(i)<0)
           shift_x(i)=abs(fengPaiXu_x(i)-aim);
        end
        feng(i) = buHuo(fengPaiXu_x(i));
        Conj(i) = conj(feng(i));
      end      
      path_1 = rx;
      path_2 = rx;
      path_3 = rx;
%       path_1(1,shift_x(1)+1:end) = path_1(1,1:end-shift_x(1));%右移，得到不同延时的三径
%       path_2(1,shift_x(2)+1:end) = path_2(1,1:end-shift_x(2));
%       path_3(1,shift_x(3)+1:end) = path_3(1,1:end-shift_x(3));
      path_1(1,1:end-shift_x(1)) = path_1(1,shift_x(1)+1:end);%左移
      path_2(1,1:end-shift_x(2)) = path_2(1,shift_x(2)+1:end);
      path_3(1,1:end-shift_x(3)) = path_3(1,shift_x(3)+1:end);
      merge = path_1 * conj(1) + path_2 * Conj(2) + path_3 * Conj(3);
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
%     path_1(1,shift_x(1)+1:end) = path_1(1,1:end-shift_x(1)); %右移，得到不同延时的两径
%     path_2(1,shift_x(2)+1:end) = path_2(1,1:end-shift_x(2));
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
  %%
    %%软解扩
    data_out = merge(L_m+1:end); %去掉前缀部分
    despreadData = m_despread.*data_out;%调制后的扩频m序列解扩
    receive = ones(1,Tlen);
    k=1;
    for i = 1:T_sample:length(despreadData)-T_sample+1
        receive(k) = 0;
        for j = 0:T_sample-1
          receive(k) = receive(k) + despreadData(i+j);
        end
        k = k+1;
    end
    receiveJudge = sign(real(receive));%对实部判决
    if(length(find(receiveJudge ~= s_initial))<Tlen/2)
       BitErrorRake = length(find(receiveJudge ~= s_initial)); 
    else
       BitErrorRake = length(find(receiveJudge == s_initial)); 
    end
    BitErrorRake_rate(snr) = BitErrorRake/(Tlen); 
    %%对比直接软解扩
    despreadData_no = m_despread.* rx(L_m+1:end);
    receive2 = ones(1,Tlen);
    k=1;
    for i = 1:T_sample:length(despreadData_no)-T_sample+1
        receive2(k) = 0;
        for j = 0:T_sample-1
          receive2(k) = receive2(k) + despreadData_no(i+j);
        end
        k = k+1;
    end 
    receiveJudge_no = sign(real(receive2));
    number_norake_different = length(find(receiveJudge_no ~= s_initial));
    if(number_norake_different<Tlen/2)
        BitErrorCompare = number_norake_different; 
    else
        BitErrorCompare = length(find(receiveJudge_no == s_initial));
    end
    BitErrorRakeCompare(snr) = BitErrorCompare/(Tlen); 
    
end


figure(111)
subplot(121),semilogy(EbNo,BitErrorRakeCompare,'b');hold on;title('单纯扩频解扩不做同步');
subplot(122),semilogy(EbNo,BitErrorRake_rate,'r');hold on;title('rake接收机');



% %观察一下叠加后的相关峰
% %merge(1,1:end-s(1)) = merge(1,s(1)+1:end);
% y = step(xcor,merge',m2'); %computes cross-correlation of x1 and x2
% y2 = abs(y);
% figure(4);plot(y2); title('Correlated output');
% axis([0 2000 0 50e10000]);



    