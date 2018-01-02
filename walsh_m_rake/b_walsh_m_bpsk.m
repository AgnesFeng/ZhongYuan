clear all;        %仅生成一个周期的m序列画图效果不好
close all;
clc;
%%
%----------------------1、walsh序列---------------------
wal2 = [ 1 1; 1 -1 ];
wal4 = [wal2 wal2; wal2 wal2*(-1)];  
wal8 = [wal4 wal4; wal4 wal4*(-1)];  %8*8
wal16 = [wal8 wal8; wal8 wal8*(-1)]; %16*16
wal32 = [wal16 wal16; wal16 wal16*(-1)]; %32*32
p_0 = wal32(5,:);
p_1 = wal32(11,:);
KP = 32; %扩频因子
Tm_sample = 4;%采样率越大，结果图越精准
pp0 = SAM_d(p_0,Tm_sample,length(p_0));%多周期，有持续时间
pp1 = SAM_d(p_1,Tm_sample,length(p_1));%多周期，有持续时间
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
% %---------------------3、信号序列，扩频,加前导序列---------------------
Tlen = 1000; 
v = 225000;
Tb = 1/v; %码元持续时间  
Tc = Tb/KP;
dt = Tc/Tm_sample;
s_initial = randsrc( 1, Tlen );
T_sample = Tm_sample * KP;  %m序列的采样点数为4
s_initial1 =  SAM_d(s_initial,T_sample,length(s_initial));
L_s_initial1 = length(s_initial1);
spreadData = ones(1,L_s_initial1);
for i = 1: T_sample : L_s_initial1-T_sample+1 %扩频
     if(s_initial1(i) == -1)
         spreadData(i : i+T_sample-1) = pp0.*s_initial1(i : i+T_sample-1);%对为-1的码元扩频
     else
         spreadData(i : i+T_sample-1) = pp1.*s_initial1(i : i+T_sample-1);%对为 1的码元扩频
     end
end
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
%txx = mo_bpsk(spreadData_m,v,Tm_sample);
%plot(txx(1:100));
txx = spreadData_m;
%%
%channel_type = input('请输入信道类型：'); 
channel_type=2;
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
EbNo = 1:1:10;
for snr = 1:length(EbNo)
   y_snr = awgn(c_out,snr,'measured');

%%
    %--------------------省略了解调-------------------
    rx = y_snr;

%%
    %------------调制前导序列和用于扩频的m序列
%     m_modu = mo_bpsk(m2,v,Tm_sample);    %调制前导序列
%     m_despread = mo_bpsk(pt,v,Tm_sample);  %调制用于扩频的m序列，去软解扩
    %plot(m_modu);
    %------------相关峰捕获
    xcor = dsp.Crosscorrelator;
    buHuo = step(xcor,rx',m2'); %computes cross-correlation of x1 and x2
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
    xcor = dsp.Crosscorrelator;
    y1 = step(xcor,spreadData',pp1');
    y1Origin = abs(y1);
    y0 = step(xcor,spreadData',pp0');
    y0Origin = abs(y0);
    xcor = dsp.Crosscorrelator;
    y3 = step(xcor,data_out',pp1');
    y1Rake = abs(y3);
    y4 = step(xcor,data_out',pp0');
    y0Rake = abs(y4);
%     figure
%     subplot(411);plot(y1Origin);title('用扩1的walsh');
%     subplot(412);plot(y0Origin);title('用扩-1的walsh');
%     subplot(413);plot(y1Rake);title('矫正后，用扩1的walsh相关');
%     subplot(414);plot(y0Rake);title('矫正后，用扩-1的walsh相关');
    
    data_out_norake = rx(L_m+1:end);
    y5 = step(xcor,data_out_norake',pp1');
    y1NoRake = abs(y5);
    y6 = step(xcor,data_out_norake',pp0');
    y0NoRake = abs(y6);
%     figure
%     subplot(411);plot(y1Origin);title('用扩1的walsh');
%     subplot(412);plot(y0Origin);title('用扩-1的walsh');
%     subplot(413);plot(y1NoRake);title('过信道，用扩1的walsh相关');
%     subplot(414);plot(y0NoRake);title('过信道，用扩-1的walsh相关');
    
    %找每个码元内的最大值
    receiveJudge = ones(1,Tlen);
    tt = 1;
    sh = 5;
    for i = T_sample-sh : T_sample : L_s_initial1-T_sample-sh+1 
         if(sum(y1Rake(i : i+2*sh)) > sum(y0Rake(i : i+2*sh)))
             receiveJudge(tt) = 1;%对为-1的码元扩频
         else
             receiveJudge(tt) = -1;
         end
         tt = tt+1;
    end
    BitErrorRake = length(find(receiveJudge ~= s_initial));
     ErrorRate(snr) = BitErrorRake/(Tlen);
    %对比不同步不矫正
    receiveJudgeNoRake = ones(1,Tlen);
    tt = 1;
    for i = T_sample-sh : T_sample : L_s_initial1-T_sample-sh+1 
         if(sum(y1NoRake(i : i+2*sh)) > sum(y0NoRake(i : i+2*sh)))
             receiveJudgeNoRake(tt) = 1;%对为-1的码元扩频
         else
             receiveJudgeNoRake(tt) = -1;
         end
         tt = tt+1;
    end
    BitErrorNoRake = length(find(receiveJudgeNoRake ~= s_initial));
    ErrorRateNoRake(snr) = BitErrorNoRake/(Tlen);    
end
    figure
    subplot(411);plot(y1Origin);title('用扩1的walsh');
    subplot(412);plot(y0Origin);title('用扩-1的walsh');
    subplot(413);plot(y1Rake);title('矫正后，用扩1的walsh相关');
    subplot(414);plot(y0Rake);title('矫正后，用扩-1的walsh相关');
    figure
    subplot(411);plot(y1Origin);title('用扩1的walsh');
    subplot(412);plot(y0Origin);title('用扩-1的walsh');
    subplot(413);plot(y1NoRake);title('过信道，用扩1的walsh相关');
    subplot(414);plot(y0NoRake);title('过信道，用扩-1的walsh相关');



    