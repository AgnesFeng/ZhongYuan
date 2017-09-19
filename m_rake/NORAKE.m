clear all;
close all;
clc;
%-----------------------1、原始信号-----------------------------------------
Tlen = 800; %数据长度  编码后是三倍加9的关系
s_initial = randsrc( 1, Tlen );
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%编码
% R = 2/3;              % code rate
% g{1} = [1 0 1 1];                % Impulse Responses _ 1
% g{2} = [1 1 0 1];                % Impulse Responses _ 2
% g{3} = [1 1 1 1];                % Impluse Responses _ 3
% %n = length(g);                       % Convolution Code (1/n) parameter,n=3
% memory_els = 3;%K=4寄存器个数,m=K-1,约束长度
% initial = round(rand(1,Tlen));            % message vector信息位长度：L
% %ENCODER: 1/3 Convolution Encoder%
% s_initial = encode_1_3(initial,g,length(g));%m是信息位(1*L),g是生成矩阵，上面已经给出，n是码率的倒数。
% for i = 1:length(s_initial)%变双极性
%     if s_initial(i) == 0;
%         s_initial(i) = -1;
%     end
% end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%编码
v = 225000; %要求应该是225KBPS,对应4.4微秒，此处是0.004微秒，所以将延时降低三个数量级，0,001微秒~0.003微秒
Tb = 1/v; %码元持续时间    (最好改成真实数据。225000)

%-----------------------4、调制--------------------------------------------
%降采样为无持续时间的
downData= zeros(1,length(spreadData)/Tm_sample);
j = 1;
for i=1:Tm_sample:length(spreadData)-9
    downData(j) = spreadData(i);
    j = j+1;
end
spread = s_initial;
%串并转换
fc=300e6; %载波
len_spread = length(spread);
tx = ones(1,len_spread/2);%串并转换后的数组
tx_re = tx; 
tx_im = tx;
k2 = 1;
for k1=1:2:len_spread-1 %两两一组 ， 得到tx_re,im_re,长度16*8000 = 64000
    tx_re(k2) = spread(k1);
    tx_im(k2) = spread(k1+1);
    k2 = k2+1;
end
I = tx_re;
Q = tx_im;
%插零过采样
supersam=5;         %对IFFT而言相当于将高频短信号置0，是低通的一部分吧。
data = length(tx);  
nnn = supersam*data;  % 采样点数目=过采样率*原码数目
for  ii=1:nnn     
    if rem(ii,supersam)==1        %取余数
        tem1 = fix((ii-1)/supersam)+1;  %fix:靠0方向取整
        tem2 = fix((ii-1)/supersam)+1;
        Iinit(ii)=I(tem1);  %Iinit长度应该等于data
        Qinit(ii)=Q(tem2);
    else
         Iinit(ii)=0;
         Qinit(ii)=0;
    end
end
%脉冲成形
NT=50;
N=2*supersam*NT;    % N=500
fs=1500e6; %采样率 为采样点数乘码速率
rf=0.2;
psf=rcosfir(rf,NT,supersam,fs,'sqrt');% psf大小为500
Ipulse=conv(Iinit,psf);
Qpulse=conv(Qinit,psf);
%频谱搬移
for i=1:supersam*data+N   %采样点数目改变 （因为卷积的缘故）
    t(i)=(i-1)/(fs);      %用载频fc乘以过采样率（每符号的采样点数）=采样率。
    Imod(i)=Ipulse(i)*sqrt(2)*cos(2*pi*fc*t(i));  %点数不变，但两点之间的距离代表了采样间隔：1/fs
    Qmod(i)=Qpulse(i)*(-sqrt(2)*sin(2*pi*fc*t(i)));
end
txx=Imod+Qmod;
%-----------------------4~5、过信道----------------------------------------
sample = 1/fs;         %采样频率100
P_OR = cha3_mrake(sample,length(txx)); %（采样频率，序列长度，序列）变三径
pp1 = P_OR(1,:); %每一阶倍数   length(rx)个点
pp2 = P_OR(2,:);
pp3 = P_OR(3,:); 
EbN0db  = 0:1:15;
for snr = 1:length( EbN0db ) 
   rx = awgn(txx,snr,'measured'); 
k11 = 1;
    tap1 = pp1.* rx;
    tap2 = pp2.* rx;
    tap3 = pp3.* rx;   

%-----------------------5、解调--------------------------------------------
len_spread = length(spread);
receive1 = demod_nojudge(len_spread,data,t,tap1);
receive2 = demod_nojudge(len_spread,data,t,tap2);
receive3 = demod_nojudge(len_spread,data,t,tap3);
%还原成有持续时间的码片,使重复
r1 = REPEAT(receive1,Tm_sample); %10
r2 = REPEAT(receive2,Tm_sample);
r3 = REPEAT(receive3,Tm_sample);
%--------------

    %--------------------解调----------------------%     注意，每径信号都要解调
    tap = tap1 + tap2 + tap3;
    receive = demod(len_spread,data,t,tap) ;
   
    m =0;m1=0;m2=0;m3=0;
    for i= 1:len_spread
          if receive(i)~= spread(i)
            m= m+1;
          end
    end
    ber(k11) = m/len_spread;
    k11 =k11 + 1; %将每一个信噪比对应的误码率写进一个矩阵
end
figure
semilogy(0:15,ber,'-m+');hold on;


