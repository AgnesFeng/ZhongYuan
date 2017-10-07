clear all;
close all;

Numusers = 1;
Nc = 16; %扩频因子
ISI_Length = 1;%每径延时为ISI_Length/2
EbN0db = [0:1:15];%信噪比，单位db
Tlen = 8000; %数据长度 （取越多越好）100000
Bit_Error_Number1 = 0;%误比特率的初始值
Bit_Error_Number2 = 0;
Bit_Error_Number3 = 0;
power_unitary_factor1 = sqrt( 5/9 ); %每径功率因子
power_unitary_factor2 = sqrt( 3/9 );
power_unitary_factor3 = sqrt( 1/9 );
s_initial = randsrc( 1, Tlen ); %数据源 (1*8000的矩阵，-1,1各占一半)

% 给码元持续时间，方便画功率密度谱函数图
v = 112500;
Ts = 1; % 码元周期(画图用)
N_sample = 100; % 单个码元抽样点数
dt = Ts / N_sample; % 抽样时间间隔
N = 50; % 码元数
t = 0 : dt : (N * N_sample - 1) * dt; % 序列传输时间
gt1 = ones(1, N_sample); % 不归零波形
% 生成随机序列
RAN = round(s_initial); % 随机0 1序列
selq = [];
for i = 1 : N % 生成序列
   if RAN(i)==1
       selq = [selq gt1]; %让100个点都是1，表示一个码元
   else
       selq = [selq -1*gt1];%让100个点都是-1，表示一个码元
   end
end
figure(111)
subplot(2, 1, 1);plot(t, selq);grid on;title('原始序列时域波形');
axis([0 10 -2 2])
t3=t;after=selq;
% 功率谱密度计算方法一
fft_se1 = fftshift(fft(selq)); % 求序列的频谱
PE1 = 10 * log10(abs(fft_se1) .^ 2 / (N * Ts)); % 公式法求概率谱密度
PEL1 = (-length(fft_se1) / 2 : length(fft_se1) / 2 - 1) / 10; % 求区间长度
%%绘制出结果
figure(222)
subplot(2, 1, 1);plot(PEL1, PE1); grid on; title('原始序列PSD');
axis([-50 50 -50 100]);xlabel('f');ylabel('幅度');

%%%%%%%%%%%%%%%%%%%%%%%%产生walsh矩阵%%%%%%%%%%%%%%%%%%%%%%
wal2 = [ 1 1; 1 -1 ];
wal4 = [wal2 wal2; wal2 wal2*(-1)];  
wal8 = [wal4 wal4; wal4 wal4*(-1)];  %8*8
wal16 = [wal8 wal8; wal8 wal8*(-1)]; %16*16
%%%%%%%%%%%%%%%%%%%%%%%%%%%扩频%%%%%%%%%%%%%%%%%%%%%%%%%%%%
s_spread = zeros( Numusers, Tlen*Nc );  %[1,8000*16]扩频后的码片 初始化
ray1 = zeros( Numusers, 2*Tlen*Nc ); 
ray2 = zeros( Numusers, 2*Tlen*Nc );
ray3 = zeros( Numusers, 2*Tlen*Nc );

for i = 1:Numusers
    x0 = s_initial( i,: ).'*wal16( 8,: );  %经过转置，行变列，乘上8*8矩阵的第八行
    x1 = x0.';%变成16*8000的矩阵
    s_spread( i,: ) = ( x1(:) ).';%得到扩频后的序列，【x1(:)从头到尾排序，即16*8000个】    到时候画图用s_spread,给码元周期，画功率谱图
end

% 给码元持续时间，方便画功率密度谱函数图
Tc = Ts/16; % 码片周期
dt = Tc / N_sample; % 抽样时间间隔
N = N*16; % 码片数
t = 0 : dt : (N * N_sample - 1) * dt; % 序列传输时间
gt1 = ones(1, N_sample); % 不归零波形
% 生成随机序列
RAN = round(s_initial); % 随机0 1序列
selq = [];
for i = 1 : N % 生成序列
   if RAN(i)==1
       selq = [selq gt1]; %让100个点都是1，表示一个码元
   else
       selq = [selq -1*gt1];%让100个点都是-1，表示一个码元
   end
end
figure(111)
subplot(2, 1, 2);plot(t, selq);grid on;title('扩频后序列时域波形');xlabel('t');ylabel('幅度');
axis([0 10 -2 2]);
t2=t; before=selq;
N2 = 800;
RAN2 = round(s_spread); % 随机0 1序列
selq2 = [];
for i = 1 : N2 % 生成序列
   if RAN2(i)==1
       selq2 = [selq2 gt1]; %让100个点都是1，表示一个码元
   else
       selq2 = [selq2 -1*gt1];%让100个点都是-1，表示一个码元
   end
end
% 功率谱密度计算方法一
fft_se12 = fftshift(fft(selq2)); % 求序列的频谱
PE12 = 10 * log10(abs(fft_se12) .^ 2 / (N2 * Ts)); % 公式法求概率谱密度
PEL12 = (-length(fft_se12) / 2 : length(fft_se12) / 2 - 1) / 10; % 求区间长度
%%绘制出结果
figure(222);subplot(2, 1, 2);plot(PEL12, PE12); grid on; title('扩频后功率密度谱');
axis([-1000 1000 -350 70]);xlabel('f');ylabel('幅度');

figure(2111)
subplot(211);selq1 = SAM(s_initial,N_sample,N);%axis([-1000 1000 -350 70])
sing1=PSD2(1024,N_sample,selq1,N,Ts);title('原序列单码元功率谱');
subplot(212);selq2 = SAM(s_spread,N_sample,32*N);%axis([-1000 1000 -350 70])
sing2=PSD2(1024,N_sample,selq2,N,Ts/16);title('扩频后单码元功率谱');


%将每个扩频后的输出重复两次，有利于下面的延迟（延迟半个码元）
ray1( 1:2:2*Tlen*Nc - 1 ) = s_spread( 1:Tlen*Nc );         %相邻的即是重复的
ray1( 2:2:2*Tlen*Nc ) = ray1( 1:2:2*Tlen*Nc - 1 );         %ray1是第一径信号
spread = ray1;                                            %扩频后的发送信号
%--------------------调制----------------------%
%（省略了不同的信噪比）
fc=5e6; %载波频率
len_spread = length(spread);
tx = ones(1,len_spread/2);%调制后的数组
% receive1 = ones(1,len_spread);%解调后的接收序列   ，%避免分配空间给函数后的输出
% receive2 = ones(1,len_spread);
% receive3 = ones(1,len_spread);
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
% zero insertion，（零插入，方便直观画图）此过程称为成形。成形的意思就是实现由消息到波形的转换，以便发射，脉冲“成形”应该是在基带调制之后。
supersam=5;         %sampling  rate  25M HZ  ,supersam为过采样率。它等于 采样率fs/码速率。
data = length(tx);   % data=128000
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
%脉冲成形滤波器， 接着，将进行低通滤波，因为 随着传输速率的增大，基带脉冲的频谱将变宽
%如果不滤波（如升余弦滤波）进行低通滤波，后面加载频的时候可能会出现困难。
%平方根升余弦滤波器
%psf=rcosfir(rf,n_t,rate,fs,'sqrt')   rate:过采样率，rf:滚降因子，n_t:滤波器阶数，fs:采样率
%用在调制或发送之前，用在解调或接受之后，用来降低过采样符号流带宽并不引发ISI（码间串扰）
NT=50;
N=2*supersam*NT;    % N=500
fs=25e6;
rf=0.3;
psf=rcosfir(rf,NT,supersam,fs,'sqrt');% psf大小为500
Ipulse=conv(Iinit,psf);
Qpulse=conv(Qinit,psf);
%modulation
for i=1:supersam*data+N   %采样点数目改变 （因为卷积的缘故）
    t(i)=(i-1)/(fs);      %这里因为假设载频与码速率大小相等，所以用载频fc乘以过采样率=采样率。
    Imod(i)=Ipulse(i)*sqrt(2)*cos(2*pi*fc*t(i));
    Qmod(i)=Qpulse(i)*(-sqrt(2)*sin(2*pi*fc*t(i)));
end
tx=Imod+Qmod;
%画出成形的波
% x1 = 0:1:4999;
% figure
% txx = tx(1:5000);      %节选5000个点
% plot(x1,txx)
% xlabel('节选前5000点');ylabel('信号幅度')
% title('调制后的成形波')
%    figure
%    plot(20*log(abs(fft(s_initial))));
%    axis([0  data  -40  100]);
%    grid on;

ttt = 0:length(tx);
nfft = 1024;
cxn = xcorr(s_initial,'unbiased');%计算序列的自相关函数
CXk = fftshift(fft(cxn,nfft));
Pxx = abs(CXk);
%index = 0:round(nfft/2-1);
%k = index*1/nfft;
%plotPxx = 10*log10(Pxx(index+1));
in = (-length(CXk) / 2 : length(CXk) / 2 - 1) / 10; % 求区间长度
plotPxx = 10*log10(Pxx);
% figure(111)
% subplot(313)
% plot(in,plotPxx);%变平滑一点就好？
% xlabel('频率');ylabel('幅度');
% title('第一径调制后功率密度谱函数')


%-------------------过信道----------------------%
sample = 1/fs;         %采样频率100
P_OR = cha_3(sample,length(tx),tx);          %（采样频率，序列长度，序列）变三径
pp1 = abs(P_OR(1,:)); %每一阶倍数   length(rx)个点
pp2 = abs(P_OR(2,:));
pp3 = abs(P_OR(3,:)); 
k11 = 1;
for snr = 1:length( EbN0db ) 
    rx = awgn(tx,snr,'measured');            %高斯信道(只考虑高斯噪声) 加性噪声
    tap1 = pp1.* rx;
    tap2 = pp2.* rx;
    tap3 = pp3.* rx;   
%【【【【【【【【【【单个信噪比时画图用-----------------------
% tt1 = tap1(1:5000);      %节选5000个点
% tt2 = tap2(1:5000);
% tt3 = tap3(1:5000);
% figure
% subplot(311)
% plot(x1,tt1);
% ylabel('第1径信号幅度')
% title('信号经过衰落后的传送序列');
% subplot(312)
% plot(x1,tt2);ylabel('第2径信号幅度')
% subplot(313)
% plot(x1,tt3);xlabel('前5000个点');ylabel('第3径信号幅度')
%--------------------单个信噪比时画图用】】】】】】】】】】】】

    %--------------------解调----------------------%     注意，每径信号都要解调
    receive1 = demod(len_spread,data,t,tap1) ;
    receive_2 = demod(len_spread,data,t,tap2) ;
    receive2( ISI_Length + 1:len_spread) = receive_2( 1:len_spread - ISI_Length ); %…对第一径信号进行一次延迟，模拟延时后的第二径
    %receive2 = receive_2;%毕设时用的上面一句
    
    receive_3 = demod(len_spread,data,t,tap3) ;
    receive_3( ISI_Length + 1:len_spread ) = receive_3( 1:len_spread - ISI_Length );%…对第一径信号进行两次延迟，模拟延时后的第二径和第三经
    receive3( 2*ISI_Length + 1:len_spread ) = receive_3( 1:len_spread - 2*ISI_Length );
    %receive3 = receive_3;%毕设时用的上面一句,误码率反而变大了，说明接收端延时对齐的不对
    
%     m1=0;m2=0;m3=0;             %算三径的误码率
%     for i= 1:len_spread
%       if receive1(i)~= spread(i)
%         m1= m1+1;
%       end
%       if receive2(i)~= spread(i)
%         m2= m2+1;
%       end
%       if receive3(i)~= spread(i)
%         m3= m3+1;
%       end
%     end
%     ber1(k11) = m1/len_spread;  %针对一个信噪比的误码率
%     ber2(k11) = m2/len_spread;
%     ber3(k11) = m3/len_spread;
%     k11 =k11 + 1; %将每一个信噪比对应的误码率写进一个矩阵

% figure
% semilogy(0:10,ber1,'-b*');grid on; hold on;
% semilogy(0:10,ber2,'-g*');hold on;
% semilogy(0:10,ber3,'-r*');hold on;

ttt = 0:length(receive1);
nfft = 1024;
cxn = xcorr(receive1,'unbiased');%计算序列的自相关函数
CXk = fftshift(fft(cxn,nfft));
Pxx = abs(CXk);
%index = 0:round(nfft/2-1);
%k = index*1/nfft;
%plotPxx = 10*log10(Pxx(index+1));
in = (-length(CXk) / 2 : length(CXk) / 2 - 1) / 10; % 求区间长度
plotPxx = 10*log10(Pxx);
% figure(2)
% subplot(211)
% plot(in,plotPxx);%变平滑一点就好？
% xlabel('频率');ylabel('幅度');
% title('第一径解调后功率密度谱函数')

%--------------------解扩----------------------
%for nEN = 1:length( EbN0db )    
    en = 10^( EbN0db(snr)/10 ); %将Eb/N0的db值转化成十进制数值
    sigma = sqrt(32/(2*en));
    %sigma = sqrt( 32/(2*en) );%根据初始状态设定的信噪比，计算出AWGN方差，32为比特能量，与扩频调制方式有关，扩频后每径32bit对应一个原始bit，比特能能量是原来的32倍
    %接收到的信号demp
    demp = power_unitary_factor1*receive1+...           %每径乘功率因子再叠加
    power_unitary_factor2*receive2+...                  %注意下一行rand和randn的区别
    power_unitary_factor3*receive3 + ( 0.2*rand( 1,len_spread )+ 0.8*randn( 1,len_spread )*Numusers )*sigma;          %最后一项为噪声     %rand:0~1之间随机数             randn:正态分布的随机数
    dt = reshape( demp,32,Tlen )';   %dt为8000*32的demp值的矩阵
    %将walsh码重复两次
    wal16_d(1:2:31) = wal16(8,1:16); %用刚才扩频的序列填充
    wal16_d(2:2:32) = wal16(8,1:16);
    
    %解扩后rdata1为第一径输出
    rdata1 = dt*wal16_d(1,:).';  
    %将walsh码延迟半个码片
    wal16_delay1(1,2:32) = wal16_d(1,1:31);
    wal16_delay1(1) = wal16_d(32);%自己加的，避免有0，尔出现三但结果偶条线交叉，其实前面是不是0无所谓，因为取得后面的数相乘，后来发现不是去后面的数，而是算个总功率，
    %解扩后rdata2为第二径输出
    rdata2 = dt*wal16_delay1(1,:).';
    %将walsh码延迟一个码片
    wal16_delay2(1,3:32) = wal16_d(1,1:30);%（矩阵的第一个下标为1不是0）
    wal16_delay2(1,1:2) = wal16_d(1,31:32);%（最后两个值放到最前，避免有0）
    %解扩后rdata3为第三径输出
    rdata3 = dt*wal16_delay2(1,:).';
    
    p1 = rdata1'*rdata1;   %计算每一径的功率
    p2 = rdata2'*rdata2;
    p3 = rdata3'*rdata3;
    p = p1 + p2 + p3;
    u1 = p1/p;             %算每一径的加权系数
    u2 = p2/p;
    u3 = p3/p; 
    %最大比合并
    rd_m1 = real( rdata1*u1+rdata2*u2+rdata3*u3);
    %等增益合并
    rd_m2 = (real(rdata1+rdata2+rdata3))/3;
    %选择式合并
    u = [u1,u2,u3];   %（循环检测）
    maxu = max(u);
    if(maxu==u1)
        rd_m3 = real(rdata1);
      else if(maxu==u2)
           rd_m3 = real(rdata2);
      else
           rd_m3 = real(rdata3);
      end
    end %三种方法判决输出
    r_Data1 = sign(rd_m1)'; %大于0为1，小于0为-1
    r_Data2 = sign(rd_m2)';
    r_Data3 = sign(rd_m3)';
    %计算误比特率
    Bit_Error_Number1 = length(find(r_Data1(1:Tlen) ~= s_initial(1:Tlen)));     %与原始信号比较
    Bit_Error_Rate1(snr) = Bit_Error_Number1/Tlen;
    Bit_Error_Number2 = length(find(r_Data2(1:Tlen) ~= s_initial(1:Tlen)));
    Bit_Error_Rate2(snr) = Bit_Error_Number2/Tlen;
    Bit_Error_Number3 = length(find(r_Data3(1:Tlen) ~= s_initial(1:Tlen)));
    Bit_Error_Rate3(snr) = Bit_Error_Number3/Tlen;

end
    
figure(3333)
subplot(2, 1, 1);plot(t2, before);grid on;title('判决后序列时域波形');
axis([0 20 -2 2]);
subplot(2, 1, 2);plot(t3, after);grid on;title('解扩后序列时域波形');
axis([0 20 -2 2]);
 
figure
semilogy(EbN0db,Bit_Error_Rate1,'ro-');hold on; 
semilogy(EbN0db,Bit_Error_Rate2,'bo-');hold on;
semilogy(EbN0db,Bit_Error_Rate3,'go-');hold on;
legend('最大比合并','等增益合并','选择式合并');
xlabel('信噪比');ylabel('误信率');
title('三种主要分集合并方式性能比较');
grid  on;