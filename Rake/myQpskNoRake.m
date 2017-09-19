
Tlen = 2000;
s_spread = randsrc( 1, Tlen ); %数据源 (1*8000的矩阵，-1,1各占一半)
%--------------------1采样,前50个脉冲----------------------%
N_sample=100; %画图时采样率大一点好，是码速率的100倍，后面用滤波时，采样为5倍就可以
N = 50;
selq = SAM(s_spread,N_sample,N);%这个函数是画图用的，采样点多，后面和滤波器相乘时用的另一种方法（插0），采样点数也较少

%--------------------2分成两通道----------------------%
%调制之后再采样效果跟采样后调制一样，但是不单纯为画图了，而是方便与载波相乘，注意相乘时采样点数要一样
%（省略了不同的信噪比）
len_spread = length(s_spread);
receive1 = ones(1,len_spread);%解调后的接收序列
tx_re = ones(1,len_spread/2); 
tx_im = ones(1,len_spread/2);
k2 = 1;
for k1=1:2:len_spread-1 %两两一组  
    tx_re(k2) = s_spread(k1);
    tx_im(k2) = s_spread(k1+1);
    k2 = k2+1;
end
selq_re = SAM(tx_re,N_sample,Tlen/2);
selq_im = SAM(tx_im,N_sample,Tlen/2);
figure(1)
subplot(211);PSD222(1024,N_sample,selq,N);title('原序列功率谱');%这里画图的功率谱指的都是密度，功率谱是功率，单位是瓦；功率谱密度是能量，焦耳
subplot(212);PSD222(1024,N_sample,selq_re,N);title('2进制变4进制后I通道功率谱'); 

%---------------------脉冲成形滤波，同时当低通。1、过采样；2、构造滤波器；3、卷积-------------%
%1、过采样
zero=5;         
I=tx_re;
Q=tx_im;
Izero = SymbolToWaveform(I,zero);
Qzero = SymbolToWaveform(Q,zero);

%2、构造滤波器
NT=50;           %滤波器阶数为2*NT+1,根据结果来调
N_fir=2*NT*zero;    %滤波器长度-1，与length(psf)-1一样
fc = 5e6;          %载波频率
fs = 25e6;         %采样频率
Ts = 1e-6;         %码元长度（SAM函数里）
rf=0.1;
psf=rcosfir(rf,NT,zero,Ts,'sqrt');% 滚降因子rf、滤波器阶数NT、过采样率（单个脉冲的采样点数）、码元时间
Ipulse=conv(Izero,psf);%卷积后长度改变
Qpulse=conv(Qzero,psf);
figure(2)
subplot(321);plot(psf);title('滤波器时域响应');
axis([200  300  -0.2  0.6]);grid  on;
set(gca,'XTick',0:250:501);set(gca,'XTicklabel',{'-250','0','250'});
subplot(322);plot(20*log(abs(fftshift(fft(psf)))));title('滤波器频域响应');
axis([0  N_fir  -350 50]);grid on;
set(gca,'XTick',0:250:500);set(gca,'XTicklabel',{'-250','0','250'});
subplot(323);plot(I(1:100));axis([0  100  -2 2]);title('滤波前单通道时域响应');
subplot(324);plot(Ipulse(1:600));title('滤波之后单通道时域响应');
subplot(325);plot(20*log(abs(fftshift(fft(I)))));title('滤波前单通道频域响应');
axis([0 1000 -50 100]);grid  on;
set(gca,'XTick',0:500:1000);set(gca,'XTicklabel',{'-500','0','500'});
subplot(326);plot(20*log(abs(fftshift(fft(Ipulse)))));title('滤波之后单通道频域响应');
axis([0 5500 -300 100]);grid  on;
set(gca,'XTick',0:2750:5500);set(gca,'XTicklabel',{'-2700','0','2700'});
%3、乘载波
dt = 1/fs;         %时间采样间隔  
T = 1;            %信号时长                 
tt = 0:dt:T;
t1 = tt(1:length(Ipulse));
cos1 = cos(2*pi*fc*t1); 
sin1 = sin(2*pi*fc*t1);
Imod=Ipulse.*(sqrt(2)*cos1);
Qmod=Qpulse.*(sqrt(2)*sin1);
sum=Imod+Qmod;

figure(3)%画功率谱密度
plot(20*log(abs(fftshift(fft(sum)))));title('调制后双通道能量谱(加成形和滤波)');
axis([0 5500 -200 100]);grid  on;
set(gca,'XTick',0:2750:5500);set(gca,'XTicklabel',{'-2700','0','2700'});
%[Pxx,f]=psd(sum,1024,fs,window,noverlap,dflag);画出的频谱不好看

%--------------------乘上正交的载波(无滤波)---------------------%
% t2 = tt(1:length(selq_re));
% cos2 = cos(2*pi*fc*t2); 
% sin2 = -sin(2*pi*fc*t2);
% tx_re1 = selq_re .* cos2;   %不加滤波
% tx_im1 = selq_im .* sin2;
% tx1 = tx_re1 + tx_im1;
%--------------------------------------------------------------------%

% %-------------------过信道----------------------%
P_OR = cha_3(1/Ts,length(sum),sum);          %（采样频率，序列长度，序列）变三径
pp1 = abs(P_OR(1,:)); %每一阶倍数   length(rx)个点 小尺度加大尺度的结果
pp2 = abs(P_OR(2,:));
pp3 = abs(P_OR(3,:));
k11 = 1;
j=1;
snrn = 15;
for snr = 0:snrn
     %rx = awgn(tx1,snr,'measured'); %高斯信道(只考虑高斯噪声) 加性噪声
     rx = awgn(sum,snr,'measured'); %高斯信道(只考虑高斯噪声) 加性噪声 
     tap1 = pp1.* rx;
     tap2 = pp2.* rx;
     tap3 = pp3.* rx; 
%     %（（（（（（（（--解调（无滤波）----------------------%    
%     re_2=rx.* (sqrt(2)*cos2);
%     im_2=rx.* (sqrt(2)*sin2);
%     re_21=rx.* re_2;
%     im_21=rx.* im_2;
%     for k3 = 1:length(rx)
%         if re_2(k3) > 0;
%             re_2(k3) = 1;
%         else
%             re_2(k3) = -1;
%         end
%         if im_2(k3) > 0;
%             im_2(k3) = 1;
%         else
%             im_2(k3) = -1;
%         end
%     end
%        %抽样
%        for i=1:length(I)                       
%        I2(i)=re_2((i-1)*N_sample+1);
%        Q2(i)=im_2((i-1)*N_sample+1);
%        end
%     k4 = 1;
%     for k1=1:2:len_spread-1 %实部虚部放到一个数组receive中  
%         receive(k1) = I2(k4);
%         receive(k1+1) = Q2(k4);
%         k4 = k4+1;
%     end
%     %-----------采样并画图-------------------------%
%     selq_rece = SAM(receive,N_sample,N);
%     % figure(2)
%     subplot(211);PSD222(1024,N_sample,tx1,Tlen/2);
%     axis([-4000, 4000, -60, 20]);title('调制后信号功率谱(未加滤波)')
%     subplot(212);PSD222(1024,N_sample,selq_rece,N);
%     axis([-4000, 4000, -60, 20]);title('解调后功率谱（不加滤波的）')
%     m=0;
%     for i= 1:length(I)
%       if (receive~=s_spread)
%         m= m+1;
%       end
%     end
%     ps(j) = m/length(receive)        %针对一个信噪比的误码率
%     j =j+1;                          %将每一个信噪比对应的误码率写进一个矩阵
% end
%--------------------------------------））））））））））））））%   

%（（（（（--------------------解调（有滤波）----------------------%   
% %频谱搬移
data = length(I);
Idem=tap1.*(sqrt(2)*cos1);
Qdem=tap1.*(sqrt(2)*sin1);

%过滤波
Imat=conv(Idem,psf);
Qmat=conv(Qdem,psf);   
supersam=zero;
   for  i=1:supersam*data                 %滤波后抽样
       Isel(i)=Imat(i+N_fir);
       Qsel(i)=Qmat(i+N_fir);
   end
          
   for i=1:data                       %提取码元,去掉0
       Isam(i)=Isel((i-1)*supersam+1);
       Qsam(i)=Qsel((i-1)*supersam+1);
   end
   
   threshold=0.2;                     %判决门限
   for  i=1:data
       if Isam(i)>=threshold
           Ifinal(i)=1;
       else
           Ifinal(i)=-1;
       end
       if Qsam(i)>=threshold
           Qfinal(i)=1;
       else
           Qfinal(i)=-1;
       end
   end
    k4 = 1;
    for k1=1:2:len_spread-1 %实部虚部放到一个数组receive中  
        receive1(k1) = Ifinal(k4);
        receive1(k1+1) = Qfinal(k4);
        k4 = k4+1;
    end
    %------------------------第二径第三经的解调-----------------------%
    tap = tap1 + tap2 + tap3;
    receive_sum = demod2(psf,cos1,sin1,Tlen,Tlen/2,tap);
    receive2 = demod2(psf,cos1,sin1,Tlen,Tlen/2,tap2); %序列原长度,分通道后长度,时间轴，输入序列
    receive3 = demod2(psf,cos1,sin1,Tlen,Tlen/2,tap3); %序列原长度,分通道后长度,时间轴，输入序列

    m1=0;m2=0;m3=0;m =0;
    for i= 1:len_spread
          if receive_sum(i)~= s_spread(i)
            m= m+1;
          end
            if receive1(i)~= s_spread(i)
            m1= m1+1;
          end
          if receive2(i)~= s_spread(i)
            m2= m2+1;
          end
          if receive3(i)~= s_spread(i)
            m3= m3+1;
          end
    end
    ber(k11) = m/len_spread;
    ber1(k11) = m1/len_spread;%针对一个信噪比的误码率
    ber2(k11) = m2/len_spread;
    ber3(k11) = m3/len_spread;
    k11 =k11 + 1;    
end

figure
%semilogy(0:15,berno,'-b*');grid on; hold on;
semilogy(0:snrn,ber,'-m*');grid on; hold on;
semilogy(0:snrn,ber1,'-b*');grid on; hold on;
semilogy(0:snrn,ber2,'-go');hold on;
semilogy(0:snrn,ber3,'-r*');hold on;
legend('三径叠加','第一径','第二径','第三经');
xlabel('信噪比');ylabel('误信率');
figure
subplot(321);
plot(20*log(abs(fftshift(fft(Imod)))));title('调制后单通道能量谱');
axis([0 5500 -200 100]);grid  on;
set(gca,'XTick',0:2750:5500);set(gca,'XTicklabel',{'-2700','0','2700'});
subplot(322);
plot(20*log(abs(fftshift(fft(Idem)))));title('解调--频谱搬移--能量谱');
axis([0 5500 -50 100]);grid  on;
set(gca,'XTick',0:2750:5500);set(gca,'XTicklabel',{'-2700','0','2700'});
subplot(323);
plot(20*log(abs(fftshift(fft(Imat)))));title('解调--过滤波器--频域响应');
axis([0 6000 -300 130]);grid  on;
set(gca,'XTick',0:3000:6000);set(gca,'XTicklabel',{'-300','0','3000'});
subplot(324);
plot(20*log(abs(fftshift(fft(Isam)))));title('抽样判决后频域响应');
axis([0 1000 -50 130]);grid  on;
set(gca,'XTick',0:500:1000);set(gca,'XTicklabel',{'-500','0','500'});
subplot(325);plot(Imat(1:800));title('解调--过滤波器--时域响应');
subplot(326);plot(Ifinal(1:100));axis([0  100  -2 2]);title('抽样判决之后单通道时域响应');


