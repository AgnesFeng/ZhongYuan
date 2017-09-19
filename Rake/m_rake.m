clear all;        %仅生成一个周期的m序列画图效果不好
close all;
primpoly(3,'all');%得到所有7阶(n为7)的所有本原多项式
m = m_sequence([0,1,0,0,1]); %得到一个周期的m序列 (各寄存器的初始状态，从左到右为C1~Cn，c0默认为1) 
%m = [1 0 0 1 0 0 1 1 0 1 0 1 1 1]; %三值
%m = [1 1 1 0 0 1 0]; %是m序列
N = length(m);
v = 225;
Tb = 1/v; %码元持续时间
Tc = Tb/N;
T_sample = 310;%采样率越大，结果图越精准（注意，m序列重复后长度可能不够）
Tm_sample = T_sample/N;  %m序列的采样点数为10
m1 =  SAM_0(m,Tm_sample,length(m),Tc); %给持续时间
m2 = 1-2*m1;  %变双极性

%--------------------- 1、最原始的自相关和功率密度谱-------------------------
figure(1)
[a1,f] = xcorr(m1,'coeff');%用了coeff零频就在中间，unbaised要搬移
subplot(311);plot(f,(a1));title('m序列的自相关函数，一个周期长度，无持续时间');%
%利用Fourier变换中的卷积定理求自相关  
%Rt = fftshift(ifft(fft(m).*fft(m)));%m和m1画出来是一样的。subplot(212);plot(Rt);title('用fft求的自相关函数');%是错的

selq2 = m1;
fft_se12 = fftshift(fft(selq2));
PE12 = 10 * log10(abs(fft_se12) .^ 2 / (N * Tc)); % 公式法求概率谱密度
PEL12 = (-length(fft_se12) / 2 : length(fft_se12) / 2 - 1) / 10; % 求区间长度
%%绘制出结果
subplot(312);plot(PEL12, PE12); grid on; title('单极性m序列功率密度谱');
axis([-50 50 -50 50]);xlabel('f');ylabel('幅度');

%-------------------------2、将m序列重复NP个周期----------------------------
NP = 500;
mm = repmat(m,1,NP); %多周期
mm1 = SAM_0(mm,Tm_sample,length(mm),Tc);%多周期，有持续时间
mm2 = 1-2*mm1;       %变双极性
    %---------画频谱
    ym = fftshift(fft(mm2,4096));
    magm = abs(ym);
    fm = (1:4096)*200/2048;
    subplot(313);plot(fm,magm*2/4096);title('多个周期的31位m序列的频谱');%m序列的频谱
    %---------
L = length(mm);   %NP个周期的m的总长
dt = Tc/Tm_sample;
t=0:dt:L*Tc-dt;
rt1=conv(m2,m2(end:-1:1))/(N*Tm_sample);%最后一行到第一行换反过来,卷积,只能是一长一短
max = max(rt1);
figure(2)
%plot(t,rt1(length(m2):length(m2)+length(t)-1));title('m序列的自相关函数(一个周期，有持续时间)');
plot(rt1);
%axis([0 0.02 -0.2 1.2]);
%--------------------------[0916重新算自相关-----------------
%一个周期的m序列是m2
m2c = m2;
t = 0:dt:N*Tc-dt;
result0=m2.*m2(end:-1:1)/(N*Tm_sample);
for temp = 1:Tm_sample-1
m2(1,temp+1:end) = m2(1,1:end-temp);
m2(1,1:temp)=m2(1,end-temp+1:end);
result=m1.*m1;
result1 = sum((result))/(N*Tm_sample);
end
plot(result);

%--------------------------0916重新算自相关]-----------------

%--------------------------3、设置阈值------(应该放在扩频后的序列)-----------------
%将双极性、有持续时间的、一个周期的m序列mm2延迟Tc/2
% m2(1,Tm_sample/2+1:end) = m2(1,1:end-Tm_sample/2);
% m2(1,1:Tm_sample/2)=m2(1,end-Tm_sample/2+1:end); %m2变了，可以跟单极性的m1比是否位移
% limid=conv(mm2,m2(end:-1:1))/(N*Tm_sample);%最后一行到第一行换反过来
% figure(3),plot(t,rt1(length(m2):length(m2)+length(t)-1));title('m序列移位后的相关性');
% axis([0 200 -0.2 1.2]);
   
%-----------------------5、原始信号-------------------------------------------
% Tlen = 50; %数据长度 
% s_initial = randsrc( 1, Tlen );
% %Tb = 31*Tc; %与上面m序列的采样率没有关系，扩频时m序列先没有采样率
% s_initial1 =  SAM2(s_initial,T_sample,Tlen,Tb); 
% %----画时域特性
% dt = Tb/T_sample;   %dt是一致的
% t = 0 : dt : (Tlen * T_sample - 1) * dt; % 序列传输时间
% figure(4);subplot(211);plot(t,s_initial1);axis([0,0.1,-2,2]);title('扩频前时域');
% %----画频域特性
% selq2 = s_initial1;
% fft_se12 = fftshift(fft(selq2)); % 求序列的频谱
% PE12 = 10 * log10(abs(fft_se12) .^ 2 / (Tlen * Tb)); % 公式法求概率谱密度
% PEL12 = (-length(fft_se12) / 2 : length(fft_se12) / 2 - 1) / 10; % 求区间长度
% figure(5);subplot(2, 1, 1);plot(PEL12, PE12); grid on; title('扩频前功率密度谱');
% axis([-200 200 -50 100]);xlabel('f');ylabel('幅度');
% 
% %---------------------------6、扩频(改成了逐位相乘)------------------------------------------   
% mseq = mm2(1:length(s_initial1));
% spreadData = s_initial1.*mseq;       %扩频后的信号
% %----画时域特性
% figure(4);subplot(212);plot(t,spreadData);axis([0,0.1,-2,2]);title('扩频后时域');
% %----画频域特性
% selq2 = spreadData;
% fft_se12 = fftshift(fft(selq2)); % 求序列的频谱
% PE12 = 10 * log10(abs(fft_se12) .^ 2 / (Tlen * Tb)); % 公式法求概率谱密度
% PEL12 = (-length(fft_se12) / 2 : length(fft_se12) / 2 - 1) / 10; % 求区间长度
% figure(5);subplot(2, 1, 2);plot(PEL12, PE12); grid on; title('扩频后功率密度谱');
% axis([-200 200 -500 100]);xlabel('f');ylabel('幅度');
% 
% figure(6)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% spreadshifft1 = spreadData;    %滑动1/2码片    几乎不改变最大相关值，改成1/2码元
% spreadshifft1(1,T_sample/2+1:end) = spreadshifft1(1,1:end-T_sample/2);
% spreadshifft1(1,1:T_sample/2)=spreadshifft1(1,end-T_sample/2+1:end); %m2变了，可以跟单极性的m1比是否位移
% [a1,f] = xcorr(spreadshifft1,spreadData,'coeff');%用了coeff零频就在中间，unbaised要搬移
% subplot(312);plot(f,(a1));title('滑动1/2码元时自相关函数');axis([-500,500,-0.2,1.2])
% 
% spreadshifft2 = spreadData;    %滑动1码元   
% spreadshifft2(1,T_sample+1:end) = spreadshifft2(1,1:end-T_sample);
% spreadshifft2(1,1:T_sample)=spreadshifft2(1,end-T_sample+1:end); %m2变了，可以跟单极性的m1比是否位移
% [a2,f] = xcorr(spreadshifft2,spreadData,'coeff');%用了coeff零频就在中间，unbaised要搬移
% subplot(313),plot(f,(a2));title('滑动1码元时相关函数');axis([-500,500,-0.2,1.2])
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [a,f] = xcorr(spreadData,'coeff');%用了coeff零频就在中间，unbaised要搬移
% subplot(311);plot(f,(a));title('相位对齐时的自相关函数');
% axis([-500,500,-0.2,1.2]);
% %---------------------------6~7、频谱搬移调制----------------
% 
% 
% %---------------------------7、过信道----------------
% chanData = awgn(spreadData,15,'measured');
% %1
% %2
% %3
% %---------------------------8、解调------------------(打算在解调前抽样变成无持续时间的)----------------------
% receive = chanData;
% 
% %---------------------------9、相关峰捕获-------------
% %for tt = 0:dt:T_samople
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%已得到相位差，得到滑动时间
% %---------------------------10、解扩------------------
% %解扩码     用m序列做相关没有相关峰，要用和m序列相乘后的序列
% demseq = mseq;  
% demseq(1,T_sample+1:end) = demseq(1,1:end-T_sample);
% demseq(1,1:T_sample)=demseq(1,end-T_sample+1:end);%demseq是对齐相位后的
% demp = demseq.* spreadshifft2;
% 
% temp = s_initial1; %所以合并前要相位对齐
% temp(1,T_sample+1:end) = temp(1,1:end-T_sample);
% temp(1,1:T_sample)=temp(1,end-T_sample+1:end);%demseq是对齐相位后的
% a = demp-temp;
% [a1,f] = xcorr(demp,s_initial1,'coeff');%用了coeff零频就在中间，unbaised要搬移
% figure(8);plot(f,(a1));title('第一条经，相关性检验'); 
% 
% demseq2 = mseq;
% demp = mseq.*spreadData;
% b = demp-s_initial1;
% %下一步，加调制解调。


