clear all;
close all;
%这个示例仅在R2016b或以后运行。如果您使用的是较早的版本，则使用等效的步骤语法替换每个调用。例如，myObject(x)变成了step(myObject，x)
xcorr = dsp.Crosscorrelator;
t = 0:0.001:1;
x1 = sin(2*pi*2*t)+0.05*sin(2*pi*50*t);
x2 = sin(2*pi*2*t);
y = step(xcorr,x1,x2); %computes cross-correlation of x1 and x2
figure,plot(t,x1,'b',t, x2, 'g');
legend('Input signal 1',' Input signal 2');
figure, plot(y); title('Correlated output')

%%
%利用交叉相关来检测联合静止的白色高斯噪声输入的延迟
S = rng('default');
x = randn(100,1);%输入为高斯噪声
%Create copy delayed by 10 samples x1[n] = x[n-10]
delay = dsp.Delay(70);
x1 = step(delay,x);
xcorr = dsp.Crosscorrelator;
y = step(xcorr,x1,x);
lags = 0:99; %Positive lags
stem(lags,y(100:end),'markerfacecolor',[0 0 1]);
%axis([0 99 -125 125]);
xlabel('Lags'); 
title('Cross-Correlation of Input Noise and Delayed Version');

%%
m = m_sequence([0,0,0,1,1,1,0,1]); %得到一个周期的m序列 (各寄存器的初始状态，从左到右为C1~Cn，c0默认为1) 
N = length(m);%m长:225
v = 225000;
Tb = 1/v; %码元持续时间  
Tc = Tb/N;
Tm_sample = 4;%采样率越大，结果图越精准（注意，m序列重复后长度可能不够）
m1 =  SAM_0(m,Tm_sample,length(m),Tc); %给持续时间
m2 = 1-2*m1;  %变双极性
L = length(m2);
NP = 1;
mm = repmat(m,1,NP); %3个周期的m序列，方便做相关
mm1 = SAM_0(mm,Tm_sample,length(mm),Tc);%多周期，有持续时间
mm2 = 1-2*mm1;       %变双极性
dt = Tc/Tm_sample;
t=0:dt:Tc-dt;
n = 1:length(m2);
% rt1=conv(mm2(1:3*L),m2(end:-1:1))/(N*Tm_sample);%取中间的一部分
% figure
% index = L+Tm_sample;
% plot(rt1(index:end-index));title('单周期m序列的自相关函数');
xcor = dsp.Crosscorrelator;
x1 = mm2';
delay = dsp.Delay(300);
x2 = step(delay,x1); %根据最大点离重点的距离判断位移
y = step(xcor,x2,x1); %computes cross-correlation of x1 and x2
figure(3), subplot(311);plot(y); title('Correlated output')
%%%用自己的方法延时，不补头
x2 = x1';
x2(1,300+1:end) = x2(1,1:end-300);
x2 = x2';
y = step(xcor,x2,x1); %computes cross-correlation of x1 and x2
subplot(312); plot(y); title('Correlated output')
%%%自己的方法延时，补头
x2 = x1';
x2(1,300+1:end) = x2(1,1:end-300);x2(1,1:300) = x2(1,end-300+1:end);
x2 = x2';
y = step(xcor,x2,x1); %computes cross-correlation of x1 and x2
subplot(313); plot(y); title('Correlated output')
%结论：延时方面和delay有细微区别，但是在用dsp.Crosscorrelator相关时，将延时后的方前面，相关结果只看又半边的就一模一样。

