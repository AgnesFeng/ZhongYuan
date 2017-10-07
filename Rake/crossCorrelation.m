clear all;
close all;
%���ʾ������R2016b���Ժ����С������ʹ�õ��ǽ���İ汾����ʹ�õ�Ч�Ĳ����﷨�滻ÿ�����á����磬myObject(x)�����step(myObject��x)
xcorr = dsp.Crosscorrelator;
t = 0:0.001:1;
x1 = sin(2*pi*2*t)+0.05*sin(2*pi*50*t);
x2 = sin(2*pi*2*t);
y = step(xcorr,x1,x2); %computes cross-correlation of x1 and x2
figure,plot(t,x1,'b',t, x2, 'g');
legend('Input signal 1',' Input signal 2');
figure, plot(y); title('Correlated output')

%%
%���ý��������������Ͼ�ֹ�İ�ɫ��˹����������ӳ�
S = rng('default');
x = randn(100,1);%����Ϊ��˹����
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
m = m_sequence([0,0,0,1,1,1,0,1]); %�õ�һ�����ڵ�m���� (���Ĵ����ĳ�ʼ״̬��������ΪC1~Cn��c0Ĭ��Ϊ1) 
N = length(m);%m��:225
v = 225000;
Tb = 1/v; %��Ԫ����ʱ��  
Tc = Tb/N;
Tm_sample = 4;%������Խ�󣬽��ͼԽ��׼��ע�⣬m�����ظ��󳤶ȿ��ܲ�����
m1 =  SAM_0(m,Tm_sample,length(m),Tc); %������ʱ��
m2 = 1-2*m1;  %��˫����
L = length(m2);
NP = 1;
mm = repmat(m,1,NP); %3�����ڵ�m���У����������
mm1 = SAM_0(mm,Tm_sample,length(mm),Tc);%�����ڣ��г���ʱ��
mm2 = 1-2*mm1;       %��˫����
dt = Tc/Tm_sample;
t=0:dt:Tc-dt;
n = 1:length(m2);
% rt1=conv(mm2(1:3*L),m2(end:-1:1))/(N*Tm_sample);%ȡ�м��һ����
% figure
% index = L+Tm_sample;
% plot(rt1(index:end-index));title('������m���е�����غ���');
xcor = dsp.Crosscorrelator;
x1 = mm2';
delay = dsp.Delay(300);
x2 = step(delay,x1); %�����������ص�ľ����ж�λ��
y = step(xcor,x2,x1); %computes cross-correlation of x1 and x2
figure(3), subplot(311);plot(y); title('Correlated output')
%%%���Լ��ķ�����ʱ������ͷ
x2 = x1';
x2(1,300+1:end) = x2(1,1:end-300);
x2 = x2';
y = step(xcor,x2,x1); %computes cross-correlation of x1 and x2
subplot(312); plot(y); title('Correlated output')
%%%�Լ��ķ�����ʱ����ͷ
x2 = x1';
x2(1,300+1:end) = x2(1,1:end-300);x2(1,1:300) = x2(1,end-300+1:end);
x2 = x2';
y = step(xcor,x2,x1); %computes cross-correlation of x1 and x2
subplot(313); plot(y); title('Correlated output')
%���ۣ���ʱ�����delay��ϸ΢���𣬵�������dsp.Crosscorrelator���ʱ������ʱ��ķ�ǰ�棬��ؽ��ֻ���ְ�ߵľ�һģһ����

