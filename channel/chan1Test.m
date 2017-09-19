%function PSL = chan1Test(sample,N)
clear all;
close all;
sample = 1500;        
N = 30000;
dt = 1/sample; %ʱ�������� 
T = N*dt;
t = 0:dt:T-dt;       
fm = 300;
mt = cos(2*pi*fm*t);

ts = dt;
fd = 25;           %��������Ƶ��
k = 12;
chan = ricianchan(ts,fd,k); %ts�������źŵĲ���ʱ�䣨s����fd����������Ƶ�ƣ�pdb��ƽ��·�����棨dB��
y = filter(chan,mt);       %Ϊ������ʽ
y1 = abs(y); %˥�����ŵ�ȡģ
fad = y1./abs(mt); %˥�䱶��
fad1 = 10*log10(fad);  %dB
average = mean(fad1);
figure
subplot(211);
plot(t,mt);%axis([-2,2,0,0.00002]);
xlabel('ʱ��');ylabel('ԭʼ�źŷ���')  
subplot(212);
plot(t,y);
xlabel('ʱ��');ylabel('˥�����ŵ����� ')  
figure
plot(t,fad1);
xlabel('ʱ��');ylabel('С�߶�˥����� dB')  
%-----��С�߶�˥��-----%


%----------------�Ӵ�߶�˥��--------------------
f = 1900e6
L = 3*10^8/f;             %����  0.1579m
D = 1;                    %���ߵ����ߴ� m
df = 2*D^2/L;             %����Զ���ο�����  12m
d0 = 100;                 %���ܹ��ʵĲο���
PL_d0 = 10*log10(L^2/(16*pi^2)*d0^2);  %����·�����
%����ʱ���ŵ���·������迼��������أ����ö�����̬��Ӱ
X = normrnd(0, 1);        %��˹�������
pn = 2.7;                   %·�����ָ��
dt1 = (1000-100)/N;
d = 100:dt1:1000-dt1;
d1 = log(d);              %logĬ����log2����
PL = PL_d0 - (10*pn*log(d/d0) + X);     %�ض�����d�µ�·�����
PSL = PL + fad1;
beishu = 10.^(PSL/20);
figure
plot(d1,PL,'b');hold on;
%axis([4.5,6,-80,-30]);
xlabel('���;���log��d��')
ylabel('��·����� dB')

plot(d1,PSL,'r');
%axis([4.5,6,-80,-30]);
xlabel('���;���log��d��')
ylabel('·����ļ���Ӱ˥��ͶྶЧӦ dB')
aa=max(PSL);
bb=min(PSL);
fprintf('�ľ�˥�����ޣ� %6.3f dB\n', bb);
fprintf('�ľ�˥�����ޣ� %6.3f dB\n', aa);
