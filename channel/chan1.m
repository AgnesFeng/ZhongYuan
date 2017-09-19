% close all;
clear all;
sample = 1500e6;        
N = 3000;

dt = 1/sample; %ʱ�������� 
T = N*dt;
t = 0:dt:T-dt;                            
OR=sample;%�۲�Ƶ�ʣ�Hz��

P=0;
K=12;%ֱ�߱������˹�ֲ���K����   ��˹����
fd=25;%��������Ƶ�Ʋ���

%%%%%%%%%%%������˹�ֲ��г������ֺ�������ֵĹ���%%%%%%%%%%%%%%%%%%%%%%%
P = 10.^(P/10);%�������Զ����
s2 = P./(K+1); % ���㷽��
m2 = P.*(K./(K+1)); % ���㳣������
m = sqrt(m2); % ���㳣������

%%%%%%%%%���ض�����������˹�ŵ�%%%%%%%%%%%%%%%%
L = length(P); % ������ L=3
paths_r = sqrt(1/2)*(randn(L,N) + 1i*randn(L,N)).*((sqrt(s2))' * ones(1,N)); %L*N����ÿ�׵�������������K�й�
paths_c = m' * ones(1,N);%��������
doppler = ones(1,N);
doppler = ray_doppler(fd,dt,N);
paths_OR = paths_r.*doppler + paths_c(1,:); 
pp1 = abs(paths_OR);
pp2 = 10*log10(pp1.^2);
x = 1:N;
x1 = x/sample;
figure
plot(x1,pp2);axis([0,0.000002,-20,15])
xlabel('ʱ�䣨s���۲�Ƶ��1500MHz'),ylabel('С�߶�˥�����dB');

%�Ӵ�߶�˥��
f = 1900e6;
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
PSL = PL + pp2;
figure
plot(d1,PL,'b');hold on;
plot(d1,PSL,'r');hold on;
%axis([4.5,6,-80,-30]);
xlabel('���;���log��d��');
ylabel('��·����� dB');
