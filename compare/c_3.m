function p_sl = c_3(sample,N)
%sample = 1500e6;N = 100000;
dt = 1/sample; %ʱ�������� 
P=[0 -9 -12];%ÿһ�׵Ĺ���˥��
K=[6 0 0];%ֱ�߱������˹�ֲ���K����   ��˹����
tau=[0 1 2];%ÿ����ͷ��ʱ��(΢��)
delay = [0,0.000001,0.000002];
fd=25;%��������Ƶ�Ʋ���

%%%%%%%%%%%����ÿһ����˹�ֲ��г������ֺ�������ֵĹ���%%%%%%%%%%%%%%%%%%%%%%%
P = 10.^(P/10);%�������Զ����
s2 = P./(K+1); % ���㷽��
m2 = P.*(K./(K+1)); % ���㳣������
m = sqrt(m2); % ���㳣������
%%%%%%%%%%%������(Root Mean Square) ��ʱ%%%%%%%%%%%%%%%%
rmsdel = sqrt( sum(P.*(tau.^2))/sum(P) - (sum(P.*tau)/sum(P))^2 );
fprintf('������ʱ��rms �� %6.3f ��s\n', rmsdel);

%%%%%%%%%���ض�����������˹�ŵ�%%%%%%%%%%%%%%%%
L = length(P); % ������ L=3
paths_r = sqrt(1/2)*(randn(L,N) + 1i*randn(L,N)).*((sqrt(s2))' * ones(1,N)); %L*N����ÿ�׵�������������K�й�
paths_c = m' * ones(1,N);%��������
%λ��
t_shift=floor(delay/dt);%��һ��������ʱ  ����1��
%Ϊ�����������ź����ӳٴ���
rr1 = paths_r(1,:); 
rr2 = paths_r(2,:);
rr2(1,t_shift(2)+1:end) = rr2(1,1:end-t_shift(2));
rr3 = paths_r(3,:);
rr3(1,t_shift(3)+1:end) = rr3(1,1:end-t_shift(3));
%����
paths_sum = rr1 + rr2 + rr3;
doppler = ray_doppler(fd,dt,N); %25Hz,ȷʵ��0.4Hz��Ӱ���
paths_sum1 = paths_sum.*doppler; 
paths_OR = paths_sum1 + paths_c(1,:); %ע�⣬���Ƿֿ��ĵ����ģ���������ֻ�ӵ���һ������
pp1 = abs(paths_OR); %����
pp2 = 10*log10(pp1.^2); %dB

%�Ӵ�߶�˥��
f = 1900e6;
L = 3*10^8/f;             %����  0.1579m
D = 1;                    %���ߵ����ߴ� m
df = 2*D^2/L;             %����Զ���ο�����  12m
d0 = 100;                 %���ܹ��ʵĲο���
PL_d0 = 10*log10(L^2/(16*pi^2)*d0^2);  %����·�����
%����ʱ���ŵ���·������迼��������أ����ö�����̬��Ӱ
X = normrnd(0, 1);        %��˹�������
pn = 4;                   %·�����ָ��
d = input('�����뷢�;��룺'); 
PL = PL_d0 - (10*pn*log(d/d0) + X);     %�ض�����d�µ�·�����
PSL = PL + pp2; %dB
p_sl = 10.^(PSL./10);
end