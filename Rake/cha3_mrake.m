function [P_sl] = cha3_mrake(sample,n)
N=n;%���������ʵ�ֵ���Ŀ�� ��Ϊ������Ƿ��ȣ���������������򷽱�������

OR=sample;%���滭ͼ������ļ����1/OR
M=256;%�������˲����Ľ���

P=[0 -9 -12];%ÿһ�׵Ĺ���˥��
K=[6 0 0];%ֱ�߱������˹�ֲ���K����   ��˹����
tau=[0 1 2];%ÿ����ͷ��ʱ��(΢��)
Dop=[25 25 25];%��������Ƶ�Ʋ���

%%%%%%%%%%%����ÿһ����˹�ֲ��г������ֺ�������ֵĹ���%%%%%%%%%%%%%%%%%%%%%%%
P = 10.^(P/10);%�������Զ����
s2 = P./(K+1); % ���㷽��
m2 = P.*(K./(K+1)); % ���㳣������
m = sqrt(m2); % ���㳣������
%%%%%%%%%%%������(Root Mean Square) ��ʱ%%%%%%%%%%%%%%%%
rmsdel = sqrt( sum(P.*(tau.^2))/sum(P) - (sum(P.*tau)/sum(P))^2 );
fprintf('rms delay spread %6.3f ��s\n', rmsdel);
%%%%%%%%%���ض�����������˹�ŵ�%%%%%%%%%%%%%%%%
L = length(P); % ������ L=3
paths_r = sqrt(1/2)*(randn(L,N) + 1i*randn(L,N)).*((sqrt(s2))' * ones(1,N)); %L*N����ÿ�׵�������������K�й�
paths_c = m' * ones(1,N);%��������

for p = 1:L %�����е�����
    D = Dop(p) / max(Dop) / 2; % ��һ����������Ƶ�� 
    f0 = (0:M*D)/(M*D); % Ƶ������ 
    PSD = 0.785*f0.^4 - 1.72*f0.^2 + 1.0; % PSD���Ʒ�
    filt = [ PSD(1:end-1) zeros(1,floor(M-2*M*D)) PSD(end:-1:2) ]; % S(f) 
    
    filt = sqrt(filt); %��S(f)��|H(f)| 
    filt = ifftshift(ifft(filt)); % ���������Ӧ 
    filt = real(filt); % Ѱ��ʵ���˲���
    filt = filt / sqrt(sum(filt.^2)); %��һ���˲���
    path = fftfilt(filt, [ paths_r(p,:) zeros(1,M) ]); 
    paths_r(p,:) = path(1+M/2:end-M/2); %�������
end; 
paths = paths_r + paths_c;  %����˥������(���SUI�����������һ������)    ��һ�����ӵ����壿������������������
Pest = mean(abs(paths).^2, 2);%ƽ����ͷ�ܹ���mean��A��2���������A���еľ�ֵ��mean��A,1��=mean(A)
fprintf('tap mean power level: %0.2f dB\n', 10*log10(Pest));%��ӡÿ�е�ƽ��ֵ

paths_OR = paths;
paths_OR1=10*log10(paths_OR); %����۲�������
%NN=length(paths_OR(1,:)); %���������
y1=abs(paths_OR).^2; %����   ����     ֱ��ȡģ�õ����ȣ�ƽ�����ǹ���
y2=10*log10(y1); %ת��ΪdB
%t=60; %ʱ�䳤��
x=1:n; %(����n�Ͳ���)
x1=x/OR;

%%%%%%%%%%%%%%����С�߶�˥�䱶��%%%%%%%%%%%%%
pp = paths_OR(1,:) + paths_OR(2,:) + paths_OR(3,:); %���ն�ͬһʱ�̵���ķ��ȣ���˲�����ʱֱ�ӵ��� 
pp1 = abs(pp);  %����
pp2 = 10*log10(pp1.^2); %����dB

%----------------�Ӵ�߶�˥��--------------------
f = 1900e6;
L = 3*10^8/f;             %����  0.1579m
D = 1;                    %���ߵ����ߴ� m
df = 2*D^2/L;             %����Զ���ο�����  12m
d0 = 100;                 %���ܹ��ʵĲο���
PL_d0 = 10*log10(L^2/(16*pi^2)*d0^2);  %����·�����
%����ʱ���ŵ���·������迼��������أ����ö�����̬��Ӱ
X = normrnd(0, 1);        %��˹�������
pn = 2.8;                   %·�����ָ��
d = input('�����뷢�;��룺'); 
d1 = log(d);              %logĬ����log2����
PL = PL_d0 - (10*pn*log(d/d0) + X);     %�ض�����d�µ�·�����
%PSL = PL + pp2(1:901);
pl = 10^(PL./10);         %��log���ɱ���
P_sl = pl*ones(3,n);
P_sl(1,:) = abs(paths_OR(1,:))+ P_sl(1,:);
P_sl(2,:) = abs(paths_OR(2,:))+ P_sl(2,:);
P_sl(3,:) = abs(paths_OR(3,:))+ P_sl(3,:);
end