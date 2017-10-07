clear all;
close all;
g=19; % G=10011  n=4  m���г�15
state=8; % state=1000  ��λ�Ĵ�����ʼ״̬
L=1000;
% m���в���
N=15;
mq=mgen(g,state,L); %��Ϊ1000������Ϊ15��m����
mq1 = [1 1 1 0 0 1 0];
% �����������
ms=conv(1-2*mq,1-2*mq(15:-1:1))/N;%ǰ��˳��ı�����ν������һ����������һ��������������ǡ�һ�����ڵĵ��ţ����������ھ����
[ms1,b] = xcorr(1-2*mq,'biased');
figure(1)
subplot(222)
stem(ms(15:end));% ����״ͼ (ͷβ����15��)
%axis([0 63 -0.3 1.2]); title('m���е����������')

% m���й��ɵ��źţ��������壩
N_sample=8;
Tc=1;
dt=Tc/N_sample;
t=0:dt:Tc*L-dt;  %%%%%%5
gt=ones(1,N_sample);
mt=sigexpand(1-2*mq,N_sample);%������������ʱ�䣬mq�е�1��-1��
mt=conv(mt,gt); %���γ��� ��˫���Եģ�ȫ��
figure(1)
subplot(221);
plot(t,mt(1:length(t)));
axis([0 63 -1.2 1.2]); title('m���о��γ����ź�');
% 
st=sigexpand(1-2*mq(1:15),N_sample); %˫����m���У�һ������
s=conv(st,gt);     %������ʱ��
st=s(1:length(st));%����󳤶ȱ��ˣ�ȡ��һ������
rt1=conv(mt,st(end:-1:1))/(N*N_sample);%���һ�е���һ�л�������,�����
subplot(223)
plot(t,rt1(length(st):length(st)+length(t)-1));
%plot(t,rt1(1:length(t)));
%plot(rt1);
axis([0 63 -0.2 1.2]); title('m���о��γ����źŵ������');xlabel('t');















% sinc����
% Tc=1;
% dt=Tc/N_sample;
% t=-20:dt:20;
% gt=sinc(t/Tc);
% mt=sigexpand(1-2*mq,N_sample);
% mt=conv(mt,gt);
% st2=sigexpand(1-2*mq(1:15),N_sample);
% s2=conv(st2,gt);
% st2=s2;
% rt2=conv(mt,st2(end:-1:1))/(N*N_sample);
% subplot(224)
% t1=-55+dt:dt:Tc*L-dt;
% plot(t1,rt2(1:length(t1)));
% axis([0 63 -0.5 1.2]); title('m����sinc�����źŵ������');xlabel('t');