clear all;        %������һ�����ڵ�m���л�ͼЧ������
close all;
clc;
primpoly(7,'all')%�õ�����7��(nΪ7)�����б�ԭ����ʽ
m = m_sequence([0,0,0,1,0,0,1]); %�õ�һ�����ڵ�m���� (���Ĵ����ĳ�ʼ״̬��������ΪC1~Cn��c0Ĭ��Ϊ1) 
N = length(m);%m��:127
v = 225000;
Tb = 1/v; %��Ԫ����ʱ��  
Tc = Tb/N;
Tm_sample = 4;%������Խ�󣬽��ͼԽ��׼��ע�⣬m�����ظ��󳤶ȿ��ܲ�����
T_sample = Tm_sample * N;  %m���еĲ�������Ϊ6
m1 =  SAM_0(m,Tm_sample,length(m),Tc); %������ʱ��
m2 = 1-2*m1;  %��˫����
L = length(m2);

NP = 300;
mm = repmat(m,1,NP); %������
mm1 = SAM_0(mm,Tm_sample,length(mm),Tc);%�����ڣ��г���ʱ��
mm2 = 1-2*mm1;       %��˫����
dt = Tc/Tm_sample;
t=0:dt:Tc-dt;
rt1=conv(mm2(1:3*L),m2(end:-1:1))/(N*Tm_sample);%ȡ�м��һ����
figure
index = L+Tm_sample;
plot(rt1);title('������m���е�����غ���(�г���ʱ��)');
% x1 = m2; %������ʵ��
% x2 = [1 1 1];
% [a,b] = xcorr(x1,x2);
%---------------------------------ԭʼ�ź�---------------------------------
Tlen = 100; %���ݳ���  �������������9�Ĺ�ϵ
s_initial = rands( 1, Tlen ); 
figure,plot(s_initial);
s_initial1 =  SAM2(s_initial,T_sample,length(s_initial),Tb);
mseq = mm2(1:length(s_initial1));

spreadData = mseq.*s_initial1;
first = spreadData(1:L);%һ�����ڵ�m���еĳ��ȣ���һ����Ԫ
test = conv(mm2(1:3*L),first(end:-1:1))/L;%��λ�����ֵ��仯����Ҫ�趨һ����ֵ������ʱ����

figure
plot(test);
%---------------------------------λ��---------------------------------
rr2 = spreadData;   %1us��1/4��Ԫ
shift = 3;
rr2(1,round(T_sample/shift+1):end) = rr2(1,1:end-round(T_sample/shift));    %rr2(1,1:T_sample/4)=rr2(1,end-T_sample/4+1:end);
%˵������Ԫ������ʱ�з�ֵ������Ϊ1�����Ը������ֵ��λ��ȷ��λ�ƣ�Ҳ���Ը��ݻ��������ľ���ȷ��λ�ƣ�������Ҫ��ǰ�趨��ֵ
%�����û�������ȷ�����������ֵ������
RR2 = rr2(end-L:end);
test1 = conv(mm2(1:3*L),RR2(end:-1:1))/L;
figure
plot(test1);
%�۲�m������λ������ػ᲻�ᵽ���ֵ
mafter = mm2(1:3*L);
mafter(1,round(T_sample/shift)+1:end) = mafter(1,1:end-round(T_sample/shift));   
mafter(1,1:round(T_sample/shift)) = mafter(1,end-round(T_sample/shift)+1:end);
test2 = conv(mafter,RR2(end:-1:1))/L;
figure
plot(test2);
%---------------------------------����---------------------------------

