clear all;        %������һ�����ڵ�m���л�ͼЧ������
close all;
clc;
primpoly(7,'all')%�õ�����7��(nΪ7)�����б�ԭ����ʽ
m = m_sequence([0,0,0,1,1,1,0,1]); %�õ�һ�����ڵ�m���� (���Ĵ����ĳ�ʼ״̬��������ΪC1~Cn��c0Ĭ��Ϊ1) 
N = length(m);
Tc = 1;
Tm_sample = 6;%������Խ�󣬽��ͼԽ��׼��ע�⣬m�����ظ��󳤶ȿ��ܲ�����
T_sample = 6 * N;  %m���еĲ�������Ϊ6
m1 =  SAM_0(m,Tm_sample,length(m),Tc); %������ʱ��
m2 = 1-2*m1;  %��˫����
L = length(m2);

NP = 3;
mm = repmat(m,1,NP); %������
mm1 = SAM_0(mm,Tm_sample,length(mm),Tc);%�����ڣ��г���ʱ��
mm2 = 1-2*mm1;       %��˫����
dt = Tc/Tm_sample;
t=0:dt:Tc-dt;
rt1=conv(mm2,m2(end:-1:1))/(N*Tm_sample);%ȡ�м��һ����
figure(2)
index = L+Tm_sample;
plot(rt1(index:end-index));title('������m���е�����غ���(�г���ʱ��)');
% x1 = m2; %������ʵ��
% x2 = [1 1 1];
% [a,b] = xcorr(x1,x2);
%--------------------------3��������ֵ------(Ӧ�÷�����Ƶ�������)-----------------

