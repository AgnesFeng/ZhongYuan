clear all;        %������һ�����ڵ�m���л�ͼЧ������
close all;
primpoly(3,'all');%�õ�����7��(nΪ7)�����б�ԭ����ʽ
m = m_sequence([0,1,0,0,1]); %�õ�һ�����ڵ�m���� (���Ĵ����ĳ�ʼ״̬��������ΪC1~Cn��c0Ĭ��Ϊ1) 
%m = [1 0 0 1 0 0 1 1 0 1 0 1 1 1]; %��ֵ
%m = [1 1 1 0 0 1 0]; %��m����
N = length(m);
v = 225;
Tb = 1/v; %��Ԫ����ʱ��
Tc = Tb/N;
T_sample = 310;%������Խ�󣬽��ͼԽ��׼��ע�⣬m�����ظ��󳤶ȿ��ܲ�����
Tm_sample = T_sample/N;  %m���еĲ�������Ϊ10
m1 =  SAM_0(m,Tm_sample,length(m),Tc); %������ʱ��
m2 = 1-2*m1;  %��˫����

%--------------------- 1����ԭʼ������غ͹����ܶ���-------------------------
figure(1)
[a1,f] = xcorr(m1,'coeff');%����coeff��Ƶ�����м䣬unbaisedҪ����
subplot(311);plot(f,(a1));title('m���е�����غ�����һ�����ڳ��ȣ��޳���ʱ��');%
%����Fourier�任�еľ�������������  
%Rt = fftshift(ifft(fft(m).*fft(m)));%m��m1��������һ���ġ�subplot(212);plot(Rt);title('��fft�������غ���');%�Ǵ��

selq2 = m1;
fft_se12 = fftshift(fft(selq2));
PE12 = 10 * log10(abs(fft_se12) .^ 2 / (N * Tc)); % ��ʽ����������ܶ�
PEL12 = (-length(fft_se12) / 2 : length(fft_se12) / 2 - 1) / 10; % �����䳤��
%%���Ƴ����
subplot(312);plot(PEL12, PE12); grid on; title('������m���й����ܶ���');
axis([-50 50 -50 50]);xlabel('f');ylabel('����');

%-------------------------2����m�����ظ�NP������----------------------------
NP = 500;
mm = repmat(m,1,NP); %������
mm1 = SAM_0(mm,Tm_sample,length(mm),Tc);%�����ڣ��г���ʱ��
mm2 = 1-2*mm1;       %��˫����
    %---------��Ƶ��
    ym = fftshift(fft(mm2,4096));
    magm = abs(ym);
    fm = (1:4096)*200/2048;
    subplot(313);plot(fm,magm*2/4096);title('������ڵ�31λm���е�Ƶ��');%m���е�Ƶ��
    %---------
L = length(mm);   %NP�����ڵ�m���ܳ�
dt = Tc/Tm_sample;
t=0:dt:L*Tc-dt;
rt1=conv(m2,m2(end:-1:1))/(N*Tm_sample);%���һ�е���һ�л�������,���,ֻ����һ��һ��
max = max(rt1);
figure(2)
%plot(t,rt1(length(m2):length(m2)+length(t)-1));title('m���е�����غ���(һ�����ڣ��г���ʱ��)');
plot(rt1);
%axis([0 0.02 -0.2 1.2]);
%--------------------------[0916�����������-----------------
%һ�����ڵ�m������m2
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

%--------------------------0916�����������]-----------------

%--------------------------3��������ֵ------(Ӧ�÷�����Ƶ�������)-----------------
%��˫���ԡ��г���ʱ��ġ�һ�����ڵ�m����mm2�ӳ�Tc/2
% m2(1,Tm_sample/2+1:end) = m2(1,1:end-Tm_sample/2);
% m2(1,1:Tm_sample/2)=m2(1,end-Tm_sample/2+1:end); %m2���ˣ����Ը������Ե�m1���Ƿ�λ��
% limid=conv(mm2,m2(end:-1:1))/(N*Tm_sample);%���һ�е���һ�л�������
% figure(3),plot(t,rt1(length(m2):length(m2)+length(t)-1));title('m������λ��������');
% axis([0 200 -0.2 1.2]);
   
%-----------------------5��ԭʼ�ź�-------------------------------------------
% Tlen = 50; %���ݳ��� 
% s_initial = randsrc( 1, Tlen );
% %Tb = 31*Tc; %������m���еĲ�����û�й�ϵ����Ƶʱm������û�в�����
% s_initial1 =  SAM2(s_initial,T_sample,Tlen,Tb); 
% %----��ʱ������
% dt = Tb/T_sample;   %dt��һ�µ�
% t = 0 : dt : (Tlen * T_sample - 1) * dt; % ���д���ʱ��
% figure(4);subplot(211);plot(t,s_initial1);axis([0,0.1,-2,2]);title('��Ƶǰʱ��');
% %----��Ƶ������
% selq2 = s_initial1;
% fft_se12 = fftshift(fft(selq2)); % �����е�Ƶ��
% PE12 = 10 * log10(abs(fft_se12) .^ 2 / (Tlen * Tb)); % ��ʽ����������ܶ�
% PEL12 = (-length(fft_se12) / 2 : length(fft_se12) / 2 - 1) / 10; % �����䳤��
% figure(5);subplot(2, 1, 1);plot(PEL12, PE12); grid on; title('��Ƶǰ�����ܶ���');
% axis([-200 200 -50 100]);xlabel('f');ylabel('����');
% 
% %---------------------------6����Ƶ(�ĳ�����λ���)------------------------------------------   
% mseq = mm2(1:length(s_initial1));
% spreadData = s_initial1.*mseq;       %��Ƶ����ź�
% %----��ʱ������
% figure(4);subplot(212);plot(t,spreadData);axis([0,0.1,-2,2]);title('��Ƶ��ʱ��');
% %----��Ƶ������
% selq2 = spreadData;
% fft_se12 = fftshift(fft(selq2)); % �����е�Ƶ��
% PE12 = 10 * log10(abs(fft_se12) .^ 2 / (Tlen * Tb)); % ��ʽ����������ܶ�
% PEL12 = (-length(fft_se12) / 2 : length(fft_se12) / 2 - 1) / 10; % �����䳤��
% figure(5);subplot(2, 1, 2);plot(PEL12, PE12); grid on; title('��Ƶ�����ܶ���');
% axis([-200 200 -500 100]);xlabel('f');ylabel('����');
% 
% figure(6)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% spreadshifft1 = spreadData;    %����1/2��Ƭ    �������ı�������ֵ���ĳ�1/2��Ԫ
% spreadshifft1(1,T_sample/2+1:end) = spreadshifft1(1,1:end-T_sample/2);
% spreadshifft1(1,1:T_sample/2)=spreadshifft1(1,end-T_sample/2+1:end); %m2���ˣ����Ը������Ե�m1���Ƿ�λ��
% [a1,f] = xcorr(spreadshifft1,spreadData,'coeff');%����coeff��Ƶ�����м䣬unbaisedҪ����
% subplot(312);plot(f,(a1));title('����1/2��Ԫʱ����غ���');axis([-500,500,-0.2,1.2])
% 
% spreadshifft2 = spreadData;    %����1��Ԫ   
% spreadshifft2(1,T_sample+1:end) = spreadshifft2(1,1:end-T_sample);
% spreadshifft2(1,1:T_sample)=spreadshifft2(1,end-T_sample+1:end); %m2���ˣ����Ը������Ե�m1���Ƿ�λ��
% [a2,f] = xcorr(spreadshifft2,spreadData,'coeff');%����coeff��Ƶ�����м䣬unbaisedҪ����
% subplot(313),plot(f,(a2));title('����1��Ԫʱ��غ���');axis([-500,500,-0.2,1.2])
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [a,f] = xcorr(spreadData,'coeff');%����coeff��Ƶ�����м䣬unbaisedҪ����
% subplot(311);plot(f,(a));title('��λ����ʱ������غ���');
% axis([-500,500,-0.2,1.2]);
% %---------------------------6~7��Ƶ�װ��Ƶ���----------------
% 
% 
% %---------------------------7�����ŵ�----------------
% chanData = awgn(spreadData,15,'measured');
% %1
% %2
% %3
% %---------------------------8�����------------------(�����ڽ��ǰ��������޳���ʱ���)----------------------
% receive = chanData;
% 
% %---------------------------9����ط岶��-------------
% %for tt = 0:dt:T_samople
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%�ѵõ���λ��õ�����ʱ��
% %---------------------------10������------------------
% %������     ��m���������û����ط壬Ҫ�ú�m������˺������
% demseq = mseq;  
% demseq(1,T_sample+1:end) = demseq(1,1:end-T_sample);
% demseq(1,1:T_sample)=demseq(1,end-T_sample+1:end);%demseq�Ƕ�����λ���
% demp = demseq.* spreadshifft2;
% 
% temp = s_initial1; %���Ժϲ�ǰҪ��λ����
% temp(1,T_sample+1:end) = temp(1,1:end-T_sample);
% temp(1,1:T_sample)=temp(1,end-T_sample+1:end);%demseq�Ƕ�����λ���
% a = demp-temp;
% [a1,f] = xcorr(demp,s_initial1,'coeff');%����coeff��Ƶ�����м䣬unbaisedҪ����
% figure(8);plot(f,(a1));title('��һ����������Լ���'); 
% 
% demseq2 = mseq;
% demp = mseq.*spreadData;
% b = demp-s_initial1;
% %��һ�����ӵ��ƽ����


