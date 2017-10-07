clear all;        %������һ�����ڵ�m���л�ͼЧ������
close all;
clc;

%%--------------------------1��m����----------------
m = m_sequence([0,0,0,1,1,1,0,1]); %�õ�һ�����ڵ�m���� (���Ĵ����ĳ�ʼ״̬��������ΪC1~Cn��c0Ĭ��Ϊ1) 
N = length(m);%m��:225
v = 225000;
Tb = 1/v; %��Ԫ����ʱ��  
Tc = Tb/N;
Tm_sample = 40;%������Խ�󣬽��ͼԽ��׼��ע�⣬m�����ظ��󳤶ȿ��ܲ�����
m1 =  SAM_0(m,Tm_sample,length(m),Tc); %������ʱ��
m2 = 1-2*m1;  %��˫����
L = length(m2);
NP = 3;
mm = repmat(m,1,NP); %3�����ڵ�m���У����������
mm1 = SAM_0(mm,Tm_sample,length(mm),Tc);%�����ڣ��г���ʱ��
mm2 = 1-2*mm1;       %��˫����
dt = Tc/Tm_sample;
t=0:dt:Tc-dt;
n = 1:length(m2);
rt1=conv(mm2(1:3*L),m2(end:-1:1))/(N*Tm_sample);%ȡ�м��һ����
figure
index = L+Tm_sample;
figure(1);plot(rt1(index:end-index));title('������m���е�����غ���');


%%
%----------------------2��walsh����---------------------
wal2 = [ 1 1; 1 -1 ];
wal4 = [wal2 wal2; wal2 wal2*(-1)];  
wal8 = [wal4 wal4; wal4 wal4*(-1)];  %8*8
wal16 = [wal8 wal8; wal8 wal8*(-1)]; %16*16
p = wal16(5,:);
pp = SAM2(p,Tm_sample,length(p),Tc);%�����ڣ��г���ʱ��
%%
%---------------------3���ź����У���Ƶ,��ǰ������---------------------
Tlen = 6000; 
s_initial = randsrc( 1, Tlen );
T_sample = Tm_sample * length(p);  %m���еĲ�������Ϊ4
s_initial1 =  SAM2(s_initial,T_sample,length(s_initial),Tb);
pt = repmat(pp,1,Tlen);
spreadDate = pt.*s_initial1;%��Ƶ
spreadDate_m = ones(1,length(m2)+length(s_initial1));%��ǰ������
spreadDate_m(1:length(m2)) = m2;
spreadDate_m(length(m2)+1:end) = s_initial1;
len_spread = length(spreadDate_m);
%%
%�ҵ���λ����ʱ���±�
x1 = spreadDate_m(1:length(m2));
x2 = m2;
xcor = dsp.Crosscorrelator;
% delay = dsp.Delay(30);
% x2 = step(delay,m2');
y = step(xcor,x2',x1'); %computes cross-correlation of x1 and x2
figure(2), subplot(211);plot(y); title('Correlated output')


%%
%--------------------����-------------------
%spreadDate_m = rand(1,128);len_spread = length(spreadDate_m);
fs = 1/dt;
tx = ones(1,len_spread/2);%����ת���������
tx_re = tx; 
tx_im = tx;
k2 = 1;
for k1=1:Tm_sample*2:len_spread-2*Tm_sample+1 
    tx_re(k2:k2+Tm_sample-1) = spreadDate_m(k1:k1+Tm_sample-1);
    tx_im(k2:k2+Tm_sample-1) = spreadDate_m(k1+Tm_sample:k1+2*Tm_sample-1);
    k2 = k2+Tm_sample;
end

txx = tx_re + 1i*tx_im;
%--------------------���ŵ�-----------------
fd = 25;           %��������Ƶ�� 
k = 6;
tau = [0 0.000001 0.000002];
pdb = [0,-9,-12];
ts = dt;
chan = ricianchan(ts,fd,k,tau,pdb); %�ྶ��Ҫ��ÿһ����������˥�䡣ts�������źŵĲ���ʱ�䣨s����fd����������Ƶ�ƣ�pdb��ƽ��·�����棨dB��
y = filter(chan,txx); 

EbNo = 0:1:25;
for snr = 1:length(EbNo)
   y_snr = awgn(y,snr,'measured');

    %--------------------���-------------------
    re_2=real(y_snr);
    im_2=imag(y_snr);
    rx = ones(1,len_spread);
    k2 = 1; 
    for k1=1:Tm_sample*2:len_spread-2*Tm_sample+1 
        rx(k1:k1+Tm_sample-1) = re_2(k2:k2+Tm_sample-1);
        rx(k1+Tm_sample:k1+2*Tm_sample-1) = im_2(k2:k2+Tm_sample-1);
        k2 = k2+Tm_sample;
    end
   %�о�
   receive = sign(rx);
Bit_error = length(find(receive ~= spreadDate_m)); 
error_rate(snr) = Bit_error/len_spread; 
end

%%
%------------��ط岶��
x2 = rx(1:length(m2));
xcor = dsp.Crosscorrelator;
y = step(xcor,x2',x1'); %computes cross-correlation of x1 and x2
figure(2), subplot(212);plot(y); title('Correlated output')







