
clear all;
close all;clc;
dt = 0.0000001;
Tm_sample = 4;
spreadDate_m = randsrc(1,128);len_spread = length(spreadDate_m);
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

