clear all;
close all;
clc;
spread = randsrc(1,20000);
Tm_sample = 4;
N = length(spread);
len_spread = length(spread);
tx = ones(1,len_spread/2);%����ת���������
tx_re = tx; 
tx_im = tx;
k2 = 1;
for k1=1:2:len_spread-1 %����һ�� �� �õ�tx_re,im_re,����16*8000 = 64000
    tx_re(k2) = spread(k1);
    tx_im(k2) = spread(k1+1);
    k2 = k2+1;
end
I = tx_re;
Q = tx_im;
txx = I + 1i*Q;
rx = ones(1,N);

%���ƣ�4����
tx = ones(1,N/2);%����ת���������
tx_re = tx; 
tx_im = tx;
k2 = 1;
for k1=1:Tm_sample*2:len_spread-2*Tm_sample+1 
    tx_re(k2:k2+Tm_sample-1) = spread(k1:k1+Tm_sample-1);
    tx_im(k2:k2+Tm_sample-1) = spread(k1+Tm_sample:k1+2*Tm_sample-1);
    k2 = k2+Tm_sample;
end
txx_4 = tx_re + 1i*tx_im;

x = 1;
for snr  = 1:0.2:10
    out_noqpsk = awgn(spread,8,'measured');
    out = awgn(txx,snr,'measured');
    out_4 = awgn(txx_4,snr,'measured');
    
    %���
    re =real(out);
    im = imag(out);
    k4 = 1;
    for k1=1:2:len_spread-1 %ʵ���鲿�ŵ�һ������receive��  
       rx(k1) = re(k4);
       rx(k1+1) = im(k4);
       k4 = k4+1;
    end
    
    %���4����
    re_2=real(out_4);
    im_2=imag(out_4);
    rx_4 = ones(1,N);
    k2 = 1; 
    for k1=1:Tm_sample*2:len_spread-2*Tm_sample+1 
        rx_4(k1:k1+Tm_sample-1) = re_2(k2:k2+Tm_sample-1);
        rx_4(k1+Tm_sample:k1+2*Tm_sample-1) = im_2(k2:k2+Tm_sample-1);
        k2 = k2+Tm_sample;
    end     
     
    %���������� no qpsk
    receive_no = sign(out_noqpsk);
    Bit_error_no = length(find(receive_no ~= spread)); 
    error_rate_no(x) = Bit_error_no/N;  
    %����������
    receive = sign(rx);
    Bit_error = length(find(receive ~= spread)); 
    error_rate(x) = Bit_error/N; 
    %����������4����
    receive_4 = sign(rx_4);
    Bit_error_4 = length(find(receive_4 ~= spread)); 
    error_rate_4(x) = Bit_error_4/N; 
    x = x+1;
end
figure
a = 1:x-1;
semilogy(a,error_rate_no,'b');hold on;
semilogy(a,error_rate,'r');hold on;
semilogy(a,error_rate_4,'m');hold on;




