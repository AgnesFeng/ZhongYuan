clear all;
close all;
clc;
spread = randsrc(1,20000);
Tm_sample = 4;
N = length(spread);
len_spread = length(spread);
tx = ones(1,len_spread/2);%串并转换后的数组
tx_re = tx; 
tx_im = tx;
k2 = 1;
for k1=1:2:len_spread-1 %两两一组 ， 得到tx_re,im_re,长度16*8000 = 64000
    tx_re(k2) = spread(k1);
    tx_im(k2) = spread(k1+1);
    k2 = k2+1;
end
I = tx_re;
Q = tx_im;
txx = I + 1i*Q;
rx = ones(1,N);

%调制：4进制
tx = ones(1,N/2);%串并转换后的数组
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
    
    %解调
    re =real(out);
    im = imag(out);
    k4 = 1;
    for k1=1:2:len_spread-1 %实部虚部放到一个数组receive中  
       rx(k1) = re(k4);
       rx(k1+1) = im(k4);
       k4 = k4+1;
    end
    
    %解调4进制
    re_2=real(out_4);
    im_2=imag(out_4);
    rx_4 = ones(1,N);
    k2 = 1; 
    for k1=1:Tm_sample*2:len_spread-2*Tm_sample+1 
        rx_4(k1:k1+Tm_sample-1) = re_2(k2:k2+Tm_sample-1);
        rx_4(k1+Tm_sample:k1+2*Tm_sample-1) = im_2(k2:k2+Tm_sample-1);
        k2 = k2+Tm_sample;
    end     
     
    %计算误码率 no qpsk
    receive_no = sign(out_noqpsk);
    Bit_error_no = length(find(receive_no ~= spread)); 
    error_rate_no(x) = Bit_error_no/N;  
    %计算误码率
    receive = sign(rx);
    Bit_error = length(find(receive ~= spread)); 
    error_rate(x) = Bit_error/N; 
    %计算误码率4进制
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




