function [sing]=PSD2(N_FFT,N_sample1,selq,NUM,Ts)
%�������ܶȼ��㷽����,����ȡƽ��.fft������ԭ�źŲ����ʣ�ԭ�źţ������������

NN=N_FFT;%��NN��ĸ���Ҷ�任
N_sample = N_sample1;
N = NUM;
sumsel = zeros(1,NN);
for k=1:N_sample:length(selq)-N_sample
    fftsel = fftshift(fft(selq(k:k+N_sample-1),NN));
    a=10 * log10(abs(fftsel).^ 2 / (N_sample));   %ԭ��������N_sample��NUM*Ts
    sumsel = sumsel + a;%���㹦��
end
sing = (sumsel/N);

x=linspace(-12*pi*NUM/Ts,12*pi*NUM/Ts,NN);
plot(x,sing); grid on; %title('˫���Բ������źŹ�����');
end