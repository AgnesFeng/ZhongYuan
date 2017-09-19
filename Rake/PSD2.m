function [sing]=PSD2(N_FFT,N_sample1,selq,NUM,Ts)
%功率谱密度计算方法二,叠加取平均.fft点数，原信号采样率，原信号，叠加脉冲个数

NN=N_FFT;%做NN点的傅立叶变换
N_sample = N_sample1;
N = NUM;
sumsel = zeros(1,NN);
for k=1:N_sample:length(selq)-N_sample
    fftsel = fftshift(fft(selq(k:k+N_sample-1),NN));
    a=10 * log10(abs(fftsel).^ 2 / (N_sample));   %原理这里是N_sample，NUM*Ts
    sumsel = sumsel + a;%计算功率
end
sing = (sumsel/N);

x=linspace(-12*pi*NUM/Ts,12*pi*NUM/Ts,NN);
plot(x,sing); grid on; %title('双极性不归零信号功率谱');
end