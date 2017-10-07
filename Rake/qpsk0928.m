close all;
clear all;clc;
%创建一个平方根升余弦滤波器
txfilter = comm.RaisedCosineTransmitFilter;
%创建一个相位频率偏移量对象，其中频率偏移量设置为-250赫兹，使用名值对将采样值设置为4000赫兹。
pfo = comm.PhaseFrequencyOffset(...
    'FrequencyOffset',-250,'SampleRate',4000);
%创建一个PSK粗糙频率估计系统对象，采样率为4千赫，频率分辨率为1赫兹。
frequencyEst = comm.PSKCoarseFrequencyEstimator(...
    'SampleRate',4000,'FrequencyResolution',1);
%创建第二个相位频率偏移量来校正偏移量。将频率设置的属性设置为输入端口，以便频率校正估计是输入参数。
pfoCorrect = comm.PhaseFrequencyOffset(...
    'FrequencyOffsetSource','Input port', ...
    'SampleRate',4000);
%生成QPSK信号，过滤信号，应用频率偏移，并通过AWGN信道传递信号。
modData = pskmod(randi([0 3],4096,1),4,pi/4);     % Generate QPSK signal
%txFiltData = txfilter(modData);                   % Apply Tx filter
txFiltData = step(txfilter,modData);
%offsetData = pfo(txFiltData);                     % Apply frequency offset
offsetData = step(pfo,txFiltData);
noisyData = awgn(offsetData,25);                  % Pass through AWGN channel
%用频率来估计频率偏移量。观察到这个估计值接近-250赫兹目标。
%estFreqOffset = frequencyEst(noisyData);
estFreqOffset = step(frequencyEst,noisyData);
%使用pfoCorrect和估计频率偏移量的倒数来校正频率偏移量。
%compensatedData = pfoCorrect(noisyData,-estFreqOffset);
compensatedData = step(pfoCorrect,noisyData ,-estFreqOffset);
%创建一个频谱分析仪对象来查看信号的频率响应。
spectrum = dsp.SpectrumAnalyzer('SampleRate',4000, 'ShowLegend',true, ...
    'ChannelNames',{'Received Signal' ,'Compensated Signal'});
%将接收到的信号的频率响应绘制到左边，并使用频谱分析仪对补偿信号进行偏移。补偿信号现在以正确的中心为中心。
spectrum([noisyData compensatedData]);