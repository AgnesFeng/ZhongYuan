close all;
clear all;clc;
%����һ��ƽ�����������˲���
txfilter = comm.RaisedCosineTransmitFilter;
%����һ����λƵ��ƫ������������Ƶ��ƫ��������Ϊ-250���ȣ�ʹ����ֵ�Խ�����ֵ����Ϊ4000���ȡ�
pfo = comm.PhaseFrequencyOffset(...
    'FrequencyOffset',-250,'SampleRate',4000);
%����һ��PSK�ֲ�Ƶ�ʹ���ϵͳ���󣬲�����Ϊ4ǧ�գ�Ƶ�ʷֱ���Ϊ1���ȡ�
frequencyEst = comm.PSKCoarseFrequencyEstimator(...
    'SampleRate',4000,'FrequencyResolution',1);
%�����ڶ�����λƵ��ƫ������У��ƫ��������Ƶ�����õ���������Ϊ����˿ڣ��Ա�Ƶ��У�����������������
pfoCorrect = comm.PhaseFrequencyOffset(...
    'FrequencyOffsetSource','Input port', ...
    'SampleRate',4000);
%����QPSK�źţ������źţ�Ӧ��Ƶ��ƫ�ƣ���ͨ��AWGN�ŵ������źš�
modData = pskmod(randi([0 3],4096,1),4,pi/4);     % Generate QPSK signal
%txFiltData = txfilter(modData);                   % Apply Tx filter
txFiltData = step(txfilter,modData);
%offsetData = pfo(txFiltData);                     % Apply frequency offset
offsetData = step(pfo,txFiltData);
noisyData = awgn(offsetData,25);                  % Pass through AWGN channel
%��Ƶ��������Ƶ��ƫ�������۲쵽�������ֵ�ӽ�-250����Ŀ�ꡣ
%estFreqOffset = frequencyEst(noisyData);
estFreqOffset = step(frequencyEst,noisyData);
%ʹ��pfoCorrect�͹���Ƶ��ƫ�����ĵ�����У��Ƶ��ƫ������
%compensatedData = pfoCorrect(noisyData,-estFreqOffset);
compensatedData = step(pfoCorrect,noisyData ,-estFreqOffset);
%����һ��Ƶ�׷����Ƕ������鿴�źŵ�Ƶ����Ӧ��
spectrum = dsp.SpectrumAnalyzer('SampleRate',4000, 'ShowLegend',true, ...
    'ChannelNames',{'Received Signal' ,'Compensated Signal'});
%�����յ����źŵ�Ƶ����Ӧ���Ƶ���ߣ���ʹ��Ƶ�׷����ǶԲ����źŽ���ƫ�ơ������ź���������ȷ������Ϊ���ġ�
spectrum([noisyData compensatedData]);