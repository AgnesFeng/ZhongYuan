clear all;
close all;
%��16-qam�źŽ�����λƫ�ƣ����۲����������Ӱ�졣
%����һ��QAM��������һ����λƵ��ƫ��ϵͳ���󡣽���λƫ����Ϊ30�ȡ�
qamModulator = comm.RectangularQAMModulator('ModulationOrder',16);
pfo = comm.PhaseFrequencyOffset('PhaseOffset',30);
%����������ţ���Ӧ��16-qam���ơ�
data = (0:15)';
modData = step(qamModulator,data);
scatterplot(modData);
title(' Original Constellation')
xlim([-5 5])
ylim([-5 5])
%
impairedData = step(pfo,modData);
scatterplot(impairedData);
title('Constellation after phase offset')
xlim([-5 5])
ylim([-5 5])