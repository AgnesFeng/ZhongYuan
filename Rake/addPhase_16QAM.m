clear all;
close all;
%对16-qam信号进行相位偏移，并观察其对星座的影响。
%创建一个QAM调制器和一个相位频率偏移系统对象。将相位偏移设为30度。
qamModulator = comm.RectangularQAMModulator('ModulationOrder',16);
pfo = comm.PhaseFrequencyOffset('PhaseOffset',30);
%生成随机符号，并应用16-qam调制。
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