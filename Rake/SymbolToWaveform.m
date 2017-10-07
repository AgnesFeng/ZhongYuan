function [Waveform] = SymbolToWaveform(SymbolStream,numSamplesPerBit) 
%SymbolStream：原数字信号   numSamplesPerBit：单个脉冲的采样点数
lenWaveform = length(SymbolStream)*numSamplesPerBit; 
Waveform = zeros(1,lenWaveform); 
Waveform(1:numSamplesPerBit:lenWaveform) = SymbolStream; 
end