function [Waveform] = SymbolToWaveform(SymbolStream,numSamplesPerBit) 
%SymbolStream��ԭ�����ź�   numSamplesPerBit����������Ĳ�������
lenWaveform = length(SymbolStream)*numSamplesPerBit; 
Waveform = zeros(1,lenWaveform); 
Waveform(1:numSamplesPerBit:lenWaveform) = SymbolStream; 
end