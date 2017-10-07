function [c] = REPEAT(mt,sample)
%方波采样函数。 原序列，单个脉冲采样点数，码元（脉冲）数,码元周期
% clear all;
% close all;
% mt = [1,3,4,2,1,3,4];
% sample = 5;
N = length(mt); % 码元数
a = ones(1,sample);
b = ones(N,sample);
for i = 1:N
     b(i,:)= mt(i)*a;
end
c = reshape(b',1,N*sample);