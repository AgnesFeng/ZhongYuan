function [c] = REPEAT(mt,sample)
%�������������� ԭ���У��������������������Ԫ�����壩��,��Ԫ����
% clear all;
% close all;
% mt = [1,3,4,2,1,3,4];
% sample = 5;
N = length(mt); % ��Ԫ��
a = ones(1,sample);
b = ones(N,sample);
for i = 1:N
     b(i,:)= mt(i)*a;
end
c = reshape(b',1,N*sample);