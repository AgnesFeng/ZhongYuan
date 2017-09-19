%-----------------------1、5编码-------------------------------------------
clc; close all; clear all;
Tlen = 500;              % message length
R = 1/3;              % code rate

g{1} = [1 0 1 1];                % Impulse Responses _ 1
g{2} = [1 1 0 1];                % Impulse Responses _ 2
g{3} = [1 1 1 1];                % Impluse Responses _ 3
%n = length(g);                       % Convolution Code (1/n) parameter,n=3
memory_els = 3;%K=4寄存器个数,m=K-1,约束长度


initial = round(rand(1,Tlen));            % message vector信息位长度：L

%--------- ENCODER: 1/3 Convolution Encoder -------%
s_initial = encode_1_3(initial,g,length(g));%m是信息位(1*L),g是生成矩阵，上面已经给出，n是码率的倒数。
for i = 1:length(s_initial)%变双极性
    if s_initial(i) == 0;
        s_initial(i) = -1;
    end
end
%-----------------------------------------
% for i = 1:length(s_initial)%变单极性
%     if s_initial(i) == -1;
%         s_initial(i) = 0;
%     end
% end
r = s_initial;
[mhat,node] = decode_1_3(r,length(g),memory_els,Tlen,0);%r是接收到的序列
Bit_Error_Number1 = length(find(mhat(1:Tlen) ~= initial(1:Tlen)));     %与原始信号比较
%     Bit_Error_Rate1(snr) = Bit_Error_Number1/Tlen;