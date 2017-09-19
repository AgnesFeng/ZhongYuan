%-----------------------1��5����-------------------------------------------
clc; close all; clear all;
Tlen = 500;              % message length
R = 1/3;              % code rate

g{1} = [1 0 1 1];                % Impulse Responses _ 1
g{2} = [1 1 0 1];                % Impulse Responses _ 2
g{3} = [1 1 1 1];                % Impluse Responses _ 3
%n = length(g);                       % Convolution Code (1/n) parameter,n=3
memory_els = 3;%K=4�Ĵ�������,m=K-1,Լ������


initial = round(rand(1,Tlen));            % message vector��Ϣλ���ȣ�L

%--------- ENCODER: 1/3 Convolution Encoder -------%
s_initial = encode_1_3(initial,g,length(g));%m����Ϣλ(1*L),g�����ɾ��������Ѿ�������n�����ʵĵ�����
for i = 1:length(s_initial)%��˫����
    if s_initial(i) == 0;
        s_initial(i) = -1;
    end
end
%-----------------------------------------
% for i = 1:length(s_initial)%�䵥����
%     if s_initial(i) == -1;
%         s_initial(i) = 0;
%     end
% end
r = s_initial;
[mhat,node] = decode_1_3(r,length(g),memory_els,Tlen,0);%r�ǽ��յ�������
Bit_Error_Number1 = length(find(mhat(1:Tlen) ~= initial(1:Tlen)));     %��ԭʼ�źűȽ�
%     Bit_Error_Rate1(snr) = Bit_Error_Number1/Tlen;