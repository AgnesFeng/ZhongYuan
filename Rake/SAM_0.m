function [selq] = SAM_0(mt,sample,NUM,Ts)
%�������������� ԭ���У��������������������Ԫ�����壩��,��Ԫ����

N_sample = sample; % ������Ԫ��������
dt = Ts / N_sample; % ����ʱ����
N = NUM; % ��Ԫ��
t = 0 : dt : (N * N_sample - 1) * dt; % ���д���ʱ��
gt1 = ones(1, N_sample); % �����㲨��
% �����������
RAN = round(mt); % ���0 1����
selq = [];
for i = 1 : N % ��������
   if RAN(i)==1
       selq = [selq gt1]; 
   else
       selq = [selq 0*gt1];
   end
end

end