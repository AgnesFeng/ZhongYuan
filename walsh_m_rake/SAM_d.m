function [selq] = SAM_d(mt,sample,NUM)
%�������������� ԭ���У��������������������Ԫ�����壩��,��Ԫ����
N = NUM; % ��Ԫ��
gt1 = ones(1, sample); % �����㲨��
% �����������
RAN = round(mt); % ���0 1����
selq = [];
for i = 1 : N % ��������
   if RAN(i)==1
       selq = [selq gt1]; 
   else
       selq = [selq -1*gt1];
   end
end

end