function  [receive] = demod2(psf,cos1,sin1,len_spread,data,tap) %����ԭ����,��ͨ���󳤶�,ʱ���ᣬ��������
supersam = 5;
NT=50;
N=2*NT*supersam;    %�˲�������-1����length(psf)-1һ��
mtf=psf;
%���ز�
Idem=tap.*(sqrt(2)*cos1);
Qdem=tap.*(sqrt(2)*sin1);

%���˲�
Imat=conv(Idem,mtf);
Qmat=conv(Qdem,mtf);   
 
   for  i=1:supersam*data                 %data selection
       Isel(i)=Imat(i+N);
       Qsel(i)=Qmat(i+N);
   end
              
   for i=1:data                       %��ȡ��Ԫ   
       Isam(i)=Isel((i-1)*supersam+1);
       Qsam(i)=Qsel((i-1)*supersam+1);
   end
   
   threshold=0.2;                     %�о�����
   for  i=1:data
       if Isam(i)>=threshold
           Ifinal(i)=1;
       else
           Ifinal(i)=-1;
       end
       if Qsam(i)>=threshold
           Qfinal(i)=1;
       else
           Qfinal(i)=-1;
       end
   end
   
k4 = 1;
for k1=1:2:len_spread-1 %ʵ���鲿�ŵ�һ������receive��  
    receive(k1) = Ifinal(k4);
    receive(k1+1) = Qfinal(k4);
    k4 = k4+1;
end