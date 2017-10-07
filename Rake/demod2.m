function  [receive] = demod2(psf,cos1,sin1,len_spread,data,tap) %序列原长度,分通道后长度,时间轴，输入序列
supersam = 5;
NT=50;
N=2*NT*supersam;    %滤波器长度-1，与length(psf)-1一样
mtf=psf;
%乘载波
Idem=tap.*(sqrt(2)*cos1);
Qdem=tap.*(sqrt(2)*sin1);

%过滤波
Imat=conv(Idem,mtf);
Qmat=conv(Qdem,mtf);   
 
   for  i=1:supersam*data                 %data selection
       Isel(i)=Imat(i+N);
       Qsel(i)=Qmat(i+N);
   end
              
   for i=1:data                       %提取码元   
       Isam(i)=Isel((i-1)*supersam+1);
       Qsam(i)=Qsel((i-1)*supersam+1);
   end
   
   threshold=0.2;                     %判决门限
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
for k1=1:2:len_spread-1 %实部虚部放到一个数组receive中  
    receive(k1) = Ifinal(k4);
    receive(k1+1) = Qfinal(k4);
    k4 = k4+1;
end