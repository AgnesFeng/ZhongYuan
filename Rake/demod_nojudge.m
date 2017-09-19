function  [receive] = demod_nojudge(len_spread,data,t,tap)
supersam = 5;
NT = 50;
N = 500;
fc = 5e6;
fs = 25e6;
rf = 0.1; 
Idem = ones(1,supersam*data+N);
Qdem = Idem;
Isel = ones(1,supersam*data);
Qsel = Isel;
Isam = ones(1,data);
Qsam = Isam;
Ifinal = ones(1,data);
Qfinal = Ifinal;
receive = len_spread;
for i=1:supersam*data+N
       Idem(i)=tap(i)*sqrt(2)*cos(2*pi*fc*t(i));
       Qdem(i)=tap(i)*(-sqrt(2)*sin(2*pi*fc*t(i)));
   end
   mtf=rcosfir(rf,NT,supersam,fs,'sqrt'); %匹配滤波器
   Imat=conv(Idem,mtf);
   Qmat=conv(Qdem,mtf);
 
   for  i=1:supersam*data                 %data selection
       Isel(i)=Imat(i+N);
       Qsel(i)=Qmat(i+N);
   end
              
   for i=1:data                       %提取码元   星座图
       Isam(i)=Isel((i-1)*supersam+1);
       Qsam(i)=Qsel((i-1)*supersam+1);
   end
   
k4 = 1;
for k1=1:2:len_spread-1 %实部虚部放到一个数组receive中  
    receive(k1) = Isam(k4);
    receive(k1+1) = Qsam(k4);
    k4 = k4+1;
end