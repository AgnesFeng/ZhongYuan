
Tlen = 2000;
s_spread = randsrc( 1, Tlen ); %����Դ (1*8000�ľ���-1,1��ռһ��)
%--------------------1����,ǰ50������----------------------%
N_sample=100; %��ͼʱ�����ʴ�һ��ã��������ʵ�100�����������˲�ʱ������Ϊ5���Ϳ���
N = 50;
selq = SAM(s_spread,N_sample,N);%��������ǻ�ͼ�õģ�������࣬������˲������ʱ�õ���һ�ַ�������0������������Ҳ����

%--------------------2�ֳ���ͨ��----------------------%
%����֮���ٲ���Ч�������������һ�������ǲ�����Ϊ��ͼ�ˣ����Ƿ������ز���ˣ�ע�����ʱ��������Ҫһ��
%��ʡ���˲�ͬ������ȣ�
len_spread = length(s_spread);
receive1 = ones(1,len_spread);%�����Ľ�������
tx_re = ones(1,len_spread/2); 
tx_im = ones(1,len_spread/2);
k2 = 1;
for k1=1:2:len_spread-1 %����һ��  
    tx_re(k2) = s_spread(k1);
    tx_im(k2) = s_spread(k1+1);
    k2 = k2+1;
end
selq_re = SAM(tx_re,N_sample,Tlen/2);
selq_im = SAM(tx_im,N_sample,Tlen/2);
figure(1)
subplot(211);PSD222(1024,N_sample,selq,N);title('ԭ���й�����');%���ﻭͼ�Ĺ�����ָ�Ķ����ܶȣ��������ǹ��ʣ���λ���ߣ��������ܶ�������������
subplot(212);PSD222(1024,N_sample,selq_re,N);title('2���Ʊ�4���ƺ�Iͨ��������'); 

%---------------------��������˲���ͬʱ����ͨ��1����������2�������˲�����3�����-------------%
%1��������
zero=5;         
I=tx_re;
Q=tx_im;
Izero = SymbolToWaveform(I,zero);
Qzero = SymbolToWaveform(Q,zero);

%2�������˲���
NT=50;           %�˲�������Ϊ2*NT+1,���ݽ������
N_fir=2*NT*zero;    %�˲�������-1����length(psf)-1һ��
fc = 5e6;          %�ز�Ƶ��
fs = 25e6;         %����Ƶ��
Ts = 1e-6;         %��Ԫ���ȣ�SAM�����
rf=0.1;
psf=rcosfir(rf,NT,zero,Ts,'sqrt');% ��������rf���˲�������NT���������ʣ���������Ĳ�������������Ԫʱ��
Ipulse=conv(Izero,psf);%����󳤶ȸı�
Qpulse=conv(Qzero,psf);
figure(2)
subplot(321);plot(psf);title('�˲���ʱ����Ӧ');
axis([200  300  -0.2  0.6]);grid  on;
set(gca,'XTick',0:250:501);set(gca,'XTicklabel',{'-250','0','250'});
subplot(322);plot(20*log(abs(fftshift(fft(psf)))));title('�˲���Ƶ����Ӧ');
axis([0  N_fir  -350 50]);grid on;
set(gca,'XTick',0:250:500);set(gca,'XTicklabel',{'-250','0','250'});
subplot(323);plot(I(1:100));axis([0  100  -2 2]);title('�˲�ǰ��ͨ��ʱ����Ӧ');
subplot(324);plot(Ipulse(1:600));title('�˲�֮��ͨ��ʱ����Ӧ');
subplot(325);plot(20*log(abs(fftshift(fft(I)))));title('�˲�ǰ��ͨ��Ƶ����Ӧ');
axis([0 1000 -50 100]);grid  on;
set(gca,'XTick',0:500:1000);set(gca,'XTicklabel',{'-500','0','500'});
subplot(326);plot(20*log(abs(fftshift(fft(Ipulse)))));title('�˲�֮��ͨ��Ƶ����Ӧ');
axis([0 5500 -300 100]);grid  on;
set(gca,'XTick',0:2750:5500);set(gca,'XTicklabel',{'-2700','0','2700'});
%3�����ز�
dt = 1/fs;         %ʱ��������  
T = 1;            %�ź�ʱ��                 
tt = 0:dt:T;
t1 = tt(1:length(Ipulse));
cos1 = cos(2*pi*fc*t1); 
sin1 = sin(2*pi*fc*t1);
Imod=Ipulse.*(sqrt(2)*cos1);
Qmod=Qpulse.*(sqrt(2)*sin1);
sum=Imod+Qmod;

figure(3)%���������ܶ�
plot(20*log(abs(fftshift(fft(sum)))));title('���ƺ�˫ͨ��������(�ӳ��κ��˲�)');
axis([0 5500 -200 100]);grid  on;
set(gca,'XTick',0:2750:5500);set(gca,'XTicklabel',{'-2700','0','2700'});
%[Pxx,f]=psd(sum,1024,fs,window,noverlap,dflag);������Ƶ�ײ��ÿ�

%--------------------�����������ز�(���˲�)---------------------%
% t2 = tt(1:length(selq_re));
% cos2 = cos(2*pi*fc*t2); 
% sin2 = -sin(2*pi*fc*t2);
% tx_re1 = selq_re .* cos2;   %�����˲�
% tx_im1 = selq_im .* sin2;
% tx1 = tx_re1 + tx_im1;
%--------------------------------------------------------------------%

% %-------------------���ŵ�----------------------%
P_OR = cha_3(1/Ts,length(sum),sum);          %������Ƶ�ʣ����г��ȣ����У�������
pp1 = abs(P_OR(1,:)); %ÿһ�ױ���   length(rx)���� С�߶ȼӴ�߶ȵĽ��
pp2 = abs(P_OR(2,:));
pp3 = abs(P_OR(3,:));
k11 = 1;
j=1;
snrn = 15;
for snr = 0:snrn
     %rx = awgn(tx1,snr,'measured'); %��˹�ŵ�(ֻ���Ǹ�˹����) ��������
     rx = awgn(sum,snr,'measured'); %��˹�ŵ�(ֻ���Ǹ�˹����) �������� 
     tap1 = pp1.* rx;
     tap2 = pp2.* rx;
     tap3 = pp3.* rx; 
%     %����������������--��������˲���----------------------%    
%     re_2=rx.* (sqrt(2)*cos2);
%     im_2=rx.* (sqrt(2)*sin2);
%     re_21=rx.* re_2;
%     im_21=rx.* im_2;
%     for k3 = 1:length(rx)
%         if re_2(k3) > 0;
%             re_2(k3) = 1;
%         else
%             re_2(k3) = -1;
%         end
%         if im_2(k3) > 0;
%             im_2(k3) = 1;
%         else
%             im_2(k3) = -1;
%         end
%     end
%        %����
%        for i=1:length(I)                       
%        I2(i)=re_2((i-1)*N_sample+1);
%        Q2(i)=im_2((i-1)*N_sample+1);
%        end
%     k4 = 1;
%     for k1=1:2:len_spread-1 %ʵ���鲿�ŵ�һ������receive��  
%         receive(k1) = I2(k4);
%         receive(k1+1) = Q2(k4);
%         k4 = k4+1;
%     end
%     %-----------��������ͼ-------------------------%
%     selq_rece = SAM(receive,N_sample,N);
%     % figure(2)
%     subplot(211);PSD222(1024,N_sample,tx1,Tlen/2);
%     axis([-4000, 4000, -60, 20]);title('���ƺ��źŹ�����(δ���˲�)')
%     subplot(212);PSD222(1024,N_sample,selq_rece,N);
%     axis([-4000, 4000, -60, 20]);title('��������ף������˲��ģ�')
%     m=0;
%     for i= 1:length(I)
%       if (receive~=s_spread)
%         m= m+1;
%       end
%     end
%     ps(j) = m/length(receive)        %���һ������ȵ�������
%     j =j+1;                          %��ÿһ������ȶ�Ӧ��������д��һ������
% end
%--------------------------------------����������������������������%   

%����������--------------------��������˲���----------------------%   
% %Ƶ�װ���
data = length(I);
Idem=tap1.*(sqrt(2)*cos1);
Qdem=tap1.*(sqrt(2)*sin1);

%���˲�
Imat=conv(Idem,psf);
Qmat=conv(Qdem,psf);   
supersam=zero;
   for  i=1:supersam*data                 %�˲������
       Isel(i)=Imat(i+N_fir);
       Qsel(i)=Qmat(i+N_fir);
   end
          
   for i=1:data                       %��ȡ��Ԫ,ȥ��0
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
        receive1(k1) = Ifinal(k4);
        receive1(k1+1) = Qfinal(k4);
        k4 = k4+1;
    end
    %------------------------�ڶ����������Ľ��-----------------------%
    tap = tap1 + tap2 + tap3;
    receive_sum = demod2(psf,cos1,sin1,Tlen,Tlen/2,tap);
    receive2 = demod2(psf,cos1,sin1,Tlen,Tlen/2,tap2); %����ԭ����,��ͨ���󳤶�,ʱ���ᣬ��������
    receive3 = demod2(psf,cos1,sin1,Tlen,Tlen/2,tap3); %����ԭ����,��ͨ���󳤶�,ʱ���ᣬ��������

    m1=0;m2=0;m3=0;m =0;
    for i= 1:len_spread
          if receive_sum(i)~= s_spread(i)
            m= m+1;
          end
            if receive1(i)~= s_spread(i)
            m1= m1+1;
          end
          if receive2(i)~= s_spread(i)
            m2= m2+1;
          end
          if receive3(i)~= s_spread(i)
            m3= m3+1;
          end
    end
    ber(k11) = m/len_spread;
    ber1(k11) = m1/len_spread;%���һ������ȵ�������
    ber2(k11) = m2/len_spread;
    ber3(k11) = m3/len_spread;
    k11 =k11 + 1;    
end

figure
%semilogy(0:15,berno,'-b*');grid on; hold on;
semilogy(0:snrn,ber,'-m*');grid on; hold on;
semilogy(0:snrn,ber1,'-b*');grid on; hold on;
semilogy(0:snrn,ber2,'-go');hold on;
semilogy(0:snrn,ber3,'-r*');hold on;
legend('��������','��һ��','�ڶ���','������');
xlabel('�����');ylabel('������');
figure
subplot(321);
plot(20*log(abs(fftshift(fft(Imod)))));title('���ƺ�ͨ��������');
axis([0 5500 -200 100]);grid  on;
set(gca,'XTick',0:2750:5500);set(gca,'XTicklabel',{'-2700','0','2700'});
subplot(322);
plot(20*log(abs(fftshift(fft(Idem)))));title('���--Ƶ�װ���--������');
axis([0 5500 -50 100]);grid  on;
set(gca,'XTick',0:2750:5500);set(gca,'XTicklabel',{'-2700','0','2700'});
subplot(323);
plot(20*log(abs(fftshift(fft(Imat)))));title('���--���˲���--Ƶ����Ӧ');
axis([0 6000 -300 130]);grid  on;
set(gca,'XTick',0:3000:6000);set(gca,'XTicklabel',{'-300','0','3000'});
subplot(324);
plot(20*log(abs(fftshift(fft(Isam)))));title('�����о���Ƶ����Ӧ');
axis([0 1000 -50 130]);grid  on;
set(gca,'XTick',0:500:1000);set(gca,'XTicklabel',{'-500','0','500'});
subplot(325);plot(Imat(1:800));title('���--���˲���--ʱ����Ӧ');
subplot(326);plot(Ifinal(1:100));axis([0  100  -2 2]);title('�����о�֮��ͨ��ʱ����Ӧ');


