clear all;
close all;

Numusers = 1;
Nc = 16; %��Ƶ����
ISI_Length = 1;%ÿ����ʱΪISI_Length/2
EbN0db = [0:1:15];%����ȣ���λdb
Tlen = 8000; %���ݳ��� ��ȡԽ��Խ�ã�100000
Bit_Error_Number1 = 0;%������ʵĳ�ʼֵ
Bit_Error_Number2 = 0;
Bit_Error_Number3 = 0;
power_unitary_factor1 = sqrt( 5/9 ); %ÿ����������
power_unitary_factor2 = sqrt( 3/9 );
power_unitary_factor3 = sqrt( 1/9 );
s_initial = randsrc( 1, Tlen ); %����Դ (1*8000�ľ���-1,1��ռһ��)

% ����Ԫ����ʱ�䣬���㻭�����ܶ��׺���ͼ
v = 112500;
Ts = 1; % ��Ԫ����(��ͼ��)
N_sample = 100; % ������Ԫ��������
dt = Ts / N_sample; % ����ʱ����
N = 50; % ��Ԫ��
t = 0 : dt : (N * N_sample - 1) * dt; % ���д���ʱ��
gt1 = ones(1, N_sample); % �����㲨��
% �����������
RAN = round(s_initial); % ���0 1����
selq = [];
for i = 1 : N % ��������
   if RAN(i)==1
       selq = [selq gt1]; %��100���㶼��1����ʾһ����Ԫ
   else
       selq = [selq -1*gt1];%��100���㶼��-1����ʾһ����Ԫ
   end
end
figure(111)
subplot(2, 1, 1);plot(t, selq);grid on;title('ԭʼ����ʱ����');
axis([0 10 -2 2])
t3=t;after=selq;
% �������ܶȼ��㷽��һ
fft_se1 = fftshift(fft(selq)); % �����е�Ƶ��
PE1 = 10 * log10(abs(fft_se1) .^ 2 / (N * Ts)); % ��ʽ����������ܶ�
PEL1 = (-length(fft_se1) / 2 : length(fft_se1) / 2 - 1) / 10; % �����䳤��
%%���Ƴ����
figure(222)
subplot(2, 1, 1);plot(PEL1, PE1); grid on; title('ԭʼ����PSD');
axis([-50 50 -50 100]);xlabel('f');ylabel('����');

%%%%%%%%%%%%%%%%%%%%%%%%����walsh����%%%%%%%%%%%%%%%%%%%%%%
wal2 = [ 1 1; 1 -1 ];
wal4 = [wal2 wal2; wal2 wal2*(-1)];  
wal8 = [wal4 wal4; wal4 wal4*(-1)];  %8*8
wal16 = [wal8 wal8; wal8 wal8*(-1)]; %16*16
%%%%%%%%%%%%%%%%%%%%%%%%%%%��Ƶ%%%%%%%%%%%%%%%%%%%%%%%%%%%%
s_spread = zeros( Numusers, Tlen*Nc );  %[1,8000*16]��Ƶ�����Ƭ ��ʼ��
ray1 = zeros( Numusers, 2*Tlen*Nc ); 
ray2 = zeros( Numusers, 2*Tlen*Nc );
ray3 = zeros( Numusers, 2*Tlen*Nc );

for i = 1:Numusers
    x0 = s_initial( i,: ).'*wal16( 8,: );  %����ת�ã��б��У�����8*8����ĵڰ���
    x1 = x0.';%���16*8000�ľ���
    s_spread( i,: ) = ( x1(:) ).';%�õ���Ƶ������У���x1(:)��ͷ��β���򣬼�16*8000����    ��ʱ��ͼ��s_spread,����Ԫ���ڣ���������ͼ
end

% ����Ԫ����ʱ�䣬���㻭�����ܶ��׺���ͼ
Tc = Ts/16; % ��Ƭ����
dt = Tc / N_sample; % ����ʱ����
N = N*16; % ��Ƭ��
t = 0 : dt : (N * N_sample - 1) * dt; % ���д���ʱ��
gt1 = ones(1, N_sample); % �����㲨��
% �����������
RAN = round(s_initial); % ���0 1����
selq = [];
for i = 1 : N % ��������
   if RAN(i)==1
       selq = [selq gt1]; %��100���㶼��1����ʾһ����Ԫ
   else
       selq = [selq -1*gt1];%��100���㶼��-1����ʾһ����Ԫ
   end
end
figure(111)
subplot(2, 1, 2);plot(t, selq);grid on;title('��Ƶ������ʱ����');xlabel('t');ylabel('����');
axis([0 10 -2 2]);
t2=t; before=selq;
N2 = 800;
RAN2 = round(s_spread); % ���0 1����
selq2 = [];
for i = 1 : N2 % ��������
   if RAN2(i)==1
       selq2 = [selq2 gt1]; %��100���㶼��1����ʾһ����Ԫ
   else
       selq2 = [selq2 -1*gt1];%��100���㶼��-1����ʾһ����Ԫ
   end
end
% �������ܶȼ��㷽��һ
fft_se12 = fftshift(fft(selq2)); % �����е�Ƶ��
PE12 = 10 * log10(abs(fft_se12) .^ 2 / (N2 * Ts)); % ��ʽ����������ܶ�
PEL12 = (-length(fft_se12) / 2 : length(fft_se12) / 2 - 1) / 10; % �����䳤��
%%���Ƴ����
figure(222);subplot(2, 1, 2);plot(PEL12, PE12); grid on; title('��Ƶ�����ܶ���');
axis([-1000 1000 -350 70]);xlabel('f');ylabel('����');

figure(2111)
subplot(211);selq1 = SAM(s_initial,N_sample,N);%axis([-1000 1000 -350 70])
sing1=PSD2(1024,N_sample,selq1,N,Ts);title('ԭ���е���Ԫ������');
subplot(212);selq2 = SAM(s_spread,N_sample,32*N);%axis([-1000 1000 -350 70])
sing2=PSD2(1024,N_sample,selq2,N,Ts/16);title('��Ƶ����Ԫ������');


%��ÿ����Ƶ�������ظ����Σ�������������ӳ٣��ӳٰ����Ԫ��
ray1( 1:2:2*Tlen*Nc - 1 ) = s_spread( 1:Tlen*Nc );         %���ڵļ����ظ���
ray1( 2:2:2*Tlen*Nc ) = ray1( 1:2:2*Tlen*Nc - 1 );         %ray1�ǵ�һ���ź�
spread = ray1;                                            %��Ƶ��ķ����ź�
%--------------------����----------------------%
%��ʡ���˲�ͬ������ȣ�
fc=5e6; %�ز�Ƶ��
len_spread = length(spread);
tx = ones(1,len_spread/2);%���ƺ������
% receive1 = ones(1,len_spread);%�����Ľ�������   ��%�������ռ������������
% receive2 = ones(1,len_spread);
% receive3 = ones(1,len_spread);
tx_re = tx; 
tx_im = tx;
k2 = 1;
for k1=1:2:len_spread-1 %����һ�� �� �õ�tx_re,im_re,����16*8000 = 64000
    tx_re(k2) = spread(k1);
    tx_im(k2) = spread(k1+1);
    k2 = k2+1;
end
I = tx_re;
Q = tx_im;
% zero insertion��������룬����ֱ�ۻ�ͼ���˹��̳�Ϊ���Ρ����ε���˼����ʵ������Ϣ�����ε�ת�����Ա㷢�䣬���塰���Ρ�Ӧ�����ڻ�������֮��
supersam=5;         %sampling  rate  25M HZ  ,supersamΪ�������ʡ������� ������fs/�����ʡ�
data = length(tx);   % data=128000
nnn = supersam*data;  % ��������Ŀ=��������*ԭ����Ŀ
for  ii=1:nnn     
    if rem(ii,supersam)==1        %ȡ����
        tem1 = fix((ii-1)/supersam)+1;  %fix:��0����ȡ��
        tem2 = fix((ii-1)/supersam)+1;
        Iinit(ii)=I(tem1);  %Iinit����Ӧ�õ���data
        Qinit(ii)=Q(tem2);
    else
         Iinit(ii)=0;
         Qinit(ii)=0;
    end
end
%��������˲����� ���ţ������е�ͨ�˲�����Ϊ ���Ŵ������ʵ����󣬻��������Ƶ�׽����
%������˲������������˲������е�ͨ�˲����������Ƶ��ʱ����ܻ�������ѡ�
%ƽ�����������˲���
%psf=rcosfir(rf,n_t,rate,fs,'sqrt')   rate:�������ʣ�rf:�������ӣ�n_t:�˲���������fs:������
%���ڵ��ƻ���֮ǰ�����ڽ�������֮���������͹���������������������ISI����䴮�ţ�
NT=50;
N=2*supersam*NT;    % N=500
fs=25e6;
rf=0.3;
psf=rcosfir(rf,NT,supersam,fs,'sqrt');% psf��СΪ500
Ipulse=conv(Iinit,psf);
Qpulse=conv(Qinit,psf);
%modulation
for i=1:supersam*data+N   %��������Ŀ�ı� ����Ϊ�����Ե�ʣ�
    t(i)=(i-1)/(fs);      %������Ϊ������Ƶ�������ʴ�С��ȣ���������Ƶfc���Թ�������=�����ʡ�
    Imod(i)=Ipulse(i)*sqrt(2)*cos(2*pi*fc*t(i));
    Qmod(i)=Qpulse(i)*(-sqrt(2)*sin(2*pi*fc*t(i)));
end
tx=Imod+Qmod;
%�������εĲ�
% x1 = 0:1:4999;
% figure
% txx = tx(1:5000);      %��ѡ5000����
% plot(x1,txx)
% xlabel('��ѡǰ5000��');ylabel('�źŷ���')
% title('���ƺ�ĳ��β�')
%    figure
%    plot(20*log(abs(fft(s_initial))));
%    axis([0  data  -40  100]);
%    grid on;

ttt = 0:length(tx);
nfft = 1024;
cxn = xcorr(s_initial,'unbiased');%�������е�����غ���
CXk = fftshift(fft(cxn,nfft));
Pxx = abs(CXk);
%index = 0:round(nfft/2-1);
%k = index*1/nfft;
%plotPxx = 10*log10(Pxx(index+1));
in = (-length(CXk) / 2 : length(CXk) / 2 - 1) / 10; % �����䳤��
plotPxx = 10*log10(Pxx);
% figure(111)
% subplot(313)
% plot(in,plotPxx);%��ƽ��һ��ͺã�
% xlabel('Ƶ��');ylabel('����');
% title('��һ�����ƺ����ܶ��׺���')


%-------------------���ŵ�----------------------%
sample = 1/fs;         %����Ƶ��100
P_OR = cha_3(sample,length(tx),tx);          %������Ƶ�ʣ����г��ȣ����У�������
pp1 = abs(P_OR(1,:)); %ÿһ�ױ���   length(rx)����
pp2 = abs(P_OR(2,:));
pp3 = abs(P_OR(3,:)); 
k11 = 1;
for snr = 1:length( EbN0db ) 
    rx = awgn(tx,snr,'measured');            %��˹�ŵ�(ֻ���Ǹ�˹����) ��������
    tap1 = pp1.* rx;
    tap2 = pp2.* rx;
    tap3 = pp3.* rx;   
%�����������������������������ʱ��ͼ��-----------------------
% tt1 = tap1(1:5000);      %��ѡ5000����
% tt2 = tap2(1:5000);
% tt3 = tap3(1:5000);
% figure
% subplot(311)
% plot(x1,tt1);
% ylabel('��1���źŷ���')
% title('�źž���˥���Ĵ�������');
% subplot(312)
% plot(x1,tt2);ylabel('��2���źŷ���')
% subplot(313)
% plot(x1,tt3);xlabel('ǰ5000����');ylabel('��3���źŷ���')
%--------------------���������ʱ��ͼ�á�����������������������

    %--------------------���----------------------%     ע�⣬ÿ���źŶ�Ҫ���
    receive1 = demod(len_spread,data,t,tap1) ;
    receive_2 = demod(len_spread,data,t,tap2) ;
    receive2( ISI_Length + 1:len_spread) = receive_2( 1:len_spread - ISI_Length ); %���Ե�һ���źŽ���һ���ӳ٣�ģ����ʱ��ĵڶ���
    %receive2 = receive_2;%����ʱ�õ�����һ��
    
    receive_3 = demod(len_spread,data,t,tap3) ;
    receive_3( ISI_Length + 1:len_spread ) = receive_3( 1:len_spread - ISI_Length );%���Ե�һ���źŽ��������ӳ٣�ģ����ʱ��ĵڶ����͵�����
    receive3( 2*ISI_Length + 1:len_spread ) = receive_3( 1:len_spread - 2*ISI_Length );
    %receive3 = receive_3;%����ʱ�õ�����һ��,�����ʷ�������ˣ�˵�����ն���ʱ����Ĳ���
    
%     m1=0;m2=0;m3=0;             %��������������
%     for i= 1:len_spread
%       if receive1(i)~= spread(i)
%         m1= m1+1;
%       end
%       if receive2(i)~= spread(i)
%         m2= m2+1;
%       end
%       if receive3(i)~= spread(i)
%         m3= m3+1;
%       end
%     end
%     ber1(k11) = m1/len_spread;  %���һ������ȵ�������
%     ber2(k11) = m2/len_spread;
%     ber3(k11) = m3/len_spread;
%     k11 =k11 + 1; %��ÿһ������ȶ�Ӧ��������д��һ������

% figure
% semilogy(0:10,ber1,'-b*');grid on; hold on;
% semilogy(0:10,ber2,'-g*');hold on;
% semilogy(0:10,ber3,'-r*');hold on;

ttt = 0:length(receive1);
nfft = 1024;
cxn = xcorr(receive1,'unbiased');%�������е�����غ���
CXk = fftshift(fft(cxn,nfft));
Pxx = abs(CXk);
%index = 0:round(nfft/2-1);
%k = index*1/nfft;
%plotPxx = 10*log10(Pxx(index+1));
in = (-length(CXk) / 2 : length(CXk) / 2 - 1) / 10; % �����䳤��
plotPxx = 10*log10(Pxx);
% figure(2)
% subplot(211)
% plot(in,plotPxx);%��ƽ��һ��ͺã�
% xlabel('Ƶ��');ylabel('����');
% title('��һ����������ܶ��׺���')

%--------------------����----------------------
%for nEN = 1:length( EbN0db )    
    en = 10^( EbN0db(snr)/10 ); %��Eb/N0��dbֵת����ʮ������ֵ
    sigma = sqrt(32/(2*en));
    %sigma = sqrt( 32/(2*en) );%���ݳ�ʼ״̬�趨������ȣ������AWGN���32Ϊ��������������Ƶ���Ʒ�ʽ�йأ���Ƶ��ÿ��32bit��Ӧһ��ԭʼbit��������������ԭ����32��
    %���յ����ź�demp
    demp = power_unitary_factor1*receive1+...           %ÿ���˹��������ٵ���
    power_unitary_factor2*receive2+...                  %ע����һ��rand��randn������
    power_unitary_factor3*receive3 + ( 0.2*rand( 1,len_spread )+ 0.8*randn( 1,len_spread )*Numusers )*sigma;          %���һ��Ϊ����     %rand:0~1֮�������             randn:��̬�ֲ��������
    dt = reshape( demp,32,Tlen )';   %dtΪ8000*32��dempֵ�ľ���
    %��walsh���ظ�����
    wal16_d(1:2:31) = wal16(8,1:16); %�øղ���Ƶ���������
    wal16_d(2:2:32) = wal16(8,1:16);
    
    %������rdata1Ϊ��һ�����
    rdata1 = dt*wal16_d(1,:).';  
    %��walsh���ӳٰ����Ƭ
    wal16_delay1(1,2:32) = wal16_d(1,1:31);
    wal16_delay1(1) = wal16_d(32);%�Լ��ӵģ�������0���������������ż���߽��棬��ʵǰ���ǲ���0����ν����Ϊȡ�ú��������ˣ��������ֲ���ȥ�����������������ܹ��ʣ�
    %������rdata2Ϊ�ڶ������
    rdata2 = dt*wal16_delay1(1,:).';
    %��walsh���ӳ�һ����Ƭ
    wal16_delay2(1,3:32) = wal16_d(1,1:30);%������ĵ�һ���±�Ϊ1����0��
    wal16_delay2(1,1:2) = wal16_d(1,31:32);%���������ֵ�ŵ���ǰ��������0��
    %������rdata3Ϊ���������
    rdata3 = dt*wal16_delay2(1,:).';
    
    p1 = rdata1'*rdata1;   %����ÿһ���Ĺ���
    p2 = rdata2'*rdata2;
    p3 = rdata3'*rdata3;
    p = p1 + p2 + p3;
    u1 = p1/p;             %��ÿһ���ļ�Ȩϵ��
    u2 = p2/p;
    u3 = p3/p; 
    %���Ⱥϲ�
    rd_m1 = real( rdata1*u1+rdata2*u2+rdata3*u3);
    %������ϲ�
    rd_m2 = (real(rdata1+rdata2+rdata3))/3;
    %ѡ��ʽ�ϲ�
    u = [u1,u2,u3];   %��ѭ����⣩
    maxu = max(u);
    if(maxu==u1)
        rd_m3 = real(rdata1);
      else if(maxu==u2)
           rd_m3 = real(rdata2);
      else
           rd_m3 = real(rdata3);
      end
    end %���ַ����о����
    r_Data1 = sign(rd_m1)'; %����0Ϊ1��С��0Ϊ-1
    r_Data2 = sign(rd_m2)';
    r_Data3 = sign(rd_m3)';
    %�����������
    Bit_Error_Number1 = length(find(r_Data1(1:Tlen) ~= s_initial(1:Tlen)));     %��ԭʼ�źűȽ�
    Bit_Error_Rate1(snr) = Bit_Error_Number1/Tlen;
    Bit_Error_Number2 = length(find(r_Data2(1:Tlen) ~= s_initial(1:Tlen)));
    Bit_Error_Rate2(snr) = Bit_Error_Number2/Tlen;
    Bit_Error_Number3 = length(find(r_Data3(1:Tlen) ~= s_initial(1:Tlen)));
    Bit_Error_Rate3(snr) = Bit_Error_Number3/Tlen;

end
    
figure(3333)
subplot(2, 1, 1);plot(t2, before);grid on;title('�о�������ʱ����');
axis([0 20 -2 2]);
subplot(2, 1, 2);plot(t3, after);grid on;title('����������ʱ����');
axis([0 20 -2 2]);
 
figure
semilogy(EbN0db,Bit_Error_Rate1,'ro-');hold on; 
semilogy(EbN0db,Bit_Error_Rate2,'bo-');hold on;
semilogy(EbN0db,Bit_Error_Rate3,'go-');hold on;
legend('���Ⱥϲ�','������ϲ�','ѡ��ʽ�ϲ�');
xlabel('�����');ylabel('������');
title('������Ҫ�ּ��ϲ���ʽ���ܱȽ�');
grid  on;