clear all;
%close all;
clc;
%-----------------------1��ԭʼ�ź�-----------------------------------------
Tlen = 40009; %���ݳ���  �������������9�Ĺ�ϵ
%s_initial = randsrc( 1, Tlen );
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%����
R = 1/3;              % code rate
g{1} = [1 0 1 1];                % Impulse Responses _ 1
g{2} = [1 1 0 1];                % Impulse Responses _ 2
g{3} = [1 1 1 1];                % Impluse Responses _ 3
%n = length(g);                       % Convolution Code (1/n) parameter,n=3
memory_els = 3;%K=4�Ĵ�������,m=K-1,Լ������
initial = round(rand(1,Tlen));            % message vector��Ϣλ���ȣ�L
%ENCODER: 1/3 Convolution Encoder%
s_initial = encode_1_3(initial,g,length(g));%m����Ϣλ(1*L),g�����ɾ��������Ѿ�������n�����ʵĵ�����
for i = 1:length(s_initial)%��˫����
    if s_initial(i) == 0;
        s_initial(i) = -1;
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%����
v = 225000; %Ҫ��Ӧ����225KBPS,��Ӧ4.4΢�룬�˴���0.004΢�룬���Խ���ʱ����������������0,001΢��~0.003΢��
Tb = 1/v; %��Ԫ����ʱ��    (��øĳ���ʵ���ݡ�225000)

%-----------------------4������--------------------------------------------

spread = s_initial;
%����ת��
fc=300e6; %�ز�
len_spread = length(spread);
tx = ones(1,len_spread/2);%����ת���������
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
%���������
supersam=5;         %��IFFT�����൱�ڽ���Ƶ���ź���0���ǵ�ͨ��һ���ְɡ�
data = length(tx);  
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
%�������
NT=50;
N=2*supersam*NT;    % N=500
fs=1500e6; %������ Ϊ����������������
rf=0.2;
psf=rcosfir(rf,NT,supersam,fs,'sqrt');% psf��СΪ500
Ipulse=conv(Iinit,psf);
Qpulse=conv(Qinit,psf);
%Ƶ�װ���
for i=1:supersam*data+N   %��������Ŀ�ı� ����Ϊ�����Ե�ʣ�
    t(i)=(i-1)/(fs);      %����Ƶfc���Թ������ʣ�ÿ���ŵĲ���������=�����ʡ�
    Imod(i)=Ipulse(i)*sqrt(2)*cos(2*pi*fc*t(i));  %�������䣬������֮��ľ�������˲��������1/fs
    Qmod(i)=Qpulse(i)*(-sqrt(2)*sin(2*pi*fc*t(i)));
end
txx=Imod+Qmod;
%-----------------------4~5�����ŵ�----------------------------------------
sample = 1/fs;         %����Ƶ��100
psl = c_4(sample,length(txx)); %������Ƶ�ʣ����г��ȣ����У�������
EbN0db  = 0:1:8;
k11 = 1;
for snr = 1:length( EbN0db ) 
   rx = awgn(txx,snr,'measured'); 
   
    tap= psl.* rx;
    %-----------------------5�����--------------------------------------------
    len_spread = length(spread);
    receive = demod_nojudge(len_spread,data,t,tap);
    %-----------------------5������--------------------------------------------
    [final1,node] = decode_1_3(receive,length(g),memory_els,Tlen,0);%���������˫����
    Bit_Error_Number = length(find(final1(1:Tlen) ~= initial(1:Tlen)));     %��ԭʼ�źűȽ�
    Bit_Error_Rate(snr) = Bit_Error_Number/Tlen;
%     %%%%�о�
%     receive = ones(1,len_spread);
%     for i= 1:len_spread     
%           if receive_sum(i)> 0
%              receive(i) = 1;
%           else
%              receive(i) = -1;
%           end
%     end
%     %û���룬��������
%     m =0;
%     for i= 1:len_spread
%           if receive(i)~= spread(i)
%             m= m+1;
%           end
%     end
%     ber(k11) = m/len_spread;
%     k11 =k11 + 1; %��ÿһ������ȶ�Ӧ��������д��һ������
end
figure
%semilogy(0:30,ber,'-m+');hold on;
%semilogy(0:30,ber,'-bo');hold on;
%semilogy(0:30,ber,'-r*');hold on;
%semilogy(EbN0db,Bit_Error_Rate,'-m+');hold on; 
semilogy(EbN0db,Bit_Error_Rate,'-bo');hold on; 
%semilogy(EbN0db,Bit_Error_Rate,'-r*');hold on; 
xlabel('SNR'),ylabel('������');




