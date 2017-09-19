clear all;
close all;
clc;
%-----------------------1��ԭʼ�ź�-----------------------------------------
Tlen = 800; %���ݳ���  �������������9�Ĺ�ϵ
s_initial = randsrc( 1, Tlen );
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%����
% R = 2/3;              % code rate
% g{1} = [1 0 1 1];                % Impulse Responses _ 1
% g{2} = [1 1 0 1];                % Impulse Responses _ 2
% g{3} = [1 1 1 1];                % Impluse Responses _ 3
% %n = length(g);                       % Convolution Code (1/n) parameter,n=3
% memory_els = 3;%K=4�Ĵ�������,m=K-1,Լ������
% initial = round(rand(1,Tlen));            % message vector��Ϣλ���ȣ�L
% %ENCODER: 1/3 Convolution Encoder%
% s_initial = encode_1_3(initial,g,length(g));%m����Ϣλ(1*L),g�����ɾ��������Ѿ�������n�����ʵĵ�����
% for i = 1:length(s_initial)%��˫����
%     if s_initial(i) == 0;
%         s_initial(i) = -1;
%     end
% end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%����
v = 225000; %Ҫ��Ӧ����225KBPS,��Ӧ4.4΢�룬�˴���0.004΢�룬���Խ���ʱ����������������0,001΢��~0.003΢��
Tb = 1/v; %��Ԫ����ʱ��    (��øĳ���ʵ���ݡ�225000)

%-----------------------4������--------------------------------------------
%������Ϊ�޳���ʱ���
downData= zeros(1,length(spreadData)/Tm_sample);
j = 1;
for i=1:Tm_sample:length(spreadData)-9
    downData(j) = spreadData(i);
    j = j+1;
end
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
P_OR = cha3_mrake(sample,length(txx)); %������Ƶ�ʣ����г��ȣ����У�������
pp1 = P_OR(1,:); %ÿһ�ױ���   length(rx)����
pp2 = P_OR(2,:);
pp3 = P_OR(3,:); 
EbN0db  = 0:1:15;
for snr = 1:length( EbN0db ) 
   rx = awgn(txx,snr,'measured'); 
k11 = 1;
    tap1 = pp1.* rx;
    tap2 = pp2.* rx;
    tap3 = pp3.* rx;   

%-----------------------5�����--------------------------------------------
len_spread = length(spread);
receive1 = demod_nojudge(len_spread,data,t,tap1);
receive2 = demod_nojudge(len_spread,data,t,tap2);
receive3 = demod_nojudge(len_spread,data,t,tap3);
%��ԭ���г���ʱ�����Ƭ,ʹ�ظ�
r1 = REPEAT(receive1,Tm_sample); %10
r2 = REPEAT(receive2,Tm_sample);
r3 = REPEAT(receive3,Tm_sample);
%--------------

    %--------------------���----------------------%     ע�⣬ÿ���źŶ�Ҫ���
    tap = tap1 + tap2 + tap3;
    receive = demod(len_spread,data,t,tap) ;
   
    m =0;m1=0;m2=0;m3=0;
    for i= 1:len_spread
          if receive(i)~= spread(i)
            m= m+1;
          end
    end
    ber(k11) = m/len_spread;
    k11 =k11 + 1; %��ÿһ������ȶ�Ӧ��������д��һ������
end
figure
semilogy(0:15,ber,'-m+');hold on;


