clear all;        %������һ�����ڵ�m���л�ͼЧ������
close all;
clc;
%%
%----------------------1��walsh����---------------------
wal2 = [ 1 1; 1 -1 ];
wal4 = [wal2 wal2; wal2 wal2*(-1)];  
wal8 = [wal4 wal4; wal4 wal4*(-1)];  %8*8
wal16 = [wal8 wal8; wal8 wal8*(-1)]; %16*16
p = wal16(5,:);
KP = length(p); %��Ƶ����
Tm_sample = 4;%������Խ�󣬽��ͼԽ��׼
pp = SAM_d(p,Tm_sample,length(p));%�����ڣ��г���ʱ��
%%
%%--------------------------2��m����----------------
m = m_sequence([0,0,0,1,1,1,0,1]); %�õ�һ�����ڵ�m���� (���Ĵ����ĳ�ʼ״̬��������ΪC1~Cn��c0Ĭ��Ϊ1) 
N_m = length(m);%m��:225
m1 =  SAM_m(m,Tm_sample,length(m)); %������ʱ��
m2 = 1-2*m1;  %��˫����
L_m = length(m2);
peri = 3;
mm = repmat(m,1,peri); %3�����ڵ�m���У����������
mm1 = SAM_m(mm,Tm_sample,length(mm));%�����ڣ��г���ʱ��
mm2 = 1-2*mm1;       %��˫����

rt1=conv(mm2(1:3*L_m),m2(end:-1:1))/(N_m*Tm_sample);%ȡ�м��һ����
index = L_m+Tm_sample;
figure(1);subplot(211);plot(rt1(index:end-index));title('������m���е�����غ���conv,���ֵ��λ�ò�����');
xcor = dsp.Crosscorrelator;
y = step(xcor,m2',m2'); 
figure(1);subplot(212);plot(y);title('������m���е�����غ��������������β�����');
%%
%---------------------3���ź����У���Ƶ,��ǰ������---------------------
Tlen = 3000; 
v = 225000;
Tb = 1/v; %��Ԫ����ʱ��  
Tc = Tb/KP;
dt = Tc/Tm_sample;
s_initial = randsrc( 1, Tlen );
T_sample = Tm_sample * KP;  %m���еĲ�������Ϊ4
s_initial1 =  SAM_d(s_initial,T_sample,length(s_initial));
pt = repmat(pp,1,Tlen);
spreadDate = pt.*s_initial1;%��Ƶ
spreadDate_m = ones(1,length(m2)+length(s_initial1));%��ǰ������
spreadDate_m(1:length(m2)) = m2;
spreadDate_m(length(m2)+1:end) = s_initial1;
len_spread = length(spreadDate_m);
%%
%�ҵ���λ����ʱ���±�
x1 = spreadDate_m(1:length(m2));
x2 = m2;
xcor = dsp.Crosscorrelator;
y = step(xcor,x2',x1'); %computes cross-correlation of x1 and x2
figure(2), subplot(211);plot(y); title('Correlated output');
aim = find(max(y) == y);

%%
%--------------------qpsk���ƣ�������ת�����޼��ر任-------------------
%spreadDate_m = rand(1,128);len_spread = length(spreadDate_m);
%fs = 1/dt;
tx = ones(1,len_spread/2);%����ת���������
tx_re = tx; 
tx_im = tx;
k2 = 1;
for k1=1:Tm_sample*2:len_spread-2*Tm_sample+1 
    tx_re(k2:k2+Tm_sample-1) = spreadDate_m(k1:k1+Tm_sample-1);
    tx_im(k2:k2+Tm_sample-1) = spreadDate_m(k1+Tm_sample:k1+2*Tm_sample-1);
    k2 = k2+Tm_sample;
end
txx = tx_re + 1i*tx_im;
%%
%-------------------�ŵ�����ѡ��-----------------
%channel_type = input('�������ŵ����ͣ�'); 
channel_type=2;

ts =dt;
if(channel_type==1)
    fd = 25;
    k = 10^(12/10);
    chan = ricianchan(ts,fd,k);  
elseif(channel_type==2)
    fd = 25;           %��������Ƶ�� 
    k = 10^(6/10);
    tau = [0 0.000001 0.000002];
    pdb = [0,-9,-12];
    chan = ricianchan(ts,fd,k,tau,pdb); %�ྶ��Ҫ��ÿһ����������˥�䡣ts�������źŵĲ���ʱ�䣨s����fd����������Ƶ�ƣ�pdb��ƽ��·�����棨dB��
elseif(channel_type==3)
    fd = 10;     %��������Ƶ�� 
    tau = [0 0.000001 0.000002 0.000003];
    pdb = [0,-3,-6,-9];
    chan = rayleighchan(ts,fd,tau,pdb);
else                                    % (channel_type~=3 &&channel_type~=1 &&channel_type~=2)
    error('Error! channel type should be one of 1 2 3!');
end 
%%
%--------------------���ŵ�-----------------
c_out = filter(chan,txx);
EbNo = 1:1:3;
for snr = 1:length(EbNo)
   y_snr = awgn(c_out,snr,'measured');
   %y_snr = c_out;
% psl = c_3(1/dt,500); %������Ƶ�ʣ����г��ȣ����У�������
% rxx = filter(psl,1,txx);
% EbN0db  = 0:1:8;
% for snr = 1:length( EbN0db ) 
%     y_snr = awgn(rxx,snr,'measured'); 
%%   
    %--------------------���������ת��-------------------
    re_2=real(y_snr);
    im_2=imag(y_snr);
    rx = ones(1,len_spread);
    k2 = 1; 
    for k1=1:Tm_sample*2:len_spread-2*Tm_sample+1 
        rx(k1:k1+Tm_sample-1) = re_2(k2:k2+Tm_sample-1);
        rx(k1+Tm_sample:k1+2*Tm_sample-1) = im_2(k2:k2+Tm_sample-1);
        k2 = k2+Tm_sample;
    end
   %�о�
   receive1 = sign(rx);%�о����������ķ�ֵ��һЩ
   BitErrorBeforeRake = length(find(receive1 ~= spreadDate_m)); 
   BitErrorBeforeRake_rate(snr) = BitErrorBeforeRake/(len_spread); 

%%
    %------------��ط岶��
    x2 = rx(1:length(m2));
    xcor = dsp.Crosscorrelator;
    buHuo = step(xcor,x2',m2'); %computes cross-correlation of x1 and x2
    y1 = abs(buHuo);
    figure(2), subplot(212);plot(y1); title('��������ط�');

    y_max = max(y1);
    id = find(y_max == y1); %�ҵ����ֵ�������ҿ���,��1��Ԫ����ʱ����4΢�����ڣ�
    L_catch = Tm_sample * KP*2 ;
    feng_x = ones(1,10);
    a = 1;
    b1 = id - L_catch/2;
    b2 = id + L_catch/2;
    for x =b1 : 1 : b2
        if y1(x) >= y_max/4 && y1(x-1) < y1(x) && y1(x) > y1(x+1) && y1(x-2) < y1(x) && y1(x) > y1(x+2)
            feng_x(a) = x;
            a = a+1;
        end
    end
    fengPaiXu_x = feng_x;
    for j = 1:a-1   %�����ֵ�ĺ����������ҳ���������ֵ�ĺ�����
        for i = 1:a-1-j
            if y1(fengPaiXu_x(i+1)) > y1(fengPaiXu_x(i))
              temp = fengPaiXu_x(i);
              fengPaiXu_x(i) = fengPaiXu_x(i+1);
              fengPaiXu_x(i+1) = temp;
            end
        end
    end
%%
    %----------------��λ����
    %%%%%%%%��׽��3��������
    if(feng_x(3)~=1)   
        shift_x = ones(1,3);  
        feng = ones(1,3);
        Conj = ones(1,3);  
        for i = 1:3
          shift_x(i) = fengPaiXu_x(i)-aim; 
          feng(i) = buHuo(fengPaiXu_x(i));
          Conj(i) = conj(feng(i));
        end
        path_1 = rx;
        path_2 = rx;
        path_3 = rx;
        path_1(1,1:end-shift_x(1)) = path_1(1,shift_x(1)+1:end);%����
        path_2(1,1:end-shift_x(2)) = path_2(1,shift_x(2)+1:end);
        path_3(1,1:end-shift_x(3)) = path_3(1,shift_x(3)+1:end);
        p = y1(fengPaiXu_x(1)) + y1(fengPaiXu_x(2)) + y1(fengPaiXu_x(3));
        %p = y1(tempxx(1)) + y1(tempxx(2)) + y1(tempxx(3))+y1(tempxx(4));   %����ÿһ���ļ�Ȩϵ��
        u1 = y1(fengPaiXu_x(1))/p;             
        u2 = y1(fengPaiXu_x(2))/p;
        u3 = y1(fengPaiXu_x(3))/p; 
        merge = path_1 * u1 + path_2 * u2 + path_3 * u3;
        %merge = path_1 * u1 + path_2 * u2 + path_3 * u3 + path_4 * u4; %Ȩ��Ҫ��Ӧ
    end
    %%%%%%%%%%%%%%%%%%��׽��2��
    if(feng_x(2)~=1 && feng_x(3)==1) 
      shift_x = ones(1,2);  
      feng = ones(1,2);
      Conj = ones(1,2);
      for i = 1:2
        shift_x(i) = fengPaiXu_x(i)-aim; 
        feng(i) = buHuo(fengPaiXu_x(i));
        Conj(i) = conj(feng(i));
      end    
      path_1 = rx;
      path_2 = rx;
      path_1(1,1:end-shift_x(1)) = path_1(1,shift_x(1)+1:end);%����
      path_2(1,1:end-shift_x(2)) = path_2(1,shift_x(2)+1:end);
      merge = path_1 * Conj(1) + path_2 * Conj(2);
    end
    %%%%%%%%%%%%%%%%%%%��׽��1��
    if(feng_x(1)~=1 && feng_x(3)==1 && feng_x(2)==1) 
      shift_x = fengPaiXu_x(1)-aim; 
      feng = buHuo(fengPaiXu_x(1));
      Conj = conj(feng(1));  
      path_1 = rx;
      path_1(1,1:end-shift_x(1)) = path_1(1,shift_x(1)+1:end);%����
      merge = path_1 * Conj(1);
    end
    %%%%����ǰ����������
    receive2 = sign(merge);
    BitErrorRake = length(find(receive2 ~= spreadDate_m)); 
    BitErrorRake_rate(snr) = BitErrorRake/(len_spread); 
end
%���Ӻ�۲�һ����ط�
y = step(xcor,merge(1:L_m)',m2'); %computes cross-correlation of x1 and x2
y2 = abs(y);
figure(4), subplot(212);plot(y2); title('Correlated output')

    