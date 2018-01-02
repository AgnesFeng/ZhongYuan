clear all;        %������һ�����ڵ�m���л�ͼЧ������
close all;
clc;
%%
%----------------------1��walsh����(����)---------------------
% wal2 = [ 1 1; 1 -1 ];
% wal4 = [wal2 wal2; wal2 wal2*(-1)];  
% wal8 = [wal4 wal4; wal4 wal4*(-1)];  %8*8
% wal16 = [wal8 wal8; wal8 wal8*(-1)]; %16*16
% p = wal16(5,:);
% KP = length(p); %��Ƶ����
Tm_sample = 4;%������Խ�󣬽��ͼԽ��׼
% pp = SAM_d(p,Tm_sample,length(p));%�����ڣ��г���ʱ��

%----------------------1����Ƶ��m����---------------------
m_spread = m_sequence([0,1,0,0,1]);%�ܵõ�31λ��m����
N_ms = length(m_spread);%m��:31
ms1 =  SAM_m(m_spread,Tm_sample,length(m_spread)); %������ʱ��
ms2 = 1-2*ms1;  %��˫����
L_ms = length(ms2);
KP = N_ms;
%%%%%%%��ͼ
mms = repmat(ms2,1,3); %3�����ڵ�m���У����������
rt1=conv(mms,ms2(end:-1:1))/(L_ms);%ȡ�м��һ����
index = L_ms+Tm_sample;
figure(11),subplot(211);plot(rt1(index:end-index));title('������m���е������');

xcor = dsp.Crosscorrelator;
y = step(xcor,ms2',ms2'); 
figure(11),subplot(212);plot(y);title('�������m���е������');
%%
%%--------------------------2��ǰ׺m����----------------
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
figure(1);subplot(311);plot(rt1);title('������m���е�����غ���conv����ֵ������ֵ��m���г���2����');
figure(1);subplot(312);plot(rt1(index:end-index));title('������m���е�����غ���conv,��ֵ��������3��');
xcor = dsp.Crosscorrelator;
y = step(xcor,m2',m2'); 
figure(1);subplot(313);plot(y);title('������m���е�����غ��������������β�����');
%%
%---------------------3���ź����У���Ƶ,��ǰ������---------------------
Tlen = 5000; 
v = 225000;
Tb = 1/v; %��Ԫ����ʱ��  
Tc = Tb/KP;
dt = Tc/Tm_sample;
s_initial = randsrc( 1, Tlen );
T_sample = Tm_sample * KP;  %m���еĲ�������Ϊ4
s_initial1 =  SAM_d(s_initial,T_sample,length(s_initial));
pt = repmat(ms2,1,Tlen);%��Ƶ������
spreadData = pt.*s_initial1;%��Ƶ
spreadData_m = ones(1,length(m2)+length(s_initial1));%��ǰ������
spreadData_m(1:length(m2)) = m2;
spreadData_m(length(m2)+1:end) = spreadData;
len_spread = length(spreadData_m);
%%
%�ҵ���λ����ʱ���±�
x1 = spreadData_m(1:length(m2));
x2 = m2;
xcor = dsp.Crosscorrelator;
y = step(xcor,x2',x1'); %computes cross-correlation of x1 and x2
figure(2), subplot(211);plot(y); title('Correlated output');
aim = find(max(y) == y);

%%
%--------------------������������-------------------
txx = mo_bpsk(spreadData_m,v,Tm_sample);
%plot(txx(1:100));
%%
%channel_type = input('�������ŵ����ͣ�'); 
channel_type=3;
%-------------------�ŵ�����ѡ��-----------------
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
c_out = filter(chan,txx);
% figure(11),
% subplot(311);plot(s_initial1);axis([0,500,-1.3,1.3]);
% subplot(312),plot(spreadData);axis([0,800,-1.3,1.3]);
% subplot(313),plot(sign(real(c_out)));axis([0,800,-1.3,1.3]);
%--------------------��˹����-----------------------
EbNo = 1:1:12;
for snr = 1:length(EbNo)
   y_snr = awgn(c_out,snr,'measured');

%%
    %--------------------ʡ���˽��-------------------
    rx = y_snr;

%%
    %------------����ǰ�����к�������Ƶ��m����
    m_modu = mo_bpsk(m2,v,Tm_sample);    %����ǰ������
    m_despread = mo_bpsk(pt,v,Tm_sample);  %����������Ƶ��m���У�ȥ������
    %plot(m_modu);
    %------------��ط岶��
    xcor = dsp.Crosscorrelator;
    buHuo = step(xcor,rx',m_modu'); %computes cross-correlation of x1 and x2
    buHuoQuMo = abs(buHuo);
    figure(2), subplot(212);plot(buHuoQuMo); title('����ص���ط�');axis([0 2000 0 600]);

    y_max = max(buHuoQuMo);
    id = find(y_max == buHuoQuMo); %�ҵ����ֵ�������ҿ���,��1��Ԫ����ʱ����4΢�����ڣ�
    L_catch_window = Tm_sample * KP*2 ;
    feng_x = ones(1,10);
    a = 1;
    b1 = id - L_catch_window/2;
    b2 = id + L_catch_window/2;
    for x =b1 : 1 : b2
        if buHuoQuMo(x) >= y_max/4 && buHuoQuMo(x-1) < buHuoQuMo(x) && buHuoQuMo(x) > buHuoQuMo(x+1) && buHuoQuMo(x-2) < buHuoQuMo(x) && buHuoQuMo(x) > buHuoQuMo(x+2)
            feng_x(a) = x;
            a = a+1;
        end
    end
    fengPaiXu_x = feng_x;
    for j = 1:a-1   %�����ֵ�ĺ����������ҳ���������ֵ�ĺ�����
        for i = 1:a-1-j
            if buHuoQuMo(fengPaiXu_x(i+1)) > buHuoQuMo(fengPaiXu_x(i))
              temp = fengPaiXu_x(i);
              fengPaiXu_x(i) = fengPaiXu_x(i+1);
              fengPaiXu_x(i+1) = temp;
            end
        end
    end
%%
    %--------------------������λ
%     jiaoZheng = step(xcor,rx(L_m+1:end)',m_despread(1:T_sample)');%��Ҫһ�����ڵ�m����
%     figure(3),plot(abs(jiaoZheng)); title('���ݲ��ֵ���ط�');%axis([0 2000 0 600]);
    
%%
    %----------------�ж�����������
    %%%%%%%%��׽��3��������
    if(feng_x(3)~=1)              
      shift_x = ones(1,3);  
      feng = ones(1,3);
      Conj = ones(1,3);
      for i = 1:3
        shift_x(i) = fengPaiXu_x(i)-aim;      %��¼λ��
        if(shift_x(i)<0)
           shift_x(i)=abs(fengPaiXu_x(i)-aim);
        end
        feng(i) = buHuo(fengPaiXu_x(i));
        Conj(i) = conj(feng(i));
      end      
      path_1 = rx;
      path_2 = rx;
      path_3 = rx;
%       path_1(1,shift_x(1)+1:end) = path_1(1,1:end-shift_x(1));%���ƣ��õ���ͬ��ʱ������
%       path_2(1,shift_x(2)+1:end) = path_2(1,1:end-shift_x(2));
%       path_3(1,shift_x(3)+1:end) = path_3(1,1:end-shift_x(3));
      path_1(1,1:end-shift_x(1)) = path_1(1,shift_x(1)+1:end);%����
      path_2(1,1:end-shift_x(2)) = path_2(1,shift_x(2)+1:end);
      path_3(1,1:end-shift_x(3)) = path_3(1,shift_x(3)+1:end);
      merge = path_1 * conj(1) + path_2 * Conj(2) + path_3 * Conj(3);
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
%     path_1(1,shift_x(1)+1:end) = path_1(1,1:end-shift_x(1)); %���ƣ��õ���ͬ��ʱ������
%     path_2(1,shift_x(2)+1:end) = path_2(1,1:end-shift_x(2));
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
  %%
    %%������
    data_out = merge(L_m+1:end); %ȥ��ǰ׺����
    despreadData = m_despread.*data_out;%���ƺ����Ƶm���н���
    receive = ones(1,Tlen);
    k=1;
    for i = 1:T_sample:length(despreadData)-T_sample+1
        receive(k) = 0;
        for j = 0:T_sample-1
          receive(k) = receive(k) + despreadData(i+j);
        end
        k = k+1;
    end
    receiveJudge = sign(real(receive));%��ʵ���о�
    if(length(find(receiveJudge ~= s_initial))<Tlen/2)
       BitErrorRake = length(find(receiveJudge ~= s_initial)); 
    else
       BitErrorRake = length(find(receiveJudge == s_initial)); 
    end
    BitErrorRake_rate(snr) = BitErrorRake/(Tlen); 
    %%�Ա�ֱ��������
    despreadData_no = m_despread.* rx(L_m+1:end);
    receive2 = ones(1,Tlen);
    k=1;
    for i = 1:T_sample:length(despreadData_no)-T_sample+1
        receive2(k) = 0;
        for j = 0:T_sample-1
          receive2(k) = receive2(k) + despreadData_no(i+j);
        end
        k = k+1;
    end 
    receiveJudge_no = sign(real(receive2));
    number_norake_different = length(find(receiveJudge_no ~= s_initial));
    if(number_norake_different<Tlen/2)
        BitErrorCompare = number_norake_different; 
    else
        BitErrorCompare = length(find(receiveJudge_no == s_initial));
    end
    BitErrorRakeCompare(snr) = BitErrorCompare/(Tlen); 
    
end


figure(111)
subplot(121),semilogy(EbNo,BitErrorRakeCompare,'b');hold on;title('������Ƶ��������ͬ��');
subplot(122),semilogy(EbNo,BitErrorRake_rate,'r');hold on;title('rake���ջ�');



% %�۲�һ�µ��Ӻ����ط�
% %merge(1,1:end-s(1)) = merge(1,s(1)+1:end);
% y = step(xcor,merge',m2'); %computes cross-correlation of x1 and x2
% y2 = abs(y);
% figure(4);plot(y2); title('Correlated output');
% axis([0 2000 0 50e10000]);



    