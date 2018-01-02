clear all;        %������һ�����ڵ�m���л�ͼЧ������
close all;
clc;
%%
number_zhen = 1;%֡��
Tlen = 1000; %ÿһ֡����Ԫ��
%channel_type = input('�������ŵ����ͣ�'); 
channel_type=2;
v = 225000; %��Ԫ����
Tm_sample = 4;%������Խ�󣬽��ͼԽ��׼

%----------------------1��walsh����---------------------
wal2 = [ 1 1; 1 -1 ];
wal4 = [wal2 wal2; wal2 wal2*(-1)];  
wal8 = [wal4 wal4; wal4 wal4*(-1)];  %8*8
wal16 = [wal8 wal8; wal8 wal8*(-1)]; %16*16
wal32 = [wal16 wal16; wal16 wal16*(-1)]; %32*32
p_0 = wal32(5,:);
p_1 = wal32(11,:);
KP = 32; %��Ƶ����
pp0 = SAM_d(p_0,Tm_sample,length(p_0));%�����ڣ��г���ʱ��
pp1 = SAM_d(p_1,Tm_sample,length(p_1));%�����ڣ��г���ʱ��
%%
%%--------------------------2��m����----------------
m = m_sequence([0,0,0,1,1,1,0,1]); %�õ�һ�����ڵ�m���� (���Ĵ����ĳ�ʼ״̬��������ΪC1~Cn��c0Ĭ��Ϊ1) 
N_m = length(m);%m��:225,8jie
m1 =  SAM_m(m,Tm_sample,length(m)); %������ʱ��
m2 = 1-2*m1;  %��˫����
L_m = length(m2);
%%
% %---------------------3���ź����У���Ƶ,��ǰ������---------------------
Tb = 1/v; %��Ԫ����ʱ��  
Tc = Tb/KP;
dt = Tc/Tm_sample;
T_sample = Tm_sample * KP;  %m���еĲ�������Ϊ4
L_s_initial1 = Tlen * T_sample;
spreadData = ones(1,L_s_initial1);
spreadData_m = ones(1,length(m2) + L_s_initial1);%��ǰ������
spreadData_m(1:length(m2)) = m2;
xx= 1;
jian = 0;%�������Ϊ����
EbNo = (1:0.5:6)-jian;
for snr = 1:0.5:6   
    BitErrorManyZhen_rake = 0;
    BitErrorManyZhen_norake = 0;
    for number_s_initial = 1 : number_zhen
        s_initial = randsrc( 1, Tlen );
        s_initial1 =  SAM_d(s_initial,T_sample,Tlen);
        for i = 1: T_sample : L_s_initial1-T_sample+1 %��Ƶ
             if(s_initial1(i) == -1)
                 spreadData(i : i+T_sample-1) = pp0;%��Ϊ-1����Ԫ��Ƶ
             else
                 spreadData(i : i+T_sample-1) = pp1;%��Ϊ 1����Ԫ��Ƶ
             end
        end
        spreadData_m(length(m2)+1:end) = spreadData;
        %len_spread = length(spreadData_m);
        %%
        %�ҵ���λ����ʱ���±�
        x1 = spreadData_m(1:length(m2));
        x2 = m2;
        xcor = dsp.Crosscorrelator;
        y = step(xcor,x2',x1'); %computes cross-correlation of x1 and x2
        aim = find(max(y) == y);

        %%
        %--------------------������������-------------------
        %txx = mo_bpsk(spreadData_m,v,Tm_sample);
        %plot(txx(1:100));
        txx = spreadData_m;
        %%
        %------------------------�ŵ�����--------------------
        ts =dt;
        if(channel_type==1)
            fd = 25;
            k = 10^(12/10);
            chan = ricianchan(ts,fd,k);  
        elseif(channel_type==2)
            %fd = 25;           %��������Ƶ�� 
            fd = 25;
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
        %--------------------��˹����-----------------------
        y_snr = awgn(c_out,snr-jian,'measured');

    %%
        %--------------------ʡ���˽��-------------------
        rx = y_snr;

    %%
%         %------------��ط岶��                           %%%%%%%%%%%%%%%%%%%%%%%%%
        xcor = dsp.Crosscorrelator;
        buHuo = step(xcor,rx',m2'); %computes cross-correlation of x1 and x2
        buHuoQuMo = abs(buHuo);
        y_max = max(buHuoQuMo);
        figure(2), subplot(211);plot(y); title('Correlated output');
        figure(2), subplot(212);plot(buHuoQuMo); title('����ص���ط�');axis([0 2000 0 900]);
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
          %merge = path_1 ;
        end
%       %%%%%%%%%%%%%%%%%%��׽��2��
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
         %merge = path_1 ;
       end
%       %%%%%%%%%%%%%%%%%%%��׽��1��
%       if(feng_x(1)~=1 && feng_x(3)==1 && feng_x(2)==1) 
%         shift_x = fengPaiXu_x(1)-aim; 
%         feng = buHuo(fengPaiXu_x(1));
%         Conj = conj(feng(1));  
%         path_1 = rx;
%         path_1(1,1:end-shift_x(1)) = path_1(1,shift_x(1)+1:end);%����
%         merge = path_1 * Conj(1);
%       end
      %%
        %%�����
        data_out = merge(L_m+1:end); %ȥ��ǰ׺����
%         xcor = dsp.Crosscorrelator;
%         y1 = step(xcor,spreadData',pp1');
%         y1Origin = abs(y1);
%         y0 = step(xcor,spreadData',pp0');
%         y0Origin = abs(y0);
        xcor = dsp.Crosscorrelator;
        y3 = step(xcor,data_out',pp1');
        y1Rake = abs(y3);
        y4 = step(xcor,data_out',pp0');
        y0Rake = abs(y4);
        xcor = dsp.Crosscorrelator;
        data_out_norake = rx(L_m+1:end);
        y5 = step(xcor,data_out_norake',pp1');
        y1NoRake = abs(y5);
        y6 = step(xcor,data_out_norake',pp0');
        y0NoRake = abs(y6);

        %��ÿ����Ԫ�ڵ����ֵ
        receiveJudge = ones(1,Tlen);            %%%%%%%%%%%%%%%%%%%%%%%%%
        tt = 1;
        sh = 20;  %���о�����Ϊ��2*sh
        for i = T_sample-sh : T_sample : L_s_initial1 - sh + 1
             if(sum(y1Rake(i : i+2*sh)) > sum(y0Rake(i : i+2*sh)))
                 receiveJudge(tt) = 1;
             else
                 receiveJudge(tt) = -1;
             end
             tt = tt+1;
        end
        BitErrorRake = length(find(receiveJudge ~= s_initial));
        BitErrorManyZhen_rake = BitErrorManyZhen_rake + BitErrorRake;
        %�ԱȲ�ͬ��������
        receiveJudgeNoRake = ones(1,Tlen);
        tt = 1;
        for i = T_sample-sh : T_sample : L_s_initial1 - sh + 1 
             if(sum(y1NoRake(i : i+2*sh)) > sum(y0NoRake(i : i+2*sh)))
                 receiveJudgeNoRake(tt) = 1;
             else
                 receiveJudgeNoRake(tt) = -1;
             end
             tt = tt+1;
        end
        BitErrorNoRake = length(find(receiveJudgeNoRake ~= s_initial));
        BitErrorManyZhen_norake = BitErrorManyZhen_norake + BitErrorNoRake;
    end
    ErrorRate(xx) = BitErrorManyZhen_rake/(Tlen * number_zhen);
    ErrorRateNoRake(xx) = BitErrorManyZhen_norake/(Tlen * number_zhen); 
    xx = xx + 1;
end
% figure(2), subplot(211);plot(y); title('Correlated output');
% figure(2), subplot(212);plot(buHuoQuMo); title('����ص���ط�');axis([0 2000 0 900]);

        xcor = dsp.Crosscorrelator;
        y1 = step(xcor,spreadData',pp1');
        y1Origin = abs(y1);
        y0 = step(xcor,spreadData',pp0');
        y0Origin = abs(y0);
figure
subplot(411);plot(y1Origin);title('δ���ŵ�������1��walsh');axis([0 1280 0 150]);
subplot(412);plot(y0Origin);title('δ���ŵ�������-1��walsh');axis([0 1280 0 150]);
subplot(413);plot(y1Rake);title('����������1��walsh���');axis([0 1280 0 200000]);
subplot(414);plot(y0Rake);title('����������-1��walsh���');axis([0 1280 0 200000]);
figure
subplot(411);plot(y1Origin);title('δ���ŵ�������1��walsh');axis([0 1280 0 150]);
subplot(412);plot(y0Origin);title('δ���ŵ�������-1��walsh');axis([0 1280 0 150]);
subplot(413);plot(y1NoRake);title('��ͬ������������1��walsh���');axis([0 1280 0 150]);
subplot(414);plot(y0NoRake);title('��ͬ��У��������-1��walsh���');axis([0 1280 0 150]);

figure(11)
semilogy(EbNo,ErrorRateNoRake,'-b');hold on;title('������Ƶ���о�');grid;
figure(22)
semilogy(EbNo,ErrorRate,'r');hold on;title('rake���ջ�ͬ������');grid;
    
 

    