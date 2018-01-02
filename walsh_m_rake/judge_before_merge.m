clear all;        %������һ�����ڵ�m���л�ͼЧ������
close all;
clc;
%%
number_zhen = 2;%֡��
Tlen = 1000; %ÿһ֡����Ԫ��
%channel_type = input('�������ŵ����ͣ�'); 
channel_type=1;
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
N_m = length(m);%m��:225
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
jian = 10;%�������Ϊ����
EbNo = (1:1:20)-jian;
for snr = 1:1:20   
    BitErrorManyZhen_rake = 0;
    BitErrorManyZhen_rake_a = 0;
    BitErrorManyZhen_rake_b = 0;
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
            %fd = 0;           %��������Ƶ�� 
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
        buHuo = step(xcor,rx(1:L_m)',m2'); %computes cross-correlation of x1 and x2
        buHuoQuMo = abs(buHuo);
        y_max = max(buHuoQuMo);
%         figure(2), subplot(211);plot(y); title('Correlated output');
%         figure(2), subplot(212);plot(buHuoQuMo); title('����ص���ط�');axis([0 2000 0 900]);
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
            rake_a = path_1 * conj(1);
            rake_b = path_2 * conj(2);
            rake_c = path_3 * conj(3);
            merge = rake_a + rake_b + rake_c;
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
            rake_a = path_1 * conj(1);
            rake_b = path_2 * conj(2);
            merge = rake_a + rake_b;
       end
      %%%%%%%%%%%%%%%%%%%��׽��1��
       if(feng_x(1)~=1 && feng_x(3)==1 && feng_x(2)==1) 
            shift_x = fengPaiXu_x(1)-aim; 
            feng = buHuo(fengPaiXu_x(1));
            Conj = conj(feng(1));  
            path_1 = rx;
            path_1(1,1:end-shift_x(1)) = path_1(1,shift_x(1)+1:end);%����
            rake_a = path_1 * conj(1);
            merge = rake_a;
       end
      %%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%�����
        data_out = merge(L_m+1:end); %ȥ��ǰ׺���֣��ϲ���
        data_out1 = rake_a(L_m+1:end); %ȥ��ǰ׺����,��һ·��
        %data_out2 = rake_b(L_m+1:end); %ȥ��ǰ׺����,�ڶ�·��
        %data_out3 = rake_c(L_m+1:end); %ȥ��ǰ׺����,����·��
        xcor = dsp.Crosscorrelator;
        y3 = step(xcor,data_out',pp1'); %�ϲ��������
        y1Rake = abs(y3);
        y4 = step(xcor,data_out',pp0');
        y0Rake = abs(y4);
        y1Rake_a = abs(step(xcor,data_out1',pp1'));%��һ·�������
        y0Rake_a = abs(step(xcor,data_out1',pp0'));
%         y1Rake_b = abs(step(xcor,data_out2',pp1'));%�ڶ�·�������
%         y0Rake_b = abs(step(xcor,data_out2',pp0'));
        xcor = dsp.Crosscorrelator;
        data_out_norake = rx(L_m+1:end);          %����rake��
        y5 = step(xcor,data_out_norake',pp1');
        y1NoRake = abs(y5);
        y6 = step(xcor,data_out_norake',pp0');
        y0NoRake = abs(y6);
       %%
        %��ÿ����Ԫ�ڵ��о���
        receiveJudge = ones(1,Tlen);              %%%%%%%%%%%%%%%%%%%%%%%%%
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
        %--��һ· %%%%%%%%%%%%%%%%%%%%%%%%%
        receiveJudge_a = ones(1,Tlen);           
        tt = 1;
        for i = T_sample-sh : T_sample : L_s_initial1 - sh + 1
             if(sum(y1Rake_a(i : i+2*sh)) > sum(y0Rake_a(i : i+2*sh)))
                 receiveJudge_a(tt) = 1;
             else
                 receiveJudge_a(tt) = -1;
             end
             tt = tt+1;
        end
        BitErrorRake_a = length(find(receiveJudge_a ~= s_initial));
        BitErrorManyZhen_rake_a = BitErrorManyZhen_rake_a + BitErrorRake_a;
        %--�ڶ�· %%%%%%%%%%%%%%%%%%%%%%%%%
%        if(feng_x(2)~=1 && feng_x(3)==1)
%             receiveJudge_b = ones(1,Tlen);           
%             tt = 1;
%             for i = T_sample-sh : T_sample : L_s_initial1 - sh + 1
%                  if(sum(y1Rake_b(i : i+2*sh)) > sum(y0Rake_b(i : i+2*sh)))
%                      receiveJudge_b(tt) = 1;
%                  else
%                      receiveJudge_b(tt) = -1;
%                  end
%                  tt = tt+1;
%             end
%             BitErrorRake_b = length(find(receiveJudge_b ~= s_initial));
%             BitErrorManyZhen_rake_b = BitErrorManyZhen_rake_b + BitErrorRake_a;
%        end
        %--����· %%%%%%%%%%%%%%%%%%%%%%%%%
        
        %�ԱȲ�ͬ��������                          %%%%%%%%%%%%%%%%%%%%%%%%%
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
    ErrorRate_a(xx) = BitErrorManyZhen_rake_a/(Tlen * number_zhen);
   % ErrorRate_b(xx) = BitErrorManyZhen_rake_b/(Tlen * number_zhen);
    
    xx = xx + 1;
end
figure(2), subplot(211);plot(y); title('Correlated output');
figure(2), subplot(212);plot(buHuoQuMo); title('����ص���ط�');axis([0 2000 0 900]);

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
semilogy(EbNo,ErrorRateNoRake,'m*-');hold on;grid;
semilogy(EbNo,ErrorRate_a,'bo-');hold on;grid;
semilogy(EbNo,ErrorRate,'r*-');hold on;title('rake���ջ�');
legend('����Ƶ�����о�','rake��һ·','rake��·�ϲ�');
% semilogy(EbNo,BitErrorRake_b,'go-');hold on;
% semilogy(EbNo,BitErrorRake_c,'mo-');hold on;  %������·��

 

    