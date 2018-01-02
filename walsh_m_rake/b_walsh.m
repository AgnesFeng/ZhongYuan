clear all;
close all;
%----------------------1��walsh����---------------------
wal2 = [ 1 1; 1 -1 ];
wal4 = [wal2 wal2; wal2 wal2*(-1)];  
wal8 = [wal4 wal4; wal4 wal4*(-1)];  %8*8
wal16 = [wal8 wal8; wal8 wal8*(-1)]; %16*16
wal32 = [wal16 wal16; wal16 wal16*(-1)]; %32*32
p_0 = wal32(2,:);
p_1 = wal32(4,:);
KP = 32; %��Ƶ����
Tm_sample = 4;%������Խ�󣬽��ͼԽ��׼
pp0 = SAM_d(p_0,Tm_sample,length(p_0));%�����ڣ��г���ʱ��
pp1 = SAM_d(p_1,Tm_sample,length(p_1));%�����ڣ��г���ʱ��

%%
% %---------------------3���ź����У���Ƶ---------------------
Tlen = 100000; 
number_zhen = 1;
v = 225000;
Tb = 1/v; %��Ԫ����ʱ��  
Tc = Tb/KP;
dt = Tc/Tm_sample;
T_sample = Tm_sample * KP;  %m���еĲ�������Ϊ4
L_s_initial1 = Tlen * T_sample;
spreadData = ones(1,L_s_initial1);

temp = 15;%�������Ϊ����
EbNo = (1:1:20)-temp;
x = 1;
for snr = 1:1 :20
    BitErrorManyZhen_norake = 0;
    for number_s_initial = 1 : number_zhen
        s_initial = randsrc( 1, Tlen );
        %s_initial = [1 -1 1 -1 ];
        s_initial1 =  SAM_d(s_initial,T_sample,length(s_initial));
        %--------------------��Ƶ
        for i = 1: T_sample : L_s_initial1-T_sample+1
             if(s_initial1(i) == -1)
                 spreadData(i : i+T_sample-1) = pp0;%��Ϊ-1����Ԫ��Ƶ
             else
                 spreadData(i : i+T_sample-1) = pp1;%��Ϊ 1����Ԫ��Ƶ
             end
        end
    %---------------------�Ӷྶ�ŵ�--------------------------------
        channel_type=2;
        %-------------------�ŵ�����ѡ��-----------
        ts =dt;
        if(channel_type==1)
            fd = 25;
            k = 10^(12/10);
            chan = ricianchan(ts,fd,k);  
        elseif(channel_type==2)
            fd = 25;           %��������Ƶ�� 
            %fd = 0;
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
        c_out_spread = filter(chan,spreadData);
         c_out_origin = filter(chan,s_initial);
        %---------------------�Ӹ�˹����
        out = awgn(c_out_spread,snr-temp,'measured');%�����㷨���������EbNo�Ĵ�С
         dataNoDispose = awgn(c_out_origin,snr-temp,'measured');
        xcor = dsp.Crosscorrelator;
        y4 = step(xcor,out',pp1');
        y4 = abs(y4);
        y5 = step(xcor,out',pp0');
        y5 = abs(y5);
%--------------------�о�
% receiveJudgeNoRake = ones(1,Tlen);
%     tt = 1;
%     for i = 1: T_sample : L_s_initial1-T_sample+1 
%          if(mean(y4(i : i+T_sample-1)) > mean(y5(i : i+T_sample-1)))
%              receiveJudgeNoRake(tt) = 1;%��Ϊ-1����Ԫ��Ƶ
%          else
%              receiveJudgeNoRake(tt) = -1;
%          end
%          tt = tt+1;
%     end
% BitErrorNoRake11 = length(find(receiveJudgeNoRake ~= s_initial));%����ķ�ʽ
         %--------------------�о�
         sh = 25;
         receiveJudgeNoRake = ones(1,Tlen);
         tt = 1;
         for i = T_sample-sh : T_sample : L_s_initial1-T_sample-sh+1+T_sample 
             if(sum(y4(i : i+2*sh)) > sum(y5(i : i+2*sh)))
                 receiveJudgeNoRake(tt) = 1;%��Ϊ-1����Ԫ��Ƶ
             else
                 receiveJudgeNoRake(tt) = -1;
             end
             tt = tt+1;
         end
         BitErrorNoRake = length(find(receiveJudgeNoRake ~= s_initial));
         BitErrorManyZhen_norake = BitErrorManyZhen_norake + BitErrorNoRake;
    end
    BitErrorAll_norake(x) = BitErrorManyZhen_norake/(Tlen * number_zhen);
    com = sign(real(dataNoDispose));
    BitErrorNoDispose(x) = length(find(com ~= s_initial))/Tlen;
    x = x+1;
end
% 
% figure
% subplot(211);plot(abs(y4));title('����1��walsh');
% subplot(212);plot(abs(y5));title('����-1��walsh');
% 
% xcor = dsp.Crosscorrelator;
% y1 = step(xcor,pp0',pp0');   
% y1 = abs(y1);
% figure,subplot(311);plot(y1);title('��ȡһ�е��������');
% y2 = step(xcor,pp1',pp1'); 
% y2 = abs(y2);
% subplot(312);plot(y2);title('��ȡ��һ�е��������');
% y3 = step(xcor,pp1',pp0');          
% subplot(313);plot(y3);title('��ͬ��Ļ������');
% y3 = abs(y3);

figure
semilogy(EbNo,BitErrorAll_norake,'-b');grid;hold on;
semilogy(EbNo,BitErrorNoDispose,'-r');
legend('��ɫ����Ƶ���о�','��ɫ����������');
xlabel('�����'),ylabel('������ ȡlog'),title('�Ӹ�˹����������˹');

