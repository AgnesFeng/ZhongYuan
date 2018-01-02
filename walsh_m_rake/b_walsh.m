clear all;
close all;
%----------------------1、walsh序列---------------------
wal2 = [ 1 1; 1 -1 ];
wal4 = [wal2 wal2; wal2 wal2*(-1)];  
wal8 = [wal4 wal4; wal4 wal4*(-1)];  %8*8
wal16 = [wal8 wal8; wal8 wal8*(-1)]; %16*16
wal32 = [wal16 wal16; wal16 wal16*(-1)]; %32*32
p_0 = wal32(2,:);
p_1 = wal32(4,:);
KP = 32; %扩频因子
Tm_sample = 4;%采样率越大，结果图越精准
pp0 = SAM_d(p_0,Tm_sample,length(p_0));%多周期，有持续时间
pp1 = SAM_d(p_1,Tm_sample,length(p_1));%多周期，有持续时间

%%
% %---------------------3、信号序列，扩频---------------------
Tlen = 100000; 
number_zhen = 1;
v = 225000;
Tb = 1/v; %码元持续时间  
Tc = Tb/KP;
dt = Tc/Tm_sample;
T_sample = Tm_sample * KP;  %m序列的采样点数为4
L_s_initial1 = Tlen * T_sample;
spreadData = ones(1,L_s_initial1);

temp = 15;%让信噪比为负的
EbNo = (1:1:20)-temp;
x = 1;
for snr = 1:1 :20
    BitErrorManyZhen_norake = 0;
    for number_s_initial = 1 : number_zhen
        s_initial = randsrc( 1, Tlen );
        %s_initial = [1 -1 1 -1 ];
        s_initial1 =  SAM_d(s_initial,T_sample,length(s_initial));
        %--------------------扩频
        for i = 1: T_sample : L_s_initial1-T_sample+1
             if(s_initial1(i) == -1)
                 spreadData(i : i+T_sample-1) = pp0;%对为-1的码元扩频
             else
                 spreadData(i : i+T_sample-1) = pp1;%对为 1的码元扩频
             end
        end
    %---------------------加多径信道--------------------------------
        channel_type=2;
        %-------------------信道类型选择-----------
        ts =dt;
        if(channel_type==1)
            fd = 25;
            k = 10^(12/10);
            chan = ricianchan(ts,fd,k);  
        elseif(channel_type==2)
            fd = 25;           %最大多普勒频移 
            %fd = 0;
            k = 10^(6/10);
            tau = [0 0.000001 0.000002];
            pdb = [0,-9,-12];
            chan = ricianchan(ts,fd,k,tau,pdb); %多径，要求每一径都是瑞丽衰落。ts：输入信号的采样时间（s）；fd：最大多普勒频移；pdb：平均路径增益（dB）
        elseif(channel_type==3)
            fd = 10;     %最大多普勒频移 
            tau = [0 0.000001 0.000002 0.000003];
            pdb = [0,-3,-6,-9];
            chan = rayleighchan(ts,fd,tau,pdb);
        else                                    % (channel_type~=3 &&channel_type~=1 &&channel_type~=2)
            error('Error! channel type should be one of 1 2 3!');
        end 
        c_out_spread = filter(chan,spreadData);
         c_out_origin = filter(chan,s_initial);
        %---------------------加高斯噪声
        out = awgn(c_out_spread,snr-temp,'measured');%参与算法的信噪比是EbNo的大小
         dataNoDispose = awgn(c_out_origin,snr-temp,'measured');
        xcor = dsp.Crosscorrelator;
        y4 = step(xcor,out',pp1');
        y4 = abs(y4);
        y5 = step(xcor,out',pp0');
        y5 = abs(y5);
%--------------------判决
% receiveJudgeNoRake = ones(1,Tlen);
%     tt = 1;
%     for i = 1: T_sample : L_s_initial1-T_sample+1 
%          if(mean(y4(i : i+T_sample-1)) > mean(y5(i : i+T_sample-1)))
%              receiveJudgeNoRake(tt) = 1;%对为-1的码元扩频
%          else
%              receiveJudgeNoRake(tt) = -1;
%          end
%          tt = tt+1;
%     end
% BitErrorNoRake11 = length(find(receiveJudgeNoRake ~= s_initial));%错误的方式
         %--------------------判决
         sh = 25;
         receiveJudgeNoRake = ones(1,Tlen);
         tt = 1;
         for i = T_sample-sh : T_sample : L_s_initial1-T_sample-sh+1+T_sample 
             if(sum(y4(i : i+2*sh)) > sum(y5(i : i+2*sh)))
                 receiveJudgeNoRake(tt) = 1;%对为-1的码元扩频
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
% subplot(211);plot(abs(y4));title('用扩1的walsh');
% subplot(212);plot(abs(y5));title('用扩-1的walsh');
% 
% xcor = dsp.Crosscorrelator;
% y1 = step(xcor,pp0',pp0');   
% y1 = abs(y1);
% figure,subplot(311);plot(y1);title('任取一行的自相关性');
% y2 = step(xcor,pp1',pp1'); 
% y2 = abs(y2);
% subplot(312);plot(y2);title('任取另一行的自相关性');
% y3 = step(xcor,pp1',pp0');          
% subplot(313);plot(y3);title('不同码的互相关性');
% y3 = abs(y3);

figure
semilogy(EbNo,BitErrorAll_norake,'-b');grid;hold on;
semilogy(EbNo,BitErrorNoDispose,'-r');
legend('蓝色：扩频软判决','红色：不做处理');
xlabel('信噪比'),ylabel('误码率 取log'),title('加高斯，过三径莱斯');

