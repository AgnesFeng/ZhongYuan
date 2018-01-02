clear all;        %仅生成一个周期的m序列画图效果不好
close all;
clc;
%%
Tlen = 1000; 
number_zhen = 10;
Tm_sample = 4;%采样率越大，结果图越精准
%----------------------1、扩频的m序列---------------------
m_spread = m_sequence([0,1,0,0,1]);%能得到31位的m序列
N_ms = length(m_spread);%m长:31
ms1 =  SAM_m(m_spread,Tm_sample,length(m_spread)); %给持续时间
ms2 = 1-2*ms1;  %变双极性
L_ms = length(ms2);
KP = N_ms;
%%
%%--------------------------2、前缀m序列----------------
m = m_sequence([0,0,0,1,1,1,0,1]); %得到一个周期的m序列 (各寄存器的初始状态，从左到右为C1~Cn，c0默认为1) 
N_m = length(m);%m长:225
m1 =  SAM_m(m,Tm_sample,length(m)); %给持续时间
m2 = 1-2*m1;  %变双极性
L_m = length(m2);
%%
%---------------------3、信号序列，扩频,加前导序列---------------------
v = 225000;
Tb = 1/v; %码元持续时间  
Tc = Tb/KP;
dt = Tc/Tm_sample;
T_sample = Tm_sample * KP;  %m序列的采样点数为4
pt = repmat(ms2,1,Tlen);%扩频的序列
L_s_initial1 = Tlen * T_sample;
spreadData_m = ones(1,L_m + L_s_initial1);%加前导序列
spreadData_m(1:L_m) = m2;
%%
%----------------------调制前导序列和用于扩频的m序列--------------------
%m_modu = mo_bpsk(m2,v,Tm_sample);    %调制前导序列
m_modu = m2;
%m_despread = mo_bpsk(pt,v,Tm_sample);  %调制用于扩频的m序列，去软解扩
m_despread = pt;
%plot(m_modu);
%%%%%%%%%%%%%开始重复
xx= 1;
jian = 15;%让信噪比为负的
EbNo = (1:1:20)-jian;
for snr = 1:1:20
    BitErrorMuchZhen_rake = 0;
    BitErrorMuchZhen_norake = 0;
    for number_s_initial = 1:number_zhen
        s_initial = randsrc( 1, Tlen );
        s_initial1 =  SAM_d(s_initial,T_sample,Tlen);
        spreadData = pt.*s_initial1;%扩频
        spreadData_m(length(m2)+1:end) = spreadData;
       %%
        %找到相位对齐时的下标
        x1 = spreadData_m(1:length(m2));
        x2 = m2;
        xcor = dsp.Crosscorrelator;
        y = step(xcor,x2',x1'); %computes cross-correlation of x1 and x2
        % figure(2), subplot(211);plot(y); title('Correlated output');
        aim = find(max(y) == y);

        %%
        %--------------------调制整个序列-------------------
        %txx = mo_bpsk(spreadData_m,v,Tm_sample);
        txx = spreadData_m;
        %plot(txx(1:100));
        %%
        %channel_type = input('请输入信道类型：'); 
        channel_type=1;
        %-------------------信道类型选择-----------------
        ts =dt;
        if(channel_type==1)
            fd = 25;
            k = 10^(12/10);
            chan = ricianchan(ts,fd,k);  
        elseif(channel_type==2)
            fd = 25;           %最大多普勒频移 
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
        c_out = filter(chan,txx);
        c_outMerge = 0.5* c_out + 0.3 * circshift(c_out,[0,29]) + 0.2 * circshift(c_out,[0,58]);
        %--------------------高斯噪声-----------------------
        y_snr = awgn(c_outMerge,snr-jian,'measured');
    %%
        %--------------------省略了解调-------------------
        rx = y_snr;

    %%
        %------------相关峰捕获
        xcor = dsp.Crosscorrelator;
        buHuo = step(xcor,rx',m_modu'); %computes cross-correlation of x1 and x2
        buHuoQuMo = abs(buHuo);
       %figure(2),subolot(212);plot(buHuoQuMo); title('软相关的相关峰');axis([0 2000 0 600]);

        y_max = max(buHuoQuMo);
        id = find(y_max == buHuoQuMo); %找到最大值处，左右开窗,各1码元（延时设在4微秒以内）
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
        for j = 1:a-1   %把最大值的横坐标排序，找出最大的三个值的横坐标
            for i = 1:a-1-j
                if buHuoQuMo(fengPaiXu_x(i+1)) > buHuoQuMo(fengPaiXu_x(i))
                  temp = fengPaiXu_x(i);
                  fengPaiXu_x(i) = fengPaiXu_x(i+1);
                  fengPaiXu_x(i+1) = temp;
                end
            end
        end

    %%
    %------------------得到矫正相位用的峰值共轭-------------
    %jiaoZheng = step(xcor,rx(L_m+1:end)',m_modu');
    %%
        %----------------判断搜索到几径
        %%%%%%%%捕捉到3个或以上
        if(feng_x(3)~=1)              
          shift_x = ones(1,3);  
          feng = ones(1,3);
          Conj = ones(1,3);
          for i = 1:3
            shift_x(i) = fengPaiXu_x(i)-aim;      %记录位移
            if(shift_x(i)<0)
               shift_x(i)=abs(fengPaiXu_x(i)-aim);
            end
            feng(i) = buHuo(fengPaiXu_x(i));
            Conj(i) = conj(feng(i));
          end      
          path_1 = rx;
          path_2 = rx;
          path_3 = rx;
          path_1(1,1:end-shift_x(1)) = path_1(1,shift_x(1)+1:end);%左移
          path_2(1,1:end-shift_x(2)) = path_2(1,shift_x(2)+1:end);
          path_3(1,1:end-shift_x(3)) = path_3(1,shift_x(3)+1:end);
          rake_a = path_1 * Conj(1);
          rake_b = path_2 * Conj(2);
          rake_c = path_3 * Conj(3);
          merge = rake_a + rake_b + rake_c;
        end
      %%%%%%%%%%%%%%%%%%捕捉到2个
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
             path_1(1,1:end-shift_x(1)) = path_1(1,shift_x(1)+1:end);%左移
             path_2(1,1:end-shift_x(2)) = path_2(1,shift_x(2)+1:end);
             rake_a = path_1 * Conj(1);
             rake_b = path_2 * Conj(2);
             merge = rake_a + rake_b;             
        end
      %%%%%%%%%%%%%%%%%%%捕捉到1个
      if(feng_x(1)~=1 && feng_x(3)==1 && feng_x(2)==1) 
        shift_x = fengPaiXu_x(1)-aim; 
        feng = buHuo(fengPaiXu_x(1));
        Conj = conj(feng(1));  
        path_1 = rx;
        path_1(1,1:end-shift_x(1)) = path_1(1,shift_x(1)+1:end);%左移
        merge = path_1 * Conj(1);
      end
      %%
        %%软解扩
        data_out = merge(L_m+1:end); %去掉前缀部分
        despreadData = m_despread.* data_out;%调制后的扩频m序列解扩
        receive_merge = ones(1,Tlen);
        k=1;
        for i = 1:T_sample:length(despreadData)-T_sample+1
            receive_merge(k) = 0;
            for j = 0:T_sample-1
              receive_merge(k) = receive_merge(k) + despreadData(i+j);
            end
            k = k+1;
        end
        receiveJudge = sign(real(receive_merge));%对实部判决
       %%
        %%%判决每一路的
%         data_out = merge(L_m+1:end); %去掉前缀部分
%         despreadData = m_despread.* data_out;%调制后的扩频m序列解扩
%         receive_merge = ones(1,Tlen);
%         k=1;
%         for i = 1:T_sample:length(despreadData)-T_sample+1
%             receive_merge(k) = 0;
%             for j = 0:T_sample-1
%               receive_merge(k) = receive_merge(k) + despreadData(i+j);
%             end
%             k = k+1;
%         end
%         receiveJudge = sign(real(receive_merge));%对实部判决
       %%
        %%%%%%%%%%%%%收集每一帧
        number_different = length(find(receiveJudge ~= s_initial));
        if(number_different<Tlen/2)
            BitErrorRake = number_different; 
        else
            BitErrorRake = length(find(receiveJudge == s_initial));
        end
        BitErrorMuchZhen_rake = BitErrorMuchZhen_rake + BitErrorRake;
        %%对比直接解扩
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
        BitErrorMuchZhen_norake = BitErrorMuchZhen_norake + BitErrorCompare;
    end
    BitErrorRake_rate(xx) = BitErrorMuchZhen_rake/(Tlen * number_zhen); 
    BitErrorRakeCompare(xx) = BitErrorMuchZhen_norake/(Tlen * number_zhen); 
    xx = xx + 1;
end
figure(11)
semilogy(EbNo,BitErrorRakeCompare,'mo-');hold on;
semilogy(EbNo,BitErrorRake_rate,'c*-');hold on;
% semilogy(EbNo,BitErrorRake_second,'-y');hold on;
% semilogy(EbNo,BitErrorRake_third,'-y');hold on;
title('m序列扩频方式对比');
legend('单纯扩频解扩','rake三路合并');

    