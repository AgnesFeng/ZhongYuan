clear all;
close all;
clc;
%-----------------------1��ԭʼ�ź�-----------------------------------------
Tlen = 1000; %���ݳ���  �������������9�Ĺ�ϵ
s_initial = randsrc( 1, Tlen );
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%����
% R = 1/3;              % code rate
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
T_sample = 310;
s_initial1 =  SAM2(s_initial,T_sample,length(s_initial),Tb); 

%-----------------------2��m����-------------------------------------------
m = m_sequence([0,1,0,0,1]); %�õ�һ�����ڵ�m���� (���Ĵ����ĳ�ʼ״̬��������ΪC1~Cn��c0Ĭ��Ϊ1) 
N = length(m);
Tc = Tb/N;
Tm_sample = T_sample/N;  %m���еĲ�������Ϊ10
m1 =  SAM_0(m,Tm_sample,length(m),Tc); %������ʱ��
m2 = 1-2*m1;  %��˫����
NP = 20000;
mm = repmat(m,1,NP); %������
mm1 = SAM_0(mm,Tm_sample,length(mm),Tc);%�����ڣ��г���ʱ��
mm2 = 1-2*mm1;       %��˫����
L = length(mm);      %NP�����ڵ�m���ܳ�
dt = Tc/Tm_sample;
mseq = mm2(1:length(s_initial1));
%-----------------------3����Ƶ--------------------------------------------
spreadData = mseq.*s_initial1;
%-----------------------4������--------------------------------------------
%������Ϊ�޳���ʱ���
downData= zeros(1,length(spreadData)/Tm_sample);
j = 1;
for i=1:Tm_sample:length(spreadData)-9
    downData(j) = spreadData(i);
    j = j+1;
end
spread = downData;
%����ת��
fc=300e6; %�ز�
len_spread = length(spread);
% %---
% if mod(len_spread,2)==1
%     len_spread = len_spread-1;
% end
% %----
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
P_OR = cha4_mrake(sample,length(txx)); %������Ƶ�ʣ����г��ȣ����У�������
pp1 = P_OR(1,:); %ÿһ�ױ���   length(rx)����
pp2 = P_OR(2,:);
pp3 = P_OR(3,:); 
pp4 = P_OR(4,:);
EbN0db  = 0:1:30;
for snr = 1:length( EbN0db ) 
       rx = awgn(txx,snr,'measured'); 
       k11 = 1;
        tap1 = pp1.* rx;
        tap2 = pp2.* rx;
        tap3 = pp3.* rx;   
        tap4 = pp4.* rx;
    %-----------------------5�����--------------------------------------------
    len_spread = length(spread);
    receive1 = demod_nojudge(len_spread,data,t,tap1);
    receive2 = demod_nojudge(len_spread,data,t,tap2);
    receive3 = demod_nojudge(len_spread,data,t,tap3);
    receive4 = demod_nojudge(len_spread,data,t,tap4);
    %��ԭ���г���ʱ�����Ƭ,ʹ�ظ�
    r1 = REPEAT(receive1,Tm_sample); %10
    r2 = REPEAT(receive2,Tm_sample);
    r3 = REPEAT(receive3,Tm_sample);
    r4 = REPEAT(receive4,Tm_sample);
    %-----------------------5-6��ģ���ŵ�����ʱ---------------------------------
    % rr3 = r3;
    % rr2 = r2;
    rr1 = r1;
    % [a1,f] = xcorr(rr1,spreadData,'coeff');
    % subplot(311);plot(f,(a1));title('��һ��������غ���');axis([-500,500,-0.2,1.2]);
    rr2 = r2;   
    rr2(1,T_sample/5+1:end) = rr2(1,1:end-T_sample/5);    %rr2(1,1:T_sample/5)=rr2(1,end-T_sample/5+1:end);
    % [a1,f] = xcorr(rr2,spreadData,'coeff');
    % subplot(312);plot(f,(a1));title('����1/4��Ԫʱ����غ���');axis([-500,500,-0.2,1.2]);
    rr3 = r3;    
    rr3(1,T_sample/2+1:end) = rr3(1,1:end-T_sample/2);    %rr3(1,1:T_sample/2)=rr3(1,end-T_sample/2+1:end);
    % [a1,f] = xcorr(rr3,spreadData,'coeff');
    % subplot(313);plot(f,(a1));title('����1/2��Ԫʱ����غ���');axis([-500,500,-0.2,1.2]);
    rr4 = r4;    %0.002us��1/2��Ԫ
    rr4(1,T_sample/2+1:end) = rr4(1,1:end-T_sample/2);
    %-----------------------6������--------------------------------------------
    %%%%%%%%����
    succ2 = 0;succ3 = 0;succ4 = 0;
    for temp = 1:5:T_sample
        comp = spreadData;
        comp(1,temp+1:end) = comp(1,1:end-temp);
        comp(1,1:temp)= comp(1,end-temp+1:end);
        [a2,f2] = xcorr(rr2,comp,'coeff');
        [a3,f3] = xcorr(rr3,comp,'coeff');
        [a4,f4] = xcorr(rr4,comp,'coeff');
        %figure;plot(f,(y));%axis([-500,500,-0.2,1.2]);
        y2 = max(a2);
        y3 = max(a3);
        y4 = max(a4);
        id2 = find(y2==a2);
        id3 = find(y3==a3);
        id4 = find(y4==a4);
        x2 = f2(id2); %�ҵ���ͼ�ж�Ӧ�����ֵ������ֵ
        x3 = f3(id3);
        x4 = f4(id4);
        if (x2<3 && x2>-3 && y2>0.8)    %��������ֵ
            succ2 = temp;
        end
        if (x3<5 && x3>-5 && y3>0.8)    %��������ֵ
            succ3 = temp;
        end
        if (x4<5 && x4>-5 && y4>0.8)    %��������ֵ
            succ4 = temp;
        end
    end
    %%%%%%%%����
    %----------------------6.2 m������Ӧλ��------------------------------------
    demseq1 = mseq;
    demseq2 = mseq;
    demseq3 = mseq;
    demseq4 = mseq;
    demseq2(1,succ2+1:end) = demseq2(1,1:end-succ2);%demseq2(1,1:succ2)=demseq2(1,end-succ2+1:end);
    demseq3(1,succ3+1:end) = demseq3(1,1:end-succ3);%demseq3(1,1:succ3)=demseq3(1,end-succ3+1:end);
    demseq4(1,succ4+1:end) = demseq4(1,1:end-succ4);
    demp1 = demseq1.*rr1;
    demp2 = demseq2.*rr2;
    demp3 = demseq3.*rr3;
    demp4 = demseq4.*rr4;

    %-----------------------6.3 ������λ����(��λ���ֻ�ԭ)������ԭ��-----------------------------------
    align1 = demp1;
    align2 = demp2;
    align3 = demp3;
    align4 = demp4;
    align2(1,1:end-succ2) = align2(1,succ2+1:end);
    align3(1,1:end-succ3) = align3(1,succ3+1:end);
    align4(1,1:end-succ4) = align4(1,succ4+1:end);

    %����
    dempsam1 = zeros(1,length(align1)/T_sample);%demp����ԭʼ�Ķ�Ӧ����310����
    dempsam2 = dempsam1;%��ʼ��
    dempsam3 = dempsam1;
    dempsam4 = dempsam1;
    j = 1;
    for i=1:T_sample:length(align1)-T_sample+1
        dempsam1(j) = align1(i);
        dempsam2(j) = align2(i);
        dempsam3(j) = align3(i);
        dempsam4(j) = align4(i);
        j = j+1;
    end

    %-----------------------6~7�����ֺϲ�--------------------------------------
        p1 = dempsam1*dempsam1';   %����ÿһ���Ĺ���
        p2 = dempsam2*dempsam2';
        p3 = dempsam3*dempsam3';
        p4 = dempsam4*dempsam4';
        p = p1 + p2 + p3 + p4;
        u1 = p1/p;             %��ÿһ���ļ�Ȩϵ��
        u2 = p2/p;
        u3 = p3/p; 
        u4 = p4/p; 
        %���Ⱥϲ�
        rd_m1 = real( dempsam1*u1+dempsam2*u2+dempsam3*u3+dempsam4*u4);
        %������ϲ�
        rd_m2 = (real(dempsam1+dempsam2+dempsam3+dempsam4))/3;
        %ѡ��ʽ�ϲ�
        u = [u1,u2,u3,u4];   %��ѭ����⣩
        maxu = max(u);
        if(maxu==u1)
            rd_m3 = real(dempsam1);end
        if(maxu==u2)
            rd_m3 = real(dempsam2);end
        if(maxu==u3)
            rd_m3 = real(dempsam3);end
        if(maxu==u4)
            rd_m3 = real(dempsam4);end
        %���ַ����о����
        r_Data1 = sign(rd_m1); %����0Ϊ1��С��0Ϊ-1  ��Ҫת���ˣ�1*������
        r_Data2 = sign(rd_m2);
        r_Data3 = sign(rd_m3);
        %�����������
        Bit_Error_Number1 = length(find(r_Data1(1:Tlen) ~= s_initial(1:Tlen)));     %��ԭʼ�źűȽ�
        Bit_Error_Rate1(snr) = Bit_Error_Number1/Tlen;
        Bit_Error_Number2 = length(find(r_Data2(1:Tlen) ~= s_initial(1:Tlen)));
        Bit_Error_Rate2(snr) = Bit_Error_Number2/Tlen;
        Bit_Error_Number3 = length(find(r_Data3(1:Tlen) ~= s_initial(1:Tlen)));
        Bit_Error_Rate3(snr) = Bit_Error_Number3/Tlen;
    %     %---------------------------����---------------------------------------
    %     [final1,node] = decode_1_3(dempsam1,length(g),memory_els,Tlen,0);%���������˫���Եģ�֮ǰ�õ�r_Data1�������ĳ��о�ǰ��rd_m1
    % %     [final2,node] = decode_1_3(r_Data2,length(g),memory_els,Tlen,0);
    % %     [final3,node] = decode_1_3(r_Data3,length(g),memory_els,Tlen,0);
    %     Bit_Error_Number11 = length(find(final1(1:Tlen) ~= initial(1:Tlen)));     %��ԭʼ�źűȽ�
    %     Bit_Error_Rate11(snr) = Bit_Error_Number11/Tlen;
    % %     Bit_Error_Number22 = length(find2(final2(1:Tlen) ~= initial(1:Tlen)));     %��ԭʼ�źűȽ�
    % %     Bit_Error_Rate22(snr) = Bit_Error_Number22/Tlen;
    % %     Bit_Error_Number33 = length(find3(fina3(1:Tlen) ~= initial(1:Tlen)));     %��ԭʼ�źűȽ�
    % %     Bit_Error_Rate33(snr) = Bit_Error_Number33/Tlen;
end

%-----------------------------7����������-----------------------------------
figure
semilogx(EbN0db,Bit_Error_Rate1,'ro-');hold on; 
semilogx(EbN0db,Bit_Error_Rate2,'bo-');hold on;
semilogx(EbN0db,Bit_Error_Rate3,'go-');hold on;
% semilogy(EbN0db,Bit_Error_Rate11,'r*-');hold on;
% semilogy(EbN0db,Bit_Error_Rate22,'b*-');hold on;
% semilogy(EbN0db,Bit_Error_Rate32,'g*-');hold on;
legend('���Ⱥϲ�','������ϲ�','ѡ��ʽ�ϲ�','*����������ֺϲ�');
xlabel('�����');ylabel('������');
title('�ľ�����Ҫ�ּ��ϲ���ʽ���ܱȽ�');
grid  on;



