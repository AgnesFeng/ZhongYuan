clear all;
Numusers = 1;
Nc = 16; %��Ƶ����
ISI_Length = 1;%ÿ����ʱΪISI_Length/2
EbN0db = [0:1:30];%����ȣ���λdb
Tlen = 20000; %���ݳ��� ��ȡԽ��Խ�ã�
Bit_Error_Number1 = 0;%������ʵĳ�ʼֵ
Bit_Error_Number2 = 0;
Bit_Error_Number3 = 0;
power_unitary_factor1 = sqrt( 5/9 ); %ÿ����������
power_unitary_factor2 = sqrt( 3/9 );
power_unitary_factor3 = sqrt( 1/9 );
%s_initial = randsrc( 1, Tlen ); %����Դ (1*8000�ľ���-1,1��ռһ��)
s_initial = randsrc( 1, Tlen );
%%%%%%%%%%%%%%%%%%%%%%%%����walsh����%%%%%%%%%%%%%%%%%%%%%%
wal2 = [ 1 1; 1 -1 ];
wal4 = [wal2 wal2; wal2 wal2*(-1)];  
wal8 = [wal4 wal4; wal4 wal4*(-1)];  %8*8
wal16 = [wal8 wal8; wal8 wal8*(-1)]; %16*16
%%%%%%%%%%%%%%%%%%%%%%%%%%%��Ƶ%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
s_spread = zeros( Numusers, Tlen*Nc );  %[1,8000*16]��Ƶ�����Ƭ ��ʼ��
ray1 = zeros( Numusers, 2*Tlen*Nc ); 
ray2 = zeros( Numusers, 2*Tlen*Nc );
ray3 = zeros( Numusers, 2*Tlen*Nc );

for i = 1:Numusers
    x0 = s_initial( i,: ).'*wal16( 8,: );  %����ת�ã��б��У�����8*8����ĵڰ���
    x1 = x0.';%���16*8000�ľ���
    s_spread( i,: ) = ( x1(:) ).';%�õ���Ƶ������У���x1(:)��ͷ��β���򣬼�16*8000����
end
%��ÿ����Ƶ�������ظ����Σ�������������ӳ٣��ӳٰ����Ԫ��
ray1( 1:2:2*Tlen*Nc - 1 ) = s_spread( 1:Tlen*Nc );%���ڵļ����ظ���
ray1( 2:2:2*Tlen*Nc ) = ray1( 1:2:2*Tlen*Nc - 1 ); 
%���Ե�һ���źŽ��������ӳ٣��ֱ�õ��ڶ����͵�������%
ray2( ISI_Length + 1:2*Tlen*Nc ) = ray1( 1:2*Tlen*Nc - ISI_Length );
ray3( 2*ISI_Length + 1:2*Tlen*Nc ) = ray1( 1:2*Tlen*Nc - 2*ISI_Length ); 
%���н���
for nEN = 1:length( EbN0db )    
    en = 10^( EbN0db(nEN)/10 ); %��Eb/N0��dbֵת����ʮ������ֵ
    sigma = sqrt( 32/(2*en) );%���ݳ�ʼ״̬�趨������ȣ������AWGN���32Ϊ��������������Ƶ���Ʒ�ʽ�йأ���Ƶ��ÿ��32bit��Ӧһ��ԭʼbit��������������ԭ����32��
    %���յ����ź�demp
    demp = power_unitary_factor1*ray1+...           %ÿ���˹��������ٵ���
    power_unitary_factor2*ray2+...                  %ע����һ��rand��randn������
    power_unitary_factor3*ray3 + ( rand( 1,2*Tlen*Nc )+randn( 1,2*Tlen*Nc )*i )*sigma;  %���һ��Ϊ����
                                   %rand:0~1֮�������             randn:��̬�ֲ��������
    dt = reshape( demp,32,Tlen )';   %dtΪ8000*32��dempֵ�ľ���
    %��walsh���ظ�����
    wal16_d(1:2:31) = wal16(8,1:16); %�øղ���Ƶ���������
    wal16_d(2:2:32) = wal16(8,1:16);
    
    %������rdata1Ϊ��һ�����
    rdata1 = dt*wal16_d(1,:).';
    %��walsh���ӳٰ����Ƭ
    wal16_delay1(1,2:32) = wal16_d(1,1:31);
    %wal16_delay1(1) =wal16_d(32);%�Լ��ӵģ�������0���������������ż���߽��棬��ʵǰ���ǲ���0����ν����Ϊȡ�ú��������ˣ��������ֲ���ȥ�����������������ܹ��ʣ�
    %������rdata2Ϊ�ڶ������
    rdata2 = dt*wal16_delay1(1,:).';
    %��walsh���ӳ�һ����Ƭ
    wal16_delay2(1,3:32) = wal16_d(1,1:30);%������ĵ�һ���±�Ϊ1����0��
    wal16_delay2(1,1:2) = wal16_d(1,31:32);%���������ֵ�ŵ���ǰ��������0��
    %������rdata3Ϊ���������
    rdata3 = dt*wal16_delay2(1,:).';
    
    p1 = rdata1'*rdata1;   %����ÿһ���Ĺ��ʣ������⣬Ϊʲô���ص�ƽ�����Ƿ��ȣ���
    p2 = rdata2'*rdata2;
    p3 = rdata3'*rdata3;
    p = p1 + p2 + p3;
    u1 = p1/p;             %��ÿһ���ļ�Ȩϵ��
    u2 = p2/p;
    u3 = p3/p; 
    %���Ⱥϲ�
    rd_m1 = real( rdata1*u1+rdata2*u2+rdata3*u3);
    %������ϲ�
    rd_m2 = (real(rdata1+rdata2+rdata3))/3;
    %ѡ��ʽ�ϲ�
    u = [u1,u2,u3];   %��ѭ����⣩
    maxu = max(u);
    if(maxu==u1)
        rd_m3 = real(rdata1);
      else if(maxu==u2)
           rd_m3 = real(rdata2);
      else
           rd_m3 = real(rdata3);
      end
    end %���ַ����о����
    r_Data1 = sign(rd_m1)'; %����0Ϊ1��С��0Ϊ-1
    r_Data2 = sign(rd_m2)';
    r_Data3 = sign(rd_m3)';
    %�����������
    Bit_Error_Number1 = length(find(r_Data1(1:Tlen) ~= s_initial(1:Tlen)));
    Bit_Error_Rate1(nEN) = Bit_Error_Number1/Tlen;
    Bit_Error_Number2 = length(find(r_Data2(1:Tlen) ~= s_initial(1:Tlen)));
    Bit_Error_Rate2(nEN) = Bit_Error_Number2/Tlen;
    Bit_Error_Number3 = length(find(r_Data3(1:Tlen) ~= s_initial(1:Tlen)));
    Bit_Error_Rate3(nEN) = Bit_Error_Number3/Tlen;
end
semilogy(EbN0db,Bit_Error_Rate1,'r*-');hold on;
semilogy(EbN0db,Bit_Error_Rate2,'bo-');hold on;
semilogy(EbN0db,Bit_Error_Rate3,'g.-');
legend('���Ⱥϲ�','������ϲ�','ѡ��ʽ�ϲ�');
xlabel('�����');
ylabel('�������');
title('������Ҫ�ּ��ϲ���ʽ���ܱȽ�');