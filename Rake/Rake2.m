clear all;
Numusers = 1;
Nc = 16; %扩频因子
ISI_Length = 1;%每径延时为ISI_Length/2
EbN0db = [0:1:30];%信噪比，单位db
Tlen = 20000; %数据长度 （取越多越好）
Bit_Error_Number1 = 0;%误比特率的初始值
Bit_Error_Number2 = 0;
Bit_Error_Number3 = 0;
power_unitary_factor1 = sqrt( 5/9 ); %每径功率因子
power_unitary_factor2 = sqrt( 3/9 );
power_unitary_factor3 = sqrt( 1/9 );
%s_initial = randsrc( 1, Tlen ); %数据源 (1*8000的矩阵，-1,1各占一半)
s_initial = randsrc( 1, Tlen );
%%%%%%%%%%%%%%%%%%%%%%%%产生walsh矩阵%%%%%%%%%%%%%%%%%%%%%%
wal2 = [ 1 1; 1 -1 ];
wal4 = [wal2 wal2; wal2 wal2*(-1)];  
wal8 = [wal4 wal4; wal4 wal4*(-1)];  %8*8
wal16 = [wal8 wal8; wal8 wal8*(-1)]; %16*16
%%%%%%%%%%%%%%%%%%%%%%%%%%%扩频%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
s_spread = zeros( Numusers, Tlen*Nc );  %[1,8000*16]扩频后的码片 初始化
ray1 = zeros( Numusers, 2*Tlen*Nc ); 
ray2 = zeros( Numusers, 2*Tlen*Nc );
ray3 = zeros( Numusers, 2*Tlen*Nc );

for i = 1:Numusers
    x0 = s_initial( i,: ).'*wal16( 8,: );  %经过转置，行变列，乘上8*8矩阵的第八行
    x1 = x0.';%变成16*8000的矩阵
    s_spread( i,: ) = ( x1(:) ).';%得到扩频后的序列，【x1(:)从头到尾排序，即16*8000个】
end
%将每个扩频后的输出重复两次，有利于下面的延迟（延迟半个码元）
ray1( 1:2:2*Tlen*Nc - 1 ) = s_spread( 1:Tlen*Nc );%相邻的即是重复的
ray1( 2:2:2*Tlen*Nc ) = ray1( 1:2:2*Tlen*Nc - 1 ); 
%…对第一径信号进行两次延迟，分别得到第二径和第三经…%
ray2( ISI_Length + 1:2*Tlen*Nc ) = ray1( 1:2*Tlen*Nc - ISI_Length );
ray3( 2*ISI_Length + 1:2*Tlen*Nc ) = ray1( 1:2*Tlen*Nc - 2*ISI_Length ); 
%进行解扩
for nEN = 1:length( EbN0db )    
    en = 10^( EbN0db(nEN)/10 ); %将Eb/N0的db值转化成十进制数值
    sigma = sqrt( 32/(2*en) );%根据初始状态设定的信噪比，计算出AWGN方差，32为比特能量，与扩频调制方式有关，扩频后每径32bit对应一个原始bit，比特能能量是原来的32倍
    %接收到的信号demp
    demp = power_unitary_factor1*ray1+...           %每径乘功率因子再叠加
    power_unitary_factor2*ray2+...                  %注意下一行rand和randn的区别
    power_unitary_factor3*ray3 + ( rand( 1,2*Tlen*Nc )+randn( 1,2*Tlen*Nc )*i )*sigma;  %最后一项为噪声
                                   %rand:0~1之间随机数             randn:正态分布的随机数
    dt = reshape( demp,32,Tlen )';   %dt为8000*32的demp值的矩阵
    %将walsh码重复两次
    wal16_d(1:2:31) = wal16(8,1:16); %用刚才扩频的序列填充
    wal16_d(2:2:32) = wal16(8,1:16);
    
    %解扩后rdata1为第一径输出
    rdata1 = dt*wal16_d(1,:).';
    %将walsh码延迟半个码片
    wal16_delay1(1,2:32) = wal16_d(1,1:31);
    %wal16_delay1(1) =wal16_d(32);%自己加的，避免有0，尔出现三但结果偶条线交叉，其实前面是不是0无所谓，因为取得后面的数相乘，后来发现不是去后面的数，而是算个总功率，
    %解扩后rdata2为第二径输出
    rdata2 = dt*wal16_delay1(1,:).';
    %将walsh码延迟一个码片
    wal16_delay2(1,3:32) = wal16_d(1,1:30);%（矩阵的第一个下标为1不是0）
    wal16_delay2(1,1:2) = wal16_d(1,31:32);%（最后两个值放到最前，避免有0）
    %解扩后rdata3为第三径输出
    rdata3 = dt*wal16_delay2(1,:).';
    
    p1 = rdata1'*rdata1;   %计算每一径的功率（有问题，为什么比特的平方就是幅度？）
    p2 = rdata2'*rdata2;
    p3 = rdata3'*rdata3;
    p = p1 + p2 + p3;
    u1 = p1/p;             %算每一径的加权系数
    u2 = p2/p;
    u3 = p3/p; 
    %最大比合并
    rd_m1 = real( rdata1*u1+rdata2*u2+rdata3*u3);
    %等增益合并
    rd_m2 = (real(rdata1+rdata2+rdata3))/3;
    %选择式合并
    u = [u1,u2,u3];   %（循环检测）
    maxu = max(u);
    if(maxu==u1)
        rd_m3 = real(rdata1);
      else if(maxu==u2)
           rd_m3 = real(rdata2);
      else
           rd_m3 = real(rdata3);
      end
    end %三种方法判决输出
    r_Data1 = sign(rd_m1)'; %大于0为1，小于0为-1
    r_Data2 = sign(rd_m2)';
    r_Data3 = sign(rd_m3)';
    %计算误比特率
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
legend('最大比合并','等增益合并','选择式合并');
xlabel('信噪比');
ylabel('误比特率');
title('三种主要分集合并方式性能比较');