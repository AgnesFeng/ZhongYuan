% 
initial = importdata('initial.mat');
m_seq = importdata('m_seq.mat');
xcor = dsp.Crosscorrelator;
feng1 = step(xcor,initial,m_seq);
feng1 = abs(feng1);
figure,plot(feng1);
% 
s1 = initial'; 
shift = 100;
s1(1,shift+1:end) = s1(1,1:end-shift);
feng2 = step(xcor,s1',m_seq);
feng2 = abs(feng2);
figure,plot(feng2);

s2 = s1;
s2(1,1:end-shift) = s2(1,shift+1:end);
feng3 = step(xcor,s2',m_seq);
feng3 = abs(feng3);
figure,plot(feng3);

tempxx = [ 1028        1032        1084        1089        1080 ];
s = ones(1,5);
for i = 1:5
    s(i) = tempxx(i)-1020;
end

path_1 = initial';
path_2 = initial';
path_3 = initial';
path_4 = initial';
path_5 = initial';
path_1(1,1:end-s(1)) = path_1(1,s(1)+1:end);
path_2(1,1:end-s(2)) = path_2(1,s(2)+1:end);
path_3(1,1:end-s(3)) = path_3(1,s(3)+1:end);
path_4(1,1:end-s(4)) = path_4(1,s(4)+1:end);
path_5(1,1:end-s(5)) = path_5(1,s(5)+1:end);
p = y1(tempxx(1)) + y1(tempxx(2)) + y1(tempxx(3))+y1(tempxx(4));   %计算每一径的加权系数
u1 = y1(tempxx(1))/p;             
u2 = y1(tempxx(2))/p;
u3 = y1(tempxx(3))/p; 
u4 = y1(tempxx(4))/p;
u5 = y1(tempxx(5))/p;
merge = path_1 * u1 + path_2 * u2 + path_3 * u3 + path_4 * u4;

y = step(xcor,merge',m_seq); %computes cross-correlation of x1 and x2
y2 = abs(y);
figure(22), plot(y2); title('4径');

merge = path_1 * u1 + path_2 * u2 + path_3 * u3 + path_4 * u4 + path_5 * u5;
y = step(xcor,merge',m_seq); %computes cross-correlation of x1 and x2
y2 = abs(y);
figure(33), plot(y2); title('5径');
% 
% 
a = [1+1i*23 2+1i*12 3+1i*4 4+1i*5];
b = [1,2,3,4];
xcor = dsp.Crosscorrelator;
feng1 = step(xcor,a ,b)

