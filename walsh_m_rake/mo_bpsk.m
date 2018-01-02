function b=mo_bpsk(code,Rb,Tm_sample)
%ÂëÔªËÙÂÊ¡¢
% code1=randsrc(1,10);
% Rb=225e3;
% Tm_sample = 4;
% code =  SAM_d(code1,Tm_sample,length(code1));

fs=Tm_sample * Rb;
dt = 1/fs;
N=length(code);
k1=0;
b=zeros(1,N);
% for i=1:k
%     for j=k1:k1+Tm_sample-1
%         if code(1,i)==1
%             b(1,j+1)=cos(2*pi*Rb*j/fs);
%         else 
%             b(1,j+1)=cos(2*pi*Rb*j/fs+pi);
%         end
%     end
%     k1=k1+Tm_sample;
% end

for i=1:Tm_sample:N-Tm_sample+1
    for j = 1:Tm_sample
        if code(i)==1
        b(k1+j)=cos(2*pi*Rb*j*dt);
        else 
        b(k1+j)=cos(2*pi*Rb*j*dt+pi);
        end
    end
    k1=k1+Tm_sample;
end
%figure,plot(b);
end