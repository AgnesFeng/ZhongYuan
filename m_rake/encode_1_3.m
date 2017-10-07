function c = encode_1_3(m,g,n)
%Function:  Encodes a Rate 1/3 Convolution Code
%m = message to encode
%g = n generators corresponding to the n outputs
%n = 1/n Convolution Encoder -- # of generators

%First Perform Convolution of Input Message for Each Generator. This
%produces n outputs... [y{1},y{2},...,y{n}]
for i = 1:n
    y{i} = mod(conv(m,g{i}),2);%卷积
end

%Initialize code word to all zeros
c = zeros(1,n*length(y{1}));

%Assemble code word from n outputs
for i =1:n
    c(i:n:end) = y{i};%并串转换后输出
end


