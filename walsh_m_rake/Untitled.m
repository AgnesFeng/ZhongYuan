a = [1,2,3];
b = [1,2,3,4,5,6,7];
xcor = dsp.Crosscorrelator;
yy = step(xcor,a',b')
    despreadData = [1,1,1,1,2,2,2,2,3,3,3,3,];
   T_sample = 4;
    receive = ones(1,3);
    k=1;
    for i = 1:T_sample:length(despreadData)-T_sample+1
        receive(k) = 0;
        for j = 0:T_sample-1
          receive(k) = receive(k) + despreadData(i+j);
        end
        k = k+1;
    end
    receive
    
a = [1,2,3,4,5,6,7];
B = circshift(a,[0,-2])%[]中的第一个数为0，表示行，第二个数为正是右移位数，为负是左移位数