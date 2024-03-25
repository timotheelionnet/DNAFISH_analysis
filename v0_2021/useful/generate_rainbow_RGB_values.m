function [ r,g,b ] = generate_rainbow_RGB_values( n )
%generates vectors r, g, b, each with n values
%each triplet r(i),g(i),b(i) encodes a separate color; the range of colors
%spans a rainbow. 

%useful to plot many curves in one plot
    
n = n-1;    
r(n,1) = 0;
g(n,1) = 0;
b(n,1) = 0;

dx = 0.8;

for f=0:(1/n):1
    gg = (6-2*dx)*f+dx;
    index = int16(f*n + 1);
    r(index,1) = max(0,(3-abs(gg-4)-abs(gg-5))/2);
    g(index,1) = max(0,(4-abs(gg-2)-abs(gg-4))/2); 
    b(index,1) = max(0,(3-abs(gg-1)-abs(gg-2))/2);

end

    

end

