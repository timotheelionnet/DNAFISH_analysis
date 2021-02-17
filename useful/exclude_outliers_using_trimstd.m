function [dataclean,idx] = exclude_outliers_using_trimstd(data)

%generate the std of the x% centermost data as a function of x = 1:99
[m,x] = trimstd_fct(data,1:99);

%fit that function to 2 connected lines to find the best separation point
%bewteen two regimes
chi2 = zeros(96,1);
x0= 2.5 : 1 : 97.5;
for i=1:96
    
    %find the best two lines fitting the data (separation at x0)
    %y = ax1+b   (x1<x0)
    %y = cx2+(a-c)x0 +b    (x1>=x0)    
    x1 = 1:floor(x0(i));
    x2 = ceil(x0(i)):99;
    y1 = m(x1);
    y2 = m(x2);
    
    Sx1y1 = sum(x1.*y1);
    Sx2y2 = sum(x2.*y2);
    Sx1x1 = sum(x1.^2);
    Sx2x2 = sum(x2.^2);
    Sx1 = sum(x1);
    Sx2 = sum(x2);
    Sy1 = sum(y1);
    n1 = numel(x1);
    n2 = numel(x2);
    
    p = [Sx1x1,Sx1,0; Sx1,n1,0; x0(i)*Sx2, Sx2, Sx2x2 - x0(i)*Sx2]^(-1) * [Sx1y1;Sy1;Sx2y2];
    
    %compute the error
    chi2(i) = sum( (y1 - p(1)*x1-p(2) ).^2 )/n1  +  sum( (y2 - p(3)*x2-(p(1)-p(3))*x0(i) -p(2) ).^2 )/n2; 
    
    %test iteratively for minimum
    if i==1 
        chi2min = chi2(1);
        a = p(1);
        b = p(2);
        c = p(3);
        m0 = floor(x0(1));
    else
        if chi2(i) < chi2min
            chi2min = chi2(i);
            a = p(1);
            b = p(2);
            c = p(3);
            m0 = floor(x0(i));
        end
    end 
end


[dataclean,idx] = sort(data,'ascend');
i1 = round((100 - m0)/200*numel(data));
i2 = round((100 + m0)/200*numel(data));
idx = idx( i1:i2 );
dataclean = dataclean(i1:i2);


end