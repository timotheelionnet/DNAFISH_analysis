function [x,y,e] = linear_interpolation_with_error(xi,yi,ei,x)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here


%INPOUT
%xi: vector array of x axis values
%yi: corresponding y axis values 
%ei: corresponding uncertainty values on yi
%desired x axis coordinates the function will be computed on 
%(values outside the initial range will be interpolated as NaN)

%OUTPUT
%x: array of x axis values sorted
%y: correspondinf y values interpolated linearly: 
    %y = y1 + (y2-y1)/(x2-x1)*(x-x1)
    %where x1 and x2 are the two points from the initial dataset flanking
    %each x value in the interpolated dataset
%e: corresponding error computed using linear propagation, assuming
%independent variables: e = sqrt( (dy/dy1)^2*e1^2 + (dy/dy2)^2*e2^2  )

%reformat as vector columns if needed
xi = reshape(xi,numel(xi),1);
yi = reshape(yi,numel(yi),1);
ei = reshape(ei,numel(ei),1);
x = reshape(x,numel(x),1);

%make sure all inputs are consistent in size
if ~isequal(size(xi),size(yi))
    disp('xi and yi arrays should have the same size');
    y = NaN*zeros(size(x));
    e = NaN*zeros(size(x));
    return
end

if ~isequal(size(xi),size(ei))
    disp('xi and ei arrays should have the same size');
    y = NaN*zeros(size(x));
    e = NaN*zeros(size(x));
    return
end

%sort input data along x
[xi,Idx] = sort(xi,'ascend');
yi = yi(Idx);
ei = ei(Idx);

%sort output coordinates along x
x = sort(x,'ascend');

%fill in values
y = zeros(size(x));
e = zeros(size(x));
for i=1:numel(x)
   xtmp = xi(xi<=x(i));
   
   if isempty(xtmp)
       y(i) = NaN;
       e(i) = NaN;
   else
       [x1,i1] = max(xtmp(:));
       y1 = yi(i1);
       e1 = ei(i1);

       if x(i) == x1
           y(i) = y1;
           e(i) = e1;
       else
           xtmp = xi;
           xtmp(xi<=x(i)) = Inf;
           [x2,i2] = min(xtmp(:));
           if isinf(xtmp)
                y(i) = NaN;
                e(i) = NaN;
           else
               y2 = yi(i2);
               e2 = ei(i2);

               y(i) = y1 + (y2-y1)/(x2-x1)*(x(i) - x1);
               e(i) = sqrt( (1 + (x1-x(i))/(x2-x1))^2*e1^2 +  ((x(i)-x1)/(x2-x1))^2*e2^2  );
           end
       end  
   end   
end

end

