function x = get_parabolic_max_from_3_pts(x1,x2,x3,y1,y2,y3)
%note: works also for a minimum

dx1 = x2-x1;
dx2 = x2 - x3;
dy1 = y2-y1;
dy2 = y2-y3;

x = x2 - 0.5 * ( dx1^2* dy2 - dx2^2*dy1 )/(dx1* dy2 - dx2*dy1);

end