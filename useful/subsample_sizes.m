function samplesize = subsample_sizes(npts)
%outputs a list of numbers ranging from npts to 2 following a decreasing logarithmic progression:
%npts,npts/2,npts/5,npts/10,npts/20,npts/50 etc

ns = 3*floor(log(npts)/log(10));
samplesize(1:3:ns-2,1) = round( npts*exp(-log(10).*(0:(ns/3-1))) );
samplesize(2:3:ns-1,1) = round( npts/2*exp(-log(10).*(0:(ns/3-1))) );
samplesize(3:3:ns,1) = round( npts/5*exp(-log(10).*(0:(ns/3-1))) );

end