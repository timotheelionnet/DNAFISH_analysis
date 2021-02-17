function l = compute_fit_likelihood(x,y,s,fun)

    if numel(y) ~= numel(x) || numel(s) ~= numel(x) || ~isa(fun,'function_handle')
        l = NaN;
        return 
    end
    
    l = prod( 1./sqrt(2*pi*s.^2) .* exp( -( y - fun(x)).^2./(2*s.^2) ) ) * prod( sqrt(2*pi*s.^2));
end





