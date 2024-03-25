function err = compute_error_on_ratio(a,b,ea,eb)

err = a./b .* sqrt( ea.^2./a.^2 + eb.^2./b.^2 );

end