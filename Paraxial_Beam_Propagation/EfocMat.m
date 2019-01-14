% calcs focal matrix for electric field of thin shell
% E0, x0 and K in "osiris units"
% r is inner radius of thin shell
% R is outer radius
function M = EfocMat(E0, x0, K, r, R)

intEdl = 2*x0*E0.*(asinh(sqrt(R.^2./x0.^2-1.0).*(x0.^2<R.^2)) ...
- asinh(sqrt(r.^2./x0.^2-1.0)).*(x0.^2<r.^2));

invf = -intEdl./(2.0.*K.*x0);
if (~isfinite(invf))
    invf = 0.0;
end

M = [1.0 0.0
     invf 1.0];

end
