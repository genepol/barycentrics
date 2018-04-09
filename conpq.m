%conpq does a least square fit of a conic through several points
xyvals = input('[x(1) x(2 ...x(t);y(1) y(2)...y(t)]');
t = length(xyvals);
M = zeros(2*t,6);
x = xyvals(1,1:t); y = xyvals(2,1:t);
for s = 1:t-1
    u(s) = (x(s)+x(s+1))/2; v(s) = (y(s)+y(s+1))/2;
    M(s,:) = 2*[1 x(s) y(s) x(s)^2 x(s)*y(s) y(s)^2];
    M(t+s,:) = [1 u(s) v(s) u(s)^2 u(s)*v(s) v(s)^2];
end
M(t,:) = 10*[1 x(t) y(t) x(t)^2 x(t)*y(t) y(t)^2];
x0= (x(1)+x(t))/2; y0 = (y(1)+y(t))/2;
M(2*t,:) = [1 x0 y0 x0^2 x0*y0 y0^2];
M(1,:) = 10*M(1,:);
rhs = [zeros(2*t-1,1);-1.0];
rr = M'*rhs;
A = M'*M;
cofs =A\rr;
syms p q
p6q = [1 p q p^2 p*q q^2];
confit = p6q*cofs;
conv = vpa(confit,6);
pretty(conv)
%When conv is a hyperbola we replace it with a parabola.
%This avoids an ill-set element due to the other branch.
disc = cofs(5)^2 - 4*cofs(4)*cofs(6);
    if  disc > 0
        tanc = (x(t)-x(1))*diff(conv,p) + (y(t)-y(1))*diff(conv,q);
        twopts = solve(tanc,conv);
        pp = double(twopts.p); qq = double(twopts.q);
            if length(pp) == 1
                x2 = pp; y2 = qq;
            else
                x1t = (x(1) + x(t))/2; y1t = (y(1)+y(t))/2;
                d1 = (pp(1)-x1t)^2 + (qq(1) - y1t)^2;
                d2 = (pp(2)-x1t)^2 + (qq(2) - y1t)^2;
                if d1 < d2
                    x2 = pp(1); y2 = qq(1);
                else
                    x2 = pp(2); y2 = qq(2);
                end               
            end
        [real(x2);real(y2)]
    end
 return