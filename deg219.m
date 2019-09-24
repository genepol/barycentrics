%deg219 follows poly2018 with computation of degree two basis functions.
syms x y adj14 adj12 adj34
si14 = zeros(2,n);
si34 = zeros(2,n);
si12 = zeros(2,n);
%We locate all side nodes.
                for k = 1:n
    kp1 = k + 1;
    km1 = k - 1;
	if k == n
		kp1 = 1;
    elseif k == 1        
		km1 = n;
    end
		x1 = vert(1,k);
		y1 = vert(2,k);
		x2 = vert(1,kp1);
		y2 = vert(2,kp1);
        if degsi(k) == 1
            si12(1,k) = (x1 + x2)/2;
            si12(2,k) = (y1 + y2)/2;
        else
            si12(1,k) = sino(1,k);
            si12(2,k) = sino(2,k);
            x12 = si12(1,k);
            y12 = si12(2,k);
%We now locate two more nodes on side k.
if x1 == x2
    y14 = y1 + .5*(y12 - y1);
    si14(2,k) = y14;
    y34 = y12 + .5*(y2 - y12);
    si34(2,k) = y34;
    siofx = subs(sixy(k),{'y'} ,{y14});
    r = double(solve(siofx,x));
    if length(r) == 1
        x14 = r;
    else
        n1 = (x12 - r(1))^2; n2 = (x12 - r(2))^2;
        if n1 < n2
            x14 = r(1);
        else
            x14 = r(2);
        end
    end
    si14(1,k) = x14;
    siofx = subs(sixy(k),{'y'} ,{y34});
    r = double(solve(siofx,x));
    if length(r) == 1
        x34 = r;
    else
        n1 = (x12 - r(1))^2; n2 = (x12 - r(2))^2;
        if n1 < n2        
            x34 = r(1);
        else
            x34 = r(2);
        end
   end
    si34(1,k) = x34;
elseif y1 == y2
    x14 = x1 + .5*(x12 - x1);
    x34 = x12 + .5*(x2 - x12);
    si14(1,k) = x14;
    si34(1,k) = x34;
    siofy = subs(sixy(k),{'x'} ,{x14});
    r = double(solve(siofy,y));
    if length(r) == 1
        y14 = r;
    else
        n1 = (y12 - r(1))^2; n2 = (y12 - r(2))^2;
        if n1 < n2
            y14 = r(1);
        else
            y14 = r(2);
        end
    end
    si14(2,k) = y14;
    siofy = subs(sixy(k),{'x'} ,{x34});
    r = double(solve(siofy,y));
    if length(r) == 1
        y34 = r;
    else
        n1 = (y12 - r(1))^2; n2 = (y12 - r(2))^2;
        if n1 < n2
           y34 = r(1);
        else
            y34 = r(2);
        end
    end
    si34(2,k) = y34;
else
    x0 = (x1 + x12)/2; y0 = (y1 + y12)/2;
    m0 = (y12 - y1)/(x12 - x1);
    perp = (y - y0) +  m0*(x - x0);
    r = solve(sixy(k),perp);
    r.x = double(r.x);
    r.y = double(r.y);
    if length(r.x) == 1
        si14(1,k) = r.x;
        si14(2,k) = r.y(1);
    else
        n1 = (x12 - r.x(1))^2 + (y12 - r.y(1))^2;
        ty = length(r.y);
        n2 = (x12 - r.x(2))^2 + (y12 - r.y(ty))^2;
        if n1 < n2
            si14(1,k) = r.x(1); 
            si14(2,k) = r.y(1);
        else
            si14(1,k) = r.x(2);
            si14(2,k) = r.y(ty);
        end
    end
    x0 = (x2 + x12)/2; y0 = (y2 +y12)/2;
    m0 = (y12 - y2)/(x12 - x2);
    perp = (y - y0) +  m0*(x - x0);
    r = solve(sixy(k),perp);
    r.x = double(r.x);
    r.y = double(r.y);
    if length(r.x) == 1
        si34(1,k) = r.x;
        si34(2,k) = r.y(1);
    else
        n1 = (x12 - r.x(1))^2 + (y12 - r.y(1))^2;
        ty = length(r.y);
        yt = r.y(ty);
        n2 = (x12 - r.x(2))^2 + (y12 - yt)^2;
        if n1 < n2
            si34(1,k) = r.x(1); 
            si34(2,k) = r.y(1);
        else
            si34(1,k) = r.x(2);
            si34(2,k) = r.y(ty);
        end
    end
end
        end
%All side nodes have been found.
%We now compute the adjacent factor for all side nodes.
    if degsi(k) == 1
        adj12(k) = 1;
        adj14(k) = 1;
        adj34(k) = 1;
    end        
            if degsi(k) == 2
x1 = si14(1,k); y1 = si14(2,k); x2 = si34(1,k); y2 = si34(2,k);
c1 = y1-y2; c2 = x2-x1; c3 = x1*y2 - y1*x2; 
c4 = c1*si12(1,k)+c2*si12(2,k) + c3;
d1 = c1/c4; d2 = c2/c4; d3 = c3/c4;
adj12(k) = d3 + d1*x + d2*y;
adj12(k) = vpa(adj12(k));
	x1 = si12(1,k); y1 = si12(2,k); x2 = si34(1,k); y2 = si34(2,k);
c1 = y1-y2; c2 = x2-x1; c3 = x1*y2 -y1*x2; 
c4 = c1*si14(1,k) + c2*si14(2,k) + c3;
d1 = c1/c4; d2 = c2/c4; d3 = c3/c4;
adj14(k) = d3 + d1*x + d2*y;
adj14(k) = vpa(adj14(k));
	x2 = si12(1,k); y2 = si12(2,k); x1 = si14(1,k); y1 = si14(2,k);
c1 = y1-y2; c2 = x2-x1; c3 = x1*y2 - y1*x2; 
c4 = c1*si34(1,k)+c2*si34(2,k) + c3;
d1 = c1/c4; d2 = c2/c4; d3 = c3/c4;
adj34(k) = d3 + d1*x + d2*y;
adj34(k) = vpa(adj34(k));
            end
                end
%We now compute the adjacent factor for vertex nodes.
                            for k = 1:n
       kp1 = k+1;
	   km1 = k-1;
		if k == 1
		   km1 = n;
		elseif k == n
		   kp1 = 1;
		end
   if degsi(k) == 1 & degsi(km1) == 1
%Both sides are linear so the adjacent factor is linear
		x1 = si12(1,k);
		y1 = si12(2,k);
		x2 = si12(1,km1);
		y2 = si12(2,km1);
c1 = y1-y2; c2 = x2-x1; c3 = x1*y2 -y1*x2; 
c4 = c1*vert(1,k) + c2*vert(2,k) + c3;
d1 = c1/c4; d2 = c2/c4; d3 = c3/c4;
va = d3 + d1*x + d2*y;
vadj2(k) = vpa(va);
   elseif degsi(k) == 1 & degsi(km1) == 2
%side k is linear but side km1 is conic.
%We locate five points on the adjacent conic
nod = [si12(:,k)  si14(:,km1) si12(:,km1) si34(:,km1) eip(2:3,3*(k-1)+1)];
nrm = vert(:,k); % This is the normalization node.
cmat = zeros(6);
rhs = [1;0;0;0;0;0];
cmat(1,:) = [1 nrm(1) nrm(2) nrm(1)^2 nrm(1)*nrm(2) nrm(2)^2];
	for m = 1:4
cmat(m+1,:) = [1 nod(1,m) nod(2,m) nod(1,m)^2 nod(1,m)*nod(2,m) nod(2,m)^2];
	end
m = 5;
		if eip(1,3*(k-1)+1) == 1
cmat(m+1,:) = [1 nod(1,m) nod(2,m) nod(1,m)^2 nod(1,m)*nod(2,m) nod(2,m)^2];
        else
cmat(m+1,:) = [0 0 0 nod(2,m)^2 nod(1,m)*nod(2,m) nod(1,m)^2];
        end
vad = cmat\rhs;
nc = norm(vad,1);
	for km = 1:m+1
if abs(vad(km)) < 1e-5*nc
	vad(km) = 0;
end
	end
va = vad(1) + vad(2)*x + vad(3)*y + vad(4)*x^2 + vad(5)*x*y + vad(6)*y^2;
vadj2(k) = vpa(va);
  elseif degsi(km1) == 1 & degsi(k) == 2
%side km1 is linear but side k is conic
%We locate five points on the adjacent conic
nod = [si12(:,km1) si14(:,k) si12(:,k) si34(:,k) eip(2:3,3*(k-1)+1)];
nrm = vert(:,k); % This is the normalization node.
cmat = zeros(6);
rhs = [1;0;0;0;0;0];
cmat(1,:) = [1 nrm(1) nrm(2) nrm(1)^2 nrm(1)*nrm(2) nrm(2)^2];
	for m = 1:4
cmat(m+1,:) = [1 nod(1,m) nod(2,m) nod(1,m)^2 nod(1,m)*nod(2,m) nod(2,m)^2];
	end
m = 5;
		if eip(1,3*(k-1)+1) == 1
cmat(m+1,:) = [1 nod(1,m) nod(2,m) nod(1,m)^2 nod(1,m)*nod(2,m) nod(2,m)^2];
        else
cmat(m+1,:) = [0 0 0 nod(2,m)^2 nod(1,m)*nod(2,m) nod(1,m)^2];                
        end
vad = cmat\rhs;
nc = norm(vad,1);
	for km = 1:m+1
if abs(vad(km)) < 1e-5*nc
	vad(km) = 0;
end
	end
va = vad(1) + vad(2)*x + vad(3)*y + vad(4)*x^2 + vad(5)*x*y + vad(6)*y^2;
vadj2(k) = vpa(va);
   else
%Both sixy are conic and we locate nine points for a cubic adjacent factor.
cumat = zeros(10);
rhs  = [1;0;0;0;0;0;0;0;0;0];
xc = vert(1,k); yc = vert(2,k);
% This is the normalization node.
cum = [1 xc yc xc^2 xc*yc yc^2];
cump = [xc^3 xc^2*yc xc*yc^2 yc^3];
cumat(1,:) = [cum cump];
xc = si12(1,k); yc = si12(2,k);
cum = [1 xc yc xc^2 xc*yc yc^2];
cump = [xc^3 xc^2*yc xc*yc^2 yc^3];
cumat(2,:) = [cum cump];
xc = si14(1,k); yc = si14(2,k);
cum = [1 xc yc xc^2 xc*yc yc^2];
cump = [xc^3 xc^2*yc xc*yc^2 yc^3];
cumat(3,:) = [cum cump];
xc = si34(1,k); yc = si34(2,k);
cum = [1 xc yc xc^2 xc*yc yc^2];
cump = [xc^3 xc^2*yc xc*yc^2 yc^3];
cumat(4,:) = [cum cump];
xc = si12(1,km1); yc = si12(2,km1);
cum = [1 xc yc xc^2 xc*yc yc^2];
cump = [xc^3 xc^2*yc xc*yc^2 yc^3];
cumat(5,:) = [cum cump];
xc = si14(1,km1); yc = si14(2,km1);
cum = [1 xc yc xc^2 xc*yc yc^2];
cump = [xc^3 xc^2*yc xc*yc^2 yc^3];
cumat(6,:) = [cum cump];
xc = si34(1,km1); yc = si34(2,km1);
cum = [1 xc yc xc^2 xc*yc yc^2];
cump = [xc^3 xc^2*yc xc*yc^2 yc^3];
cumat(7,:) = [cum cump];
xc = eip(2,3*(k-1)+1); yc = eip(3,3*(k-1)+1); wc = eip(1,3*(k-1)+1);
cum = [1 xc yc xc^2 xc*yc yc^2];
cum = [1 xc yc xc^2 xc*yc yc^2];
cump = [xc^3 xc^2*yc xc*yc^2 yc^3];
	if wc == 1
cumat(8,:) = [cum cump];
	else
cumat(8,:) = [0 0 0 0 0 0 cump];
	end
wc = eip(1,3*(k-1)+2); xc = eip(2,3*(k-1)+2); yc = eip(3,3*(k-1)+2);
cum = [1 xc yc xc^2 xc*yc yc^2];
cump = [xc^3 xc^2*yc xc*yc^2 yc^3];
		if wc == 1
cumat(9,:) = [cum cump];
		else
cumat(9,:) = [0 0 0 0 0 0 cump];
        end
xc = eip(2,3*(k-1)+3); yc = eip(3,3*(k-1)+3);wc = eip(1,3*(k-1)+3);
cum = [1 xc yc xc^2 xc*yc yc^2];
cump = [xc^3 xc^2*yc xc*yc^2 yc^3];
		if wc == 1
cumat(10,:) = [cum cump];
		else
cumat(10,:) = [0 0 0 0 0 0 cump];
        end
        if norm(cumat(8,:) - cumat(9,:)) < 1e-6 | norm(cumat(10,:) - cumat(9,:)) < 1e-6          
%The adjoint must be tangent to the side curves at the double eip
%since the side curves have a common tangent here.
wc = eip(1,3*(k-1)+2); xc = eip(2,3*(k-1)+2); yc = eip(3,3*(k-1)+2);
    ns = diff(sixy(k),x);
    nsp = subs(ns, {'x' 'y'} , {xc yc}); 
    ds = diff(sixy(k),y);
    dsp = subs(ds, {'x' 'y'} , {xc yc}); 
cumat(9,:) = double(dsp*[0 wc 0 2*xc*wc wc*yc 0 3*xc^2 2*xc*yc yc^2 0] - nsp*[0 0 wc 0 wc*xc 2*wc*yc 0 xc^2 2*xc*yc 3*yc^2]);
        elseif norm(cumat(10,:) - cumat(8,:)) < 1e-6
wc = eip(1,3*(k-1)+3); xc = eip(2,3*(k-1)+3); yc = eip(3,3*(k-1)+3);
    ns = diff(sixy(k),x);
    nsp = subs(ns, {'x' 'y'} , {xc yc}); 
    ds = diff(sixy(k),y);
    dsp = subs(ds, {'x' 'y'} , {xc yc}); 
cumat(10,:) = double(dsp*[0 wc 0 2*xc*wc wc*yc 0 3*xc^2 2*xc*yc yc^2 0] - nsp*[0 0 wc 0 wc*xc 2*wc*yc 0 xc^2 2*xc*yc 3*yc^2]);
        end        
        if cumat(9,:) == conj(cumat(8,:))
            cumat(8,:) = imag(cumat(9,:));
            cumat(9,:) = real(cumat(9,:));
        elseif cumat(9,:) == conj(cumat(10,:))
            cumat(10,:) = imag(cumat(9,:));
            cumat(9,:) = real(cumat(9,:));
        end
vad = cumat\rhs;
nc = norm(vad,1);
	for km = 1:10
if abs(vad(km)) < 1e-5*nc
	vad(km) = 0;
end
    end
va = [1 x y x^2 x*y y^2 x^3 x^2*y x*y^2 y^3]*vad;
vadj2(k) = vpa(va);
   end % On types of sixy at k
                            end % On k-loop for vertices
%The degree two basis functions are now computed.                            
for k = 1:n
    km1 = k - 1;
    km2 = k - 2;
    kp1 = k + 1;
    if k == 1
        km1 = n;
        km2 = n-1;
    elseif k == 2
        km2 = n;
    end 
    wop(k) = sixy(km1);
    for rr = 1:n
         if rr ~= k & rr ~= km1
            wop(k) = wop(k)*sixy(rr);
        end
    end
%The side node opposite factors are wop.
    w214(k) = adj14(k)*wop(k);
    w212(k) = adj12(k)*wop(k);
    w234(k) = adj34(k)*wop(k);
%The side node numerators are w214, w212, and w234.
    w2vert(k) = sixy(km2);
    for rr = 1:n
        if rr ~= k & rr ~= km1 & rr ~= km2
            w2vert(k) = w2vert(k)*sixy(rr);
        end
    end
    w2vert(k) = vadj2(k)*w2vert(k);
%The vertex numerators are w2vert.    
    nv = subs(w2vert(k), {'x' 'y'}, {vert(1,k) vert(2,k)});
    dv = subs(Q1 , {'x' 'y'}, {vert(1,k) vert(2,k)});
    numW2(k) = dv*w2vert(k)/nv;
    W2v(k) = numW2(k)/Q1;
    ns = subs(w212(k), {'x' 'y'} , {si12(1,k) si12(2,k)});
    ds = subs(Q1 , {'x' 'y'}, {si12(1,k) si12(2,k)});
    num12(k) = ds*w212(k)/ns;
    W212(k) = num12(k)/Q1;;
    if degsi(k) == 1
        W214(k) = 0;
        W234(k) = 0;
    else
    ns = subs(w214(k), {'x' 'y'} , {si14(1,k) si14(2,k)});
    ds = subs(Q1 , {'x' 'y'}, {si14(1,k) si14(2,k)});
    num14(k) = ds*w214(k)/ns;
    W214(k) = num14(k)/Q1;;
    ns = subs(w234(k), {'x' 'y'} , {si34(1,k) si34(2,k)});
    ds = subs(Q1 , {'x' 'y'}, {si34(1,k) si34(2,k)});
    num34(k) = ds*w234(k)/ns;
    W234(k) = num34(k)/Q1;
    end
end
% Num are numerators and W2 are degree 2 basis functions:
%W2v(k) are vertices and W214(k),W212(k) and W234(k) are side node bases on side k.
for k = 1:n
    W2v(k) = vpa(W2v(k));
    W214(k) = vpa(W214(k));
    W212(k) = vpa(W212(k));
    W234(k) = vpa(W234(k));
end
W2v = W2v(1:n);
W214 = W214(1:n);
W212 = W212(1:n);
W234 = W234(1:n);
return