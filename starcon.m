%starcon transforms a starpolygon into a polycon and
%constructs the polycon "Wachspress Coordinates".
%GADJ finds the common divisor from the adjacent factors.
%12/2/2016
display('x and y are vertex coordinates counter-clockwise')
display('with a = 0 for concave and 1 for convex vertices')
nodes = input('[a...; x...; y...]');
m = length(nodes);
%m is the order of the polygon
cavex = nodes(1,:);
n = sum(cavex);
vert = zeros(2,n);
sino = zeros(2,n);
vin = 1;
syms x y Padj nmwdg Qeip Qofk
    for k = 1:m
        kp1 = k+1;
        kp2 = k + 2;
        if k == m
            kp1 = 1;
            kp2 = 2;
        elseif k == 1
            km1 = m;
        elseif k == m - 1
            kp2 = 1;
        end
        if nodes(1,k) == 1
           vert(:,vin) = nodes(2:3,k);
            if nodes(1,kp1) == 1
                sino(1,vin) = .5*(nodes(2,k)+nodes(2,kp1));
                sino(2,vin) = .5*(nodes(3,k)+nodes(3,kp1));
                xk = nodes(2,k); yk = nodes(3,k);
%linear sides are generated. 
                xkp1 = nodes(2,kp1); ykp1 = nodes(3,kp1);
                sino(1,vin) = .5*(nodes(2,k)+nodes(2,kp1));
                sino(2,vin) = .5*(nodes(3,k)+nodes(3,kp1));
               xp2 = nodes(2,kp2); yp2 = nodes(3,kp2);
                x1 = xk; y1 = yk; x2 = xkp1; y2 = ykp1;                        
                sixy(vin) = x1*y2-y1*x2 + (y1-y2)*x + (x2-x1)*y;
                s4 = subs(sixy(vin),[x y],[xp2 yp2]);
                s4 = double(s4);
                sixy(vin) = sign(s4)*sixy(vin);
                degsi(vin) = 1;
            else                
%A parabola is passed through k, kp1 and kp2
xv(1) = nodes(2,k); yv(1) = nodes(3,k);
xv(3) = nodes(2,kp2); yv(3) = nodes(3,kp2);
xv(2) = nodes(2,kp1); yv(2) = nodes(3,kp1);
%Concave vertices are replaced with parabola side nodes
%where vertices are ordered counter-clockwise.      LINE 50
%Some star convex vertices become intersections of adjacent parabolas.
%These may differ slightly from star polygon vertices at adjacent V's.
l1sq = (xv(1) - xv(2))^2 + (yv(1)- yv(2))^2;
l1 = sqrt(l1sq);
l3sq = (xv(2) - xv(3))^2 + (yv(2) - yv(3))^2;
l3 = sqrt(l3sq);
dot = (xv(1) - xv(2))*(xv(3)-xv(2)) + (yv(1) - yv(2))*(yv(3)-yv(2));
dang = acos(dot/(l1*l3));
                if l1 == l3
            betd = dang/2;
                else
            be = ones(1,38);
            l1dl3 = ones(1,38);
            b0 = max(0,dang-pi/2);
    for t = 2:39
        be(t-1) = b0 +.025*t*(dang-b0);
        num = 1/cos(be(t-1)) - cos(be(t-1));
        alp = dang - be(t-1);
        den = 1/cos(alp) - cos(alp);
        l1dl3(t-1) = num/den;
    end
lrat = min(l1/l3, l3/l1);
betd = interp1(l1dl3,be,lrat);
    b1 = betd - .0125*dang;
    b2 = betd + .0125*dang;
    bet = ones(1,38);
    l1dl3 = ones(1,38);
    for s = 2:39
        bet(s-1) = b1 + .025*s*(b2-b1);
        num = 1/cos(bet(s-1)) - cos(bet(s-1));
        alp = dang - bet(s-1);
        den = 1/cos(alp) - cos(alp);
        l1dl3(s-1) = num/den;
    end
betd = interp1(l1dl3,bet,lrat);
                end
betr = betd;
% betr is the angle between the larger side and the v-axis.
if l1/l3 < 1
    betr = dang - betd;
end
%betr is the angle between (1,2) and the v-axis.
if yv(1)-yv(2) == 0 
    if xv(1) - xv(2) < 0
        ang1 = pi;
    else
        ang1 = 0;
    end
elseif xv(1) - xv(2) == 0
    if yv(1) - yv(2) > 0    %LINE 100
        ang1 = pi/2;
    else
        ang1 = 3*pi/2;
    end
else
    na = yv(1)-yv(2);
    da = xv(1)-xv(2);
    ang1 = atan(na/da);
        if na > 0 & da < 0
                ang1 = pi + ang1;
        elseif na < 0 & da < 0
                ang1 = pi + ang1;
        end        
end
ssq = (sin(betr))^2;
cpar = cos(betr)/(l1*ssq);
gam = pi/2 - betr - ang1;
syms xp yp x y
up = xp*cos(gam) - yp*sin(gam);
vp = xp*sin(gam) + yp*cos(gam);
vp = subs(vp,[xp yp],[x-xv(2) y-yv(2)]);
up = subs(up,[xp yp],[x-xv(2) y-yv(2)]);
parv = vp - cpar*up^2;
pa = vpa(parv); 
%err1 = subs(pa,[x y],[xv(1) yv(1)])
%err3 = subs(pa,[x y],[xv(3) yv(3)])
%err2 = subs(pa,[x y],[xv(2) yv(2)])
sino(1,vin) = xv(2);
sino(2,vin) = yv(2);
p1 = .5*(xv(1)+xv(3));
p2 = .5*(yv(1)+yv(3));
s4 = subs(pa,[x y],[p1 p2]);
sixy(vin) = -sign(s4)*pa;
degsi(vin) = 2;
            end
vin = vin + 1;    
                if vin == n+1
                    vin = 1;
                end
        end
    end
degsi = degsi(1:n);    
%Parabolic sides have replaced concave Vs.
sixy = sixy(1:n);
%Very small components of sixy are deleted.
R = sixy;
%Small terms due to roundoff error are eliminated from R.
M = length(R);
                    for ss = 1:M
        Q = R(ss);Q = simplify(Q); %LINE 150
    cQ = coeffs(Q);         
                if cQ == Q
                    Qred(ss) = Q;
                else
[cQ,tQ] = coeffs(Q);
mm = length(cQ);
    for s = 1:mm
        c1 = coeffs(cQ(s));
        if  length(c1) == 1 & c1 == cQ(s);
                c1= double(c1);
                if abs(c1) < 5e-6
                cQ(s) = 0*cQ(s);
                end
        else
            [ccQ,tcQ] = coeffs(cQ(s));
            ccQ = double(ccQ);
            kk = length(ccQ);
                   for k = 1:kk
                        if abs(ccQ(k)) < 5e-6
                             ccQ(k) = 0;
                        end
                   end
            cQ(s) = sum(ccQ.*tcQ);
        end
    end
                end
                    Qred(ss) = sum(cQ.*tQ);
                    end
    Qred = Qred(1:n);
%Qred is the reduced R.    
sixy = vpa(Qred);
eip = zeros(3,3*n);
    for k = 1:n
        kp1 = k+1;
        km1 = k - 1;
        if k == n
            kp1 = 1;
        elseif k == 1
            km1 = n;
        end
        Padj(kp1) = 1;
        if degsi(k) + degsi(kp1) ~= 2
            v1 = vert(1,kp1); v2 = vert(2,kp1);
            rtk = solve(sixy(k),sixy(kp1));
            xk = double(rtk.x); yk = double(rtk.y);
            lr = length(xk);
            if degsi(k) + degsi(kp1) == 3
%Padj is linear
                if lr == 2
                    e = (xk(1) - v1)^2 + (yk(1) - v2)^2;  %LINE 200
                    if abs(e) < .001            
                        vert(1,kp1) = xk(1); 
                        vert(2,kp1) = yk(1);
                        eip(1,3*kp1-2) = 1;
                        eip(2,3*kp1-2) = xk(2);
                        eip(3,3*kp1-2) = yk(2);
                    else
                        vert(1,kp1) = xk(2); 
                        vert(2,kp1) = yk(2);
                        eip(1,3*kp1-2) = 1;
                        eip(2,3*kp1-2) = xk(1);
                        eip(3,3*kp1-2) = yk(1);
                    end
                else
%The eip is on the absolute line.                            
                    vert(1,kp1) = xk; 
                    vert(2,kp1) = yk;
                    eip(2,3*kp1-2) = double(diff(up,y));
                    eip(3,3*kp1-2) = -double(diff(up,x));
                    eip(1,3*kp1-2) = 0;
                end
                pt = zeros(3,1);
                if degsi(k) == 2
                    pt(2:3,1) = sino(:,k);
                else
                    pt(2:3,1) = sino(:,kp1);
                end
                pt(:,2) = eip(1:3,3*kp1-2);
                if pt(1,2) == 1
    Padj(kp1) = (pt(3,2) - pt(3,1))*x + (pt(2,1) - pt(2,2))*y + pt(3,1)*pt(2,2) - pt(2,1)*pt(3,2);
                else
    Padj(kp1) = pt(2,2)*(y - pt(3,1)) - pt(3,2)*(x - pt(2,1));
                end
                nrp = subs(Padj(kp1),[x y],[v1 v2]);
                Padj(kp1) = Padj(kp1)/nrp;
                elseif degsi(k) + degsi(kp1) == 4
%Padj is quadratic since both sides are parabolic.
                    if lr == 4
                        ep = 0; vpp = 0;
                        for s = 1:4
                            if vpp == 0
                                e = (xk(s) - v1)^2 + (yk(s) - v2)^2;
                                if abs(e) < .001
                                    vert(1,kp1) = xk(s); 
                                    vert(2,kp1) = yk(s);
                                    vpp = 1;
                                else
                                    eip(1,3*kp1-2+ep) = 1;
                                    eip(2,3*kp1-2+ep) = xk(s);
                                    eip(3,3*kp1-2+ep) = yk(s); %LINE 250
                                    ep = ep+1;
                                end                                                                                                   
                            else
                                eip(1,3*kp1-2+ep) = 1;
                                eip(2,3*kp1-2+ep) = xk(s);
                                eip(3,3*kp1-2+ep) = yk(s);
                                ep = ep+1;
                            end
                        end
                        pt(1:2,1) = sino(:,k); pt(1:2,2) = sino(:,kp1);
                        pt(1:2,3) = eip(2:3,3*kp1-2);
                        pt(1:2,4) = eip(2:3,3*kp1-1);
                        pt(1:2,5) = eip(2:3,3*kp1);                        
                        Madj = ones(6,6);
                        v1 = vert(1,kp1); v2 = vert(2,kp1);
                            for s = 1:5
Madj(s,2:6) = [pt(1,s) pt(2,s) pt(1,s)^2 pt(1,s)*pt(2,s) pt(2,s)^2];
                            end
                        Madj(6,:) = [1 v1 v2 v1^2 v1*v2 v2^2];
                        rhsM = [0;0;0;0;0;1];
                        cfs = Madj\rhsM;
                        cfs = real(cfs);
                        Padj(kp1) = [1 x y x^2 x*y y^2]*cfs;
                    elseif lr == 1
%The adjacent parabolas intersect only at the vertex and a 
%triple point on the absolute line.
%The adjacent factor is the side parabola translated to contain the side nodes.
                        vert(1,kp1) = xk;
                        vert(2,kp1) = yk;
                        eip(1,3*kp1-1) = 0;
                        eip(2,3*kp1-1) = double(diff(up,y));
                        eip(3,3*kp1-1) = double(diff(up,x));
                        eip(1,3*kp1) = 0;
                        eip(2,3*kp1) = double(diff(up,y));
                        eip(3,3*kp1) = double(diff(up,x));
                        eip(1,3*kp1-2) = 0;
                        eip(2,3*kp1-2) = double(diff(up,y));
                        eip(3,3*kp1-2) = double(diff(up,x));
                        xs = sino(1,kp1); ys = sino(2,kp1);
                        xc = sino(1,k); yc = sino(2,k);
                        u1 = subs(up,[x y],[xc yc]);
                        u1 = double(u1);
                        u2 = subs(up,[x y],[xs ys]);
                        u2 = double(u2);
                        v1 = subs(vp,[x y],[xc yc]);
                        v1 = double(v1);
                        v2 = subs(vp,[x y],[xs ys]);
                        v2 = double(v2);
                        br = .5*(u2+u1 + (v1-v2)/(cpar*(u2-u1)));
                        ar = v1 - cpar*(u1-br)^2;     %LINE 300
                        Padj(kp1) = vp - ar - cpar*(up-br)^2;           
                        nrp = subs(Padj(kp1),[x y],[xk yk]);
                        Padj(kp1) = Padj(kp1)/nrp;
                        Padj(kp1) = subs(Padj(kp1),[x y],[x y]);
                    elseif lr == 3
                        ep = 0; vpp = 0;
                        for s = 1:3
                            if vpp == 0
                                e = (xk(s) - v1)^2 + (yk(s) - v2)^2;
                                if abs(e) < .001
                                    vert(1,kp1) = xk(s); 
                                    vert(2,kp1) = yk(s);
                                    vpp = 1;
                                else
                                    eip(1,3*kp1-2+ep) = 1;
                                    eip(2,3*kp1-2+ep) = xk(s);
                                    eip(3,3*kp1-2+ep) = yk(s);
                                    ep = ep+1;
                                end
                            else
                                eip(1,3*kp1-2+ep) = 1;
                                eip(2,3*kp1-2+ep) = xk(s);
                                eip(3,3*kp1-2+ep) = yk(s);
                                ep = ep+1;
                            end
                        end
                        if ep + vpp ~= 3
                            error('lr=3 intersections not found')
                            return
                        else
                            dux = double(diff(up,x));
                            duy = double(diff(up,y));
                            eip(1,3*kp1) = 0;
                            eip(2,3*kp1) = dux;
                            eip(3,3*kp1) = duy;
                        end
                        w = dux*x + duy*y;
                        pt(1:2,1) = sino(:,k); pt(1:2,2) = sino(:,kp1);
                        pt(1:2,3) = eip(2:3,3*kp1-2);
                        pt(1:2,4) = eip(2:3,3*kp1-1);
                        Madj = ones(5,5);
                        v1 = vert(1,kp1); v2 = vert(2,kp1);
                        for s = 1:4
                            ws = double(subs(w,[x y],[pt(1,s) pt(2,s)]));
                            wsx = ws*pt(1,s);
                            wsy = ws*pt(2,s);
                            Madj(s,2:5) = [pt(1,s) pt(2,s) wsx wsy];
                        end
                        ws = double(subs(w,[x y],[v1 v2]));
                        wsx = ws*v1;            %LINE 350
                        wsy = ws*v2;            
                        Madj(5,:) = [1 v1 v2 v1*wsx v2*wsy];
                        rhsM = [0;0;0;0;1];
                        cfs = Madj\rhsM;
                        cfs = real(cfs);
                        Padj(kp1) = [1 x y x*w y*w]*cfs;
                    elseif lr == 2
%Two of the eip are on the absolute line.
                        e = (xk(1) - v1)^2 + (yk(1) - v2)^2;
                        if abs(e) < .001
                            vert(1,kp1) = xk(1); 
                            vert(2,kp1) = yk(1);
                            eip(1,3*kp1-2) = 1;
                            eip(2,3*kp1-2) = xk(2);
                            eip(3,3*kp1-2) = yk(2);
                        else
                            vert(1,kp1) = xk(2); 
                            vert(2,kp1) = yk(2);
                            eip(1,3*kp1-2) = 1;
                            eip(2,3*kp1-2) = xk(1);
                            eip(3,3*kp1-2) = yk(1);
                        end
%The other eip are on the absolute line.
                        upx = double(diff(up,x));
                        upy = double(diff(up,y));
                        eip(1,3*kp1-1) = 0;
                        eip(2,3*kp1-1) = upx;
                        eip(3,3*kp1-1) = upy;
                        eip(1,3*kp1) = 0;
                        eip(2,3*kp1) = upx;
                        eip(3,3*kp1) = upy;
                        w = x*upx+y*upy;
                        wsq = w^2;
                        pt(1:2,1) = sino(:,k); pt(1:2,2) = sino(:,kp1);
                        pt(1:2,3) = eip(2:3,3*kp1-2);
                        Madj = ones(4,4);
                        for s = 1:3
                            ws = double(subs(wsq,[x y],[pt(1,s) pt(2,s)]));
                            Madj(2:4,s) = [pt(1,s) pt(2,s) ws];
                        end
                        wkp1 = double(subs(wsq,[x y],[v1 v2]));
                        Madj(2:4,4) = [v1 v2 wkp1];
                        rhsM = [0;0;0;1];
                        cfs = Madj\rhsM;
                        cfs = real(cfs);
                        Padj(kp1) = [1 x y w]*cfs;
                    end           
            end
        end
    end                                     %line 400
    for j = 1:n
        if eip(2,3*j) - real(eip(2,3*j)) ~= 0 | eip(3,3*j) - real(eip(3,3*j)) ~= 0
            etemp = eip(:,3*j);
            eip(:,3*j) = eip(:,3*j-2);
            eip(:,3*j-2) = etemp;
        end
    end
    eip = eip(:,1:3*n);
Padj = Padj(1:n);
for j = 1:n    
%We now eliminate small components in Padj.                
Q = Padj(j);
cQ = coeffs(Q);
    if cQ == Q
        Padj(j) = Q;
    else
        [cQ,tQ] = coeffs(Q);
        mm = length(cQ);
        for s = 1:mm
            c1 = coeffs(cQ(s));
            if  length(c1) == 1 & c1 == cQ(s);             
                c1 = double(c1);
                if abs(c1) < 1e-8
                    cQ(s) = 0*cQ(s);
                end
            else
                [ccQ,tcQ] = coeffs(cQ(s));
                ccQ = double(ccQ);
                kk = length(ccQ);
                for k = 1:kk
                    if abs(ccQ(k)) < 1e-8
                         ccQ(k) = 0;
                    end
                end
                cQ(s) = sum(ccQ.*tcQ);
            end            
        end
        Padj(j) = sum(cQ.*tQ);
    end
end
adjac = vpa(Padj,6);
%The eip and Padj adjacent factors have been found.
%We now compute the gadj recursion parameters rek and reno.
%gadj computes coefficients rek (for vertex numerators) and
%reno (for side node numerators) with
%the Dasgupta-Wachspress recursion formulas.
reno = zeros(1,n);
rek = zeros(1,n);
rek(1) = 1.0;
	for k = 1:n                 %LINE 450

			km1 = k-1;
			kp1 = k+1;
            kp2 = k + 2;
            km2 = k - 2;
	   if k == 1
			km1 = n;
            km2 = n - 1;
       elseif k == n - 1
			kp2 = 1;
       elseif k == 2 
            km2 = n;
       elseif k == n
            kp1 = 1;
            kp2 = 2;
       end
		xk = vert(1,k);
		yk = vert(2,k);
		xkp1= vert(1,kp1);
		ykp1 = vert(2,kp1);
		xkp2 = vert(1,kp2);
		ykp2 = vert(2,kp2);
		xkm1 = vert(1,km1);
		ykm1 = vert(2,km1);
        xsi = sino(1,k);
        ysi = sino(2,k);
        sk = subs(sixy(kp1),[x y],[xk yk]);
        skp1atk = double(sk);
        sk = subs(sixy(km1),[x y],[xkp1 ykp1]);
        skm1atkp1 = double(sk);
        rep = 1.0;
        pk = subs(Padj(kp1),[x y],[xk yk]);          Pkp1atk = double(pk);
        pk = subs(Padj(k),[x y],[xkp1 ykp1]);
        Pkatkp1 = double(pk);
            if degsi(k) == 1        
               nrek = skp1atk/Pkp1atk;
               drek = skm1atkp1/Pkatkp1;
               rek(kp1) = rek(k)*nrek/drek;
            else        
%Side k is a parabola
            dsitdx = diff(sixy(k),x);
            ax = subs(dsitdx,[x y],[xsi ysi]);
            dsitdy = diff(sixy(k),y);
            ay = subs(dsitdy,[x y],[xsi ysi]);
            sitz = ax*x + ay*y;
            sitzsi = subs(sitz,[x y],[xsi ysi]);
            sit = sitz - sitzsi;
%sit is the tangent to parabola k.            
            sk = subs(sit,[x y],[vert(1,k) vert(2,k)]);
            sitatk = double(sk);                        %LINE 500
            sik = subs(sit,[x y],[vert(1,kp1) vert(2,kp1)]);
            sitatkp1 = double(sik);
            p1 = subs(Padj(kp1),[x y],[vert(1,k) vert(2,k)]);
            p1 = double(p1);
            bk = skp1atk*sitatk/p1;
            p2 = subs(Padj(k),[x y],[vert(1,kp1) vert(2,kp1)]);
            p2 = double(p2);
            bkp1 = skm1atkp1*sitatkp1/p2;
            rek(kp1) = rek(k)*bk/bkp1;
%sitdpp is the ratio of sit to p12 near sino on side k.
            sx = sino(1,k) - vert(1,k); sy = sino(2,k) - vert(2,k);
                for s = 1:4
                    t = 5-s;
                    if abs(sx) > abs(sy)
                        x1(s) = sino(1,k) + 1e-5*t*sx;
                        s1 = subs(sixy(k),x,x1(s));
                        y1c = double(solve(s1,y));
                            if length(y1c) == 2
            de1 = abs(y1c(1) - sino(2,k)); de2 = abs(y1c(2) - sino(2,k));
                                if de1 < de2
                                    y1(s) = y1c(1);
                                else
                                    y1(s) = y1c(2);
                                end
                            else
                                y1(s) = y1c;
                            end
                    else
                        y1(s) = sino(2,k) + 1e-5*t*sy;
                        s1 = subs(sixy(k),y,y1(s));
                        x1c = double(solve(s1,x));      
                            if length(x1c) == 2
             de1 = abs(x1c(1) - sino(1,k)); de2 = abs(x1c(2) - sino(1,k));
                                if de1 < de2
                                   x1(s) = x1c(1);
                                else
                                   x1(s) = x1c(2);
                                end
                            else
                                x1(s) = x1c;
                            end
                    end
            numk = double(subs(sit,[x y],[x1(s) y1(s)]));
            p1 = subs(Padj(kp1),[x y],[x1(s) y1(s)]);
            p1 = double(p1);
            p2 = subs(Padj(k),[x y],[x1(s) y1(s)]);
            p2 = double(p2);
            denk = p1*p2;
            sdpp(s) = numk/denk;
                end                                 %LINE 550
            sitdpp = 4*sdpp(4) - 6*sdpp(3) + 4*sdpp(2) - sdpp(1);
            dk = sixy(km1)*sixy(kp1);
            dks = subs(dk,[x y],[sino(1,k) sino(2,k)]);
            dksi = double(dks);
            bsi = dksi*sitdpp;
            reno(k) = bk*rek(k)/bsi;
            end            
    end
    rek = real(rek(1:n));
    reno = real(reno(1:n));
%The star basis functions are now generated.
                for rit = 1:n
                    rm1 = rit - 1;
                    rp1 = rit + 1;
                    if rit == 1
                    rm1 = n;
                end
        nmwdg(rit) = 1;
	for k = 1:n
		if k ~= rit & k ~= rm1 
            nmwdg(rit) = nmwdg(rit)*sixy(k);
        end
	end
        nmwdg(rit) = nmwdg(rit)*Padj(rit);
                end
%The vertex numerators have been found.
%We now consider the side nodes.
	for rit = 1:n
        if degsi(rit) == 1
            nmwdg(n+rit) = 0;
        else                        
            nmwdg(n+rit) = 1;
	   		for k = 1:n                
                if k ~= rit
                    nmwdg(n+rit) = nmwdg(n+rit)*sixy(k);
                end
            end
		end
    end
nmwdg = nmwdg(1:2*n);
%We now determine the denominator or adjoint.
QQ = 0;
for rit = 1:n
QQ = QQ + rek(rit)*nmwdg(rit) + reno(rit)*nmwdg(rit+n);
end
QQ = simplify(QQ);
%The following reduction is actually the program 'reduce'
%operating on R with with QQ = R on entering and R = Qred on exit. 
%Small terms due to roundoff error are now eliminated from Qadj.
Q = QQ;                                     %LINE 600
    cQ = coeffs(Q);
                if cQ == Q
                    Q1 = Q;
                else
[cQ,tQ] = coeffs(Q);
mm = length(cQ);
    for s = 1:mm
        c1 = coeffs(cQ(s));
        if  length(c1) == 1 & c1 == cQ(s);
                c1= double(c1);
                if abs(c1) < 5e-3
                cQ(s) = 0*cQ(s);
                end
        else
            [ccQ,tcQ] = coeffs(cQ(s));
            ccQ = double(ccQ);
            kk = length(ccQ);
                   for k = 1:kk
                        if abs(ccQ(k)) < 5e-3
                             ccQ(k) = 0;
                        end
                   end
            cQ(s) = sum(ccQ.*tcQ);
        end            
    end
Q1 = sum(cQ.*tQ);
Q1 = vpa(Q1);
Qadj = Q1;
                end
%Q1 = Qadj is the polycon adjoint.     
for k = 1:n                          
    nu(k) = rek(k)*nmwdg(k);
    wdg(k) = nu(k)/Q1;
    nu(n+k) = reno(k)*nmwdg(n+k);
    wdg(n+k) = nu(n+k)/Q1;
end
wdg = wdg(1:2*n);
nu = vpa(nu(1:2*n));
%nu is the numerators of the barycentric coordinates
%and wdg is the coordinates (that is nu/Q1).
return                                      %LINE 641