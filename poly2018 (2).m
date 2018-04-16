%poly2018 generates polycon basis functions.
syms x y 
disp('Defaults with rep = 2, ... ,8')
rep = input('rep = 0 for new pol and rep = 2-8 for test pols.')
if rep == 0
n = input('The number of sides n is: ')
  for k = 1:n
   sixy(k) = input('side k');
   xk(k) = input('Estimated x at vertex k is')
   yk(k) = input('Estimated y at vertex k is')
   degsi(k) = input('The degree of side k is:')
  end
  vertest = [xk(1:n);yk(1:n)];  
elseif rep == 2
    n = 4;
    sixy(1) = x;
    sixy(2) = y;
    sixy(3) = 1 - x^2 - y^2;
    sixy(4) = 1 + x - 4*y^2;
    degsi = [1 1 2 2];
    vertest = [0 0 1 .75;.5 0 0 .6614];
elseif rep == 3
    n = 3;
    sixy(1) = 1 - x^2 - y^2;
    sixy(2) = -.5 + 3*x - 3*y - x^2 - y^2;
    sixy(3) = -.5 + 3*x + 3*y - x^2 - y^2;
    degsi = [2 2 2];
    vertest = [.9114 .9114 .177;-.4114 .4114 0];
elseif rep == 4
    n = 3;
    sixy(1) = x;
    sixy(2) = y;
    sixy(3) = 1 - x - y^2;
    degsi = [1 1 2];
    vertest = [0 0 1;1 0 0];
elseif rep == 5
    n =4;
    sixy(1) = y - x^2;
    sixy(2) = 3 - 2*x - y;
    sixy(3) = 4 - x^2 - y^2;
    sixy(4) = 3 + 2*x - y;
    degsi = [2 1 2 1];
    vertest = [-1 1 .5367 -.5367;1 1 1.9266 1.9266];  
elseif rep == 6
    n = 5;
    sides = zeros(n,6);
    sixy(1) =  16 -x^2 - y^2;
    sixy(2) = 3 - x;
    sixy(3) = 1- x*y;
    sixy(4) = 3 - y;                %Line 50
    sixy(5) = 4 + 7*x + y;	
    degsi = [2 1 2 1 1];
    vertest = [0 3 3 .3333 -1;-4 -2.65 .3333 3 3];
elseif rep == 7
    n = 4;
    sixy(1) =  8*x - x^2 - y^2;
    sixy(2) =  15 - 4*x + 12*y + x^2 + y^2;
    sixy(3) =  16 - x^2 - y^2;
    sixy(4) =  16 + 4*x - 8*y - x^2 - y^2;
    degsi = [2 2 2 2];
    vertest = [.4228 .2259 3.7741 3.5777;1.7885 -1.3253 -1.3253 1.7889];
elseif rep == 8
    n = 6; 
    sixy(1) = -.0001*x + y;
    sixy(2) =  2 - x;
    sixy(3) =  2 - y;
    sixy(4) =  2 + x;
    sixy(5) =  .0001*x + y;
    sixy(6) = -1 + x^2 + y^2;
    degsi = [1 1 1 1 1 2];
    vertest = [1 2 2 -2 -2 -1;.0001 .0002 2 2 .0002 .0001];
end
    for k = 1:n
        km1 = k-1;
        if k == 1
            km1 = n;
        end
        sxm1 = double(subs(sixy(k),[x y],[vertest(1,km1) vertest(2,km1)]));
        sig = sign(sxm1);
        if sig ~= 1
            sixy(k) = -sixy(k);
        end
    end
%Sides are positive within the element.
sixy = vpa(sixy,6);
    vert = zeros(2,n);
    eip = zeros(3,3*n);
    for k = 1:n
        kp1 = k+1;
        km1 = k - 1;
        if k == n
            kp1 = 1;             
        elseif k == 1
            km1 = n;                    
        end
        if degsi(k)*degsi(km1) == 4
            vkkm1 = solve(sixy(km1),sixy(k));
            xvk = double(vpa(vkkm1.x));
            yvk = double(vpa(vkkm1.y));
            lx = length(xvk);               %Line 100
            emin = 1e6;             
            for s = 1:lx
                e1(s) = abs(xvk(s) - vertest(1,k)) + abs(yvk(s) - vertest(2,k));
                if e1(s) < emin
                    emin = e1(s);
                    smin = s;
                end
            end
            vert(1,k) = xvk(smin);
            vert(2,k) = yvk(smin);
            ep = 0;
            for s = 1:lx
                if s ~= smin
                    eip(1,3*k-2+ep) = 1.;
                    eip(2,3*k-2+ep) = xvk(s);
                    eip(3,3*k-2+ep) = yvk(s);
                    ep = ep + 1;
                end
            end
            if ep == 1 %There are two eip on the absolute line
                ax = diff(sixy(k),x); ay = diff(sixy(k),y);
                axx = .5*double(diff(ax,x)); 
                axy = double(diff(ax,y));
                ayy = .5*double(diff(ay,y));
                disc = double(sqrt(axy^2 - 4*axx*ayy));
                if axx ~= 0
                        eip(3,3*k2) = 2*axx;
                        eip(2,3*k2) = -axy + disc;
                        eip(3,3*k2-1) = 2*axx;
                        eip(2,3*k2-1) = -axy - disc;                    
                elseif ayy ~= 0
                        eip(2,3*k2) = 2*ayy;
                        eip(3,3*k2) = -axy + disc;
                        eip(2,3*k2-1) = 2*ayy;
                        eip(2,3*k2-1) = -axy - disc;                          
                else
                    eip(3,3*k2) = 1.
                    eip(2,3*k2) = 0;
                    eip(3,3*k2-1) = 0;
                    eip(2,3*k2-1) = 1.
                end
            end
        else
            if degsi(k) == 1                                    
                dkx = diff(sixy(k),x); dky = diff(sixy(k),y);       
                if dkx == 0
                    ykkm1 = solve(sixy(k));
                    xkkm1 = subs(sixy(km1),y,ykkm1);
                    xkkm1 = solve(xkkm1);
                elseif dky == 0                             %Line 150
                    xkkm1 = solve(sixy(k));
                    ykkm1 = subs(sixy(km1),x,xkkm1);
                    ykkm1 = solve(ykkm1);
                else
                    xkofy = solve(sixy(k));
                    ykkm = subs(sixy(km1),x,xkofy);
                    ykkm1 = solve(ykkm);
                    xkkm1 = subs(xkofy,y,ykkm1);
                end
            else
                dkx = diff(sixy(km1),x); dky = diff(sixy(km1),y);
                if dkx == 0
                    ykkm1 = solve(sixy(km1));
                    xkkm1 = subs(sixy(k),y,ykkm1);
                    xkkm1 = solve(xkkm1);
                elseif dky == 0
                    xkkm1 = solve(sixy(km1));
                    ykkm1 = subs(sixy(k),x,xkkm1);
                    ykkm1 = solve(ykkm1);
                else
                    xkofy = solve(sixy(km1));
                    ykkm = subs(sixy(k),x,xkofy);
                    ykkm1 = solve(ykkm);
                    xkkm1 = subs(xkofy,y,ykkm1);
                end
            end
            if length(ykkm1)== 2 & length(xkkm1) == 1
                xkkm1 = [xkkm1;xkkm1];
            elseif length(ykkm1)== 1 & length(xkkm1) == 2
                ykkm1 = [ykkm1;ykkm1];
            end
            xkkm1 = double(xkkm1);
            ykkm1 = double(ykkm1);
            if degsi(k)*degsi(km1) == 1         
                vert(1,k) = xkkm1;
                vert(2,k) = ykkm1;
            else
                ly = length(ykkm1);
                if ly == 2
                    eip(1,3*k-2) = 1.;                          
                    e1 = (ykkm1(1) - vertest(2,k))^2 + (xkkm1(1) - vertest(1,k))^2;
                    e2 = (ykkm1(2) - vertest(2,k))^2 + (xkkm1(2) - vertest(1,k))^2;
                    if e1 < e2
                         vert(1,k) = xkkm1(1);              
                         vert(2,k) = ykkm1(1);                 
                         eip(2,3*k-2) = xkkm1(2);
                         eip(3,3*k-2) = ykkm1(2);
                    else
                         vert(1,k) = xkkm1(2);
                         vert(2,k) = ykkm1(2);              %Line 200
                         eip(2,3*k-2) = xkkm1(1);                        
                         eip(3,3*k-2) = ykkm1(1);
                    end
                else
%The eip is on the absolute line
                    vert(1,k) = xkkm1;
                    vert(2,k) = ykkm1;
                    if degsi(k) == 1
                         eip(2,3*k-2) = double(diff(sixy(k),x));                        
                         eip(3,3*k-2) = double(diff(sixy(k),y));
                    else
                        eip(2,3*k-2) = double(diff(sixy(km1),x));                        
                        eip(3,3*k-2) = double(diff(sixy(km1),y));
                    end
                end               
            end
        end
    end        
%The vertices and eip have been computed.
%The side nodes must now be determined.
sino = zeros(2,n);        
    for k = 1:n
        km1 = k-1;
        if k == 1
            km1 = n;
        end
        kp1 = k + 1;
        kp2 = k + 2;
        if k == n
            kp1 = 1;
            kp2 = 2;
        end
        if k == n -1
            kp2 = 1;                       
        end        
        x1 = (vert(1,k) + vert(1,kp1))/2;
        y1 = (vert(2,k) + vert(2,kp1))/2;
        x2 = (vert(1,kp1) + vert(1,kp2))/2;
        y2 = (vert(2,kp1) + vert(2,kp2))/2;
        if degsi(k) == 1
            sino(:,k) = [x1;y1];                 
        else  %the side is a conic
            pb1 = 2*(vert(1,kp1)-vert(1,k))*x;  
            pb2 = (vert(2,kp1)-vert(2,k))*y;                    
            pb3 = vert(1,k)^2 - vert(1,kp1)^2 + vert(2,k)^2 - vert(2,kp1)^2;  
            pbi = pb1 + pb2 + pb3;
            sin2 = solve(pbi,sixy(k));
            sinx = double(sin2.x);
            siny = double(sin2.y);
            sumsig = 0;                                         %Line 250
            for t = 1:n
                if t ~= k
                    si = subs(sixy(t),[x y],[sinx(1) siny(1)]);
                    si = double(si);
                    sig = sign(si);
                    sumsig = sumsig + sig;
                end
            end
            if sumsig == n-1
                sino(1,k) = sinx(1); sino(2,k) = siny(1);
            else
                sino(1,k) = sinx(2); sino(2,k) = siny(2);
            end
        end
        if degsi(kp1) == 1
            sino(:,kp1) = [x2;y2];
        else  %the side is a conic
            pb1 = 2*(vert(1,kp2)-vert(1,kp1))*x;  
            pb2 = (vert(2,kp2)-vert(2,kp1))*y;
            pb3 = vert(1,kp1)^2 - vert(1,kp2)^2 + vert(2,kp1)^2 - vert(2,kp2)^2;
            pbi = pb1 + pb2 + pb3;
            sin2 = solve(pbi,sixy(kp1));
            sinx = double(sin2.x);
            siny = double(sin2.y);
            sumsig = 0;
            for t = 1:n
                if t ~= kp1
                    si = subs(sixy(t),[x y],[sinx(1) siny(1)]);
                    si = double(si);
                    sig = sign(si);
                    sumsig = sumsig + sig;
                end
            end
            if sumsig == n-1                    
                sino(1,kp1) = sinx(1); sino(2,kp1) = siny(1);
            else
                sino(1,kp1) = sinx(2); sino(2,kp1) = siny(2);
            end
        end
    end
%The adjacent factors must now be found.
syms Padj
     for k = 1:n         
        km1 = k-1;                          
        if k == 1                
            km1 = n;
        end
        kp1 = k + 1;
        kp2 = k + 2;
        if k == n                                           %Line 300
            kp1 = 1;
            kp2 = 2;
        end
        if k == n -1
            kp2 = 1;
        end 
        if degsi(k) + degsi(kp1) == 2
            Padj(kp1) = 1.;
        elseif degsi(k) + degsi(kp1) == 3
            xe = eip(2,3*kp1-2); ye = eip(3,3*kp1-2);
%The adjacent factor is linear.
            if degsi(k) == 2
                xs = sino(1,k); ys = sino(2,k);
            else
                xs = sino(1,kp1); ys = sino(2,kp1);
            end
            if eip(1,3*kp1-2) == 0
%The eip is on the absolute line.
                Padj(kp1) = xe*(x - xs) + ye*(y - ys);
            else
                Padj(kp1) = (ye - ys)*x + (xs - xe)*y - xs*ye + xe*ys;
            end
            pnr = subs(Padj(kp1),[x y],[vert(1,kp1) vert(2,kp1)]);
            nrp = double(pnr);
            Padj(kp1) = Padj(kp1)/nrp;
        else      %Both sides are conic.
            esum = eip(1,3*kp1-2) + eip(1,3*kp1-1) + eip(1,3*kp1);
            if esum == 0
%Adjacent translated cnics must be addressed 
                syms uu vv
                fpar = subs(sixy(k),[x y],[x-uu y-vv]);
                gpar = subs(fpar,[x y],[sino(1,kp1) sino(2,kp1)]);
                hpar = subs(fpar,[x y],[sino(1,k) sino(2,k)]);
                uv = solve(hpar,gpar);
                if length(uv) ~= 1
                    error('to be programmed')
                end
                uvx = double(uv.uu);
                uvy = double(uv.vv);
                Padj(kp1) = subs(fpar,[uu vv],[uvx uvy]);
            else
fx = [sino(1,kp1) sino(1,k) eip(2,3*kp1-2) eip(2,3*kp1-1) eip(2,3*kp1)];
fy = [sino(2,kp1) sino(2,k) eip(3,3*kp1-2) eip(3,3*kp1-1) eip(3,3*kp1)];
mat = zeros(6);                                                                 
                for j = 1:5                                             
                    mat(j,1:6) = [1 fx(j) fy(j) fx(j)^2 fx(j)*fy(j) fy(j)^2];
                end
                for j = 1:3
                    if eip(1,3*kp1+1-j) == 0
                        mat(6-j,1:3) = [0 0 0];                    %Line 350
                    end
                end
                g1 = vert(1,kp1); g2 = vert(2,kp1);                 
                mat(6,:) = [1 g1 g2 g1^2 g1*g2 g2^2];
                rhs = [0 0 0 0 0 1]';
                if mat(3,:) == mat(4,:)
                    ex = eip(2,3*kp1-1); ey = eip(3,3*kp1-1);
                elseif mat(3,:) == mat(5,:)        
                    ex = eip(2,3*kp1); ey = eip(3,3*kp1);
                elseif mat(4,:) == mat(5,:)        
                    ex = eip(2,3*kp1); ey = eip(3,3*kp1);
                end
                if mat(3,:) == mat(4,:) | mat(3,:) == mat(5,:) | mat(4,:) == mat(5,:)
                    dsx = diff(sixy(kp1),x);
                    dsy = diff(sixy(kp1),y);
                    dysy = diff(dsy,y);
                    dxsy = diff(dsy,x);
                    dxsx = diff(dsx,x);
                    dsxe = double(subs(dsx,[x y],[ex ey]));
                    dsye = double(subs(dsy,[x y],[ex ey]));
                    dxsxe = double(subs(dxsx,[x y],[ex ey]));
                    dysye = double(subs(dysy,[x y],[ex ey]));        
                    dxsye = double(subs(dxsy,[x y],[ex ey]));        
                    matx = [0 1 0 2*ex ey 0];
                    maty = [0 0 1 0 ex 2*ey];
                    if mat(3,:) == mat(4,:)
                        mat(4,:) = dsye*matx-dsxe*maty;
                    else                   
                        mat(5,:) = dsye*matx-dsxe*maty;
                    end
                end
                cofs = real(mat\rhs);                   
                Padj(kp1) = cofs'*[1;x;y;x^2;x*y;y^2];
            end
        end
     end
     Padj = Padj(1:n);
    adjac = vpa(Padj,6);                   
%We now compute the gadj recursion parameters rek and reno.
%gadj computes coefficients rek (for vertex numerators) and
%reno (for side node numerators) with
%the Dasgupta-Wachspress recursion formulas.
reno = zeros(1,n);
rek = zeros(1,n);
rek(1) = 1.0;                           
	for k = 1:n        
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
        pk = subs(Padj(kp1),[x y],[xk yk]);     
        Pkp1atk = double(pk);                   
        pk = subs(Padj(k),[x y],[xkp1 ykp1]);       
        Pkatkp1 = double(pk);               
        if degsi(k) == 1        
           nrek = skp1atk/Pkp1atk;      
           drek = skm1atkp1/Pkatkp1;        
           rek(kp1) = rek(k)*nrek/drek;
        else                            
%Side k is a conic
            dsitdx = diff(sixy(k),x); dsitdy = diff(sixy(k),y);
            ax = subs(dsitdx,[x y],[xsi ysi]);                 
            ay = subs(dsitdy,[x y],[xsi ysi]);
            sitz = ax*x + ay*y;
            sitzsi = subs(sitz,[x y],[xsi ysi]);
            sit = sitz - sitzsi;
%sit is the tangent to conic k.            
            sk = subs(sit,[x y],[vert(1,k) vert(2,k)]);    
            sitatk = double(sk);
            sik = subs(sit,[x y],[vert(1,kp1) vert(2,kp1)]);
            sitatkp1 = double(sik);
            p1 = subs(Padj(kp1),[x y],[vert(1,k) vert(2,k)]);
            p1 = double(p1);
            bk = skp1atk*sitatk/p1;                             %Line 450
            p2 = subs(Padj(k),[x y],[vert(1,kp1) vert(2,kp1)]);   
            p2 = double(p2);
            bkp1 = skm1atkp1*sitatkp1/p2;
            rek(kp1) = rek(k)*bk/bkp1;
%We now determine the ratio of sit to p1*p2 at the side node
       if abs(double(ax)) < abs(double(ay))
            con2 = solve(sixy(k),y);
            if length(con2) == 1
                con1 = con2;
            else
                cc1 = double(subs(con2,x,xsi));
                if abs(ysi - cc1(1)) < abs(ysi - cc1(2))
                    con1 = con2(1);
                else                                
                    con1 = con2(2);
                end
            end            
            sitk = subs(sit,y,con1);
            pk = subs(Padj(k),y,con1);
            pkp1 = subs(Padj(kp1),y,con1);                      
            dsitk = diff(sitk,x);
            dsitk2 = diff(dsitk,x);
            sisitk = double(subs(dsitk2,x,xsi));
            dpk1 = diff(pk,x);
            dpk1si = double(subs(dpk1,x,xsi));
            dpkp1 = diff(pkp1,x);
            dpkp1si = double(subs(dpkp1,x,xsi));
       else
            con2 = solve(sixy(k),x);
            if length(con2) == 1
                con1 = con2;
            else
                cc1 = double(subs(con2,y,ysi));                
                if abs(xsi - cc1(1)) < abs(xsi - cc1(2))
                    con1 = con2(1);
                else
                    con1 = con2(2);
                end
            end            
            sitk = subs(sit,x,con1);
            pk = subs(Padj(k),x,con1);
            pkp1 = subs(Padj(kp1),x,con1);
            dsitk = diff(sitk,y);
            dsitk2 = diff(dsitk,y);
            sisitk = double(subs(dsitk2,y,ysi));
            dpk1 = diff(pk,y);
            dpk1si = double(subs(dpk1,y,ysi));
            dpkp1 = diff(pkp1,y);
            dpkp1si = double(subs(dpkp1,y,ysi));
       end                                                      %line 500
       sitdpp = sisitk/(2*dpk1si*dpkp1si);
            dk = sixy(km1)*sixy(kp1);
            dks = subs(dk,[x y],[sino(1,k) sino(2,k)]);
            dksi = double(dks);
            bsi = double(dksi*sitdpp);
            reno(k) = bk*rek(k)/bsi;
        end
    end    
    rek = real(rek(1:n));
    reno = real(reno(1:n));
%The basis functions are now generated.
syms nmwdg
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
for k = 1:n
    nu(k) = rek(k)*nmwdg(k);
    nu(n+k) = reno(k)*nmwdg(n+k);      
    QQ = QQ + nu(k) + nu(n+k);
end
Qadj = simplify(QQ);
    cQ = coeffs(Qadj);                                      %LINE 550
                if cQ ~= Qadj
                    [cQ,tQ] = coeffs(Qadj);
                    mm = length(cQ);
                    for s = 1:mm
                        c1 = coeffs(cQ(s));
                        if  length(c1) == 1 & c1 == cQ(s);
                            c1 = double(c1);
                            if abs(c1) < 1e-4
                                cQ(s) = 0*cQ(s);
                            end
                        else
                            [ccQ,tcQ] = coeffs(cQ(s));
                            ccQ = double(ccQ);
                            kk = length(ccQ);               
                             for k = 1:kk
                                if abs(ccQ(k)) < 1e-4
                                    ccQ(k) = 0;
                                end
                             end                        
                            cQ(s) = sum(ccQ.*tcQ);          
                        end
                    end
                    Qadj = sum(cQ.*tQ);
                end
Q1 = vpa(Qadj,6);
%Qadj is the polycon adjoint. 
syms wdg
for k = 1:2*n                          
    wdg(k) = nu(k)/Qadj;              
end
wdg = wdg(1:2*n);
nu = nu(1:2*n);
%nu is the numerators of the barycentric coordinates
%and wdg is the coordinates (that is nu/Q1).
return                                                     %Line 585