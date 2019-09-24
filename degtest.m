%degtest tests bases for degree accuracy.
syms x y
xav = sum(vert(1,:))/n;
yav = sum(vert(2,:))/n;
linf1 = 1;
quaf1 = x^2; 
linf = input('linf');
quaf = input('quaf');
%The degree one (barycentric basis is tested.
sl12w = zeros(1,n);
sq12w = zeros(1,n);
sl14w = zeros(1,n);
sq14w = zeros(1,n);
sl34w = zeros(1,n);
sq34w = zeros(1,n);
vlinw = zeros(1,n);
vquaw = zeros(1,n);
sl12 = zeros(1,n);
sl2qua = zeros(1,n);
verlin = zeros(1,n);
verqua = zeros(1,n);
w12    = zeros(1,n);
    wver = subs(wdg(1:n), {'x' 'y'} , {xav yav});
            for k = 1:n
    verlin(k) = double(subs(linf, {'x' 'y'}, {vert(1,k) vert(2,k)}));    
    verqua(k) = double(subs(quaf, {'x' 'y'}, {vert(1,k) vert(2,k)}));
            end
    vlinw = double(verlin.*wver);    
    vquaw = double(verqua.*wver);
        for k = 1:n
            if degsi(k) == 2
                w12(k) = subs(wdg(k+n), {'x' 'y'} , {xav yav});
                sl12(k) = double(subs(linf, {'x' 'y'}, {sino(1,k) sino(2,k)}));    
                sl2qua(k) = double(subs(quaf, {'x' 'y'}, {sino(1,k) sino(2,k)}));
            end
        end
        sl12w = double(sl12.*w12);    
        sq12w = double(sl2qua.*w12);
linone = sum(vlinw) + sum(sl12w);
quaone = sum(vquaw) + sum(sq12w);
linav = double(subs(linf, {'x' 'y'} , {xav yav}));
quafav = double(subs(quaf, {'x' 'y'} , {xav yav}));
linerr = linone - linav;
quaerr = quaone - quafav;
linerr = vpa(linerr,8);
quaerr - vpa(quaerr,8);
linone = vpa(linone,8);
quaone = vpa(quaone,8);
    for k = 1:n
    verlin(k) = double(subs(linf, {'x' 'y'}, {vert(1,k) vert(2,k)}));    
    verqua(k) = double(subs(quaf, {'x' 'y'}, {vert(1,k) vert(2,k)}));
    end
    verav = subs(W2v(1:n), {'x' 'y'} , {xav yav});
    vlinw = double(verlin.*verav);    
    vquaw = double(verqua.*verav);
    s12av = subs(W212(1:n), {'x' 'y'} , {xav yav});
    s14av = subs(W214(1:n), {'x' 'y'} , {xav yav});
    s34av = subs(W234(1:n), {'x' 'y'} , {xav yav});
    ss14 = zeros(1,n);
    ss34 = zeros(1,n);    
    ssq14 = zeros(1,n);
    ssq34 = zeros(1,n);    
    for k = 1:n
    ssl(k) = double(subs(linf, {'x' 'y'} , {si12(1,k) si12(2,k)}));
    ssq(k) = double(subs(quaf, {'x' 'y'} , {si12(1,k) si12(2,k)}));
    end
    wl12 = double(s12av.*ssl(1:n));
    wq12 = double(s12av.*ssq(1:n));
    for k = 1:n
        if degsi(k) == 1
            wl14 = zeros(1,n);    
            wl34 = zeros(1,n);
            wq14 = zeros(1,n);
            wq34 = zeros(1,n);
        else
            ss14(k) = double(subs(linf, {'x' 'y'}, {si14(1,k) si14(2,k)}));   
            ssq14(k) = double(subs(quaf, {'x' 'y'}, {si14(1,k) si14(2,k)}));
            ss34(k) = double(subs(linf, {'x' 'y'}, {si34(1,k) si34(2,k)}));   
            ssq34(k) = double(subs(quaf, {'x' 'y'}, {si34(1,k) si34(2,k)}));
        end
    end
            wl14 = double(s14av.*ss14(1:n));    
            wq14 = double(s14av.*ssq14(1:n));
            wl34 = double(s34av.*ss34(1:n));    
            wq34 = double(s34av.*ssq34(1:n));
lintwo = sum(vlinw) + sum(wl12) + sum(wl14) + sum(wl34);
quatwo = sum(vquaw) + sum(wq12) + sum(wq14) + sum(wq34);
linterr = lintwo - linav;
quaterr = quatwo - quafav;
linterr = vpa(linterr,8);
quaterr - vpa(quaterr,8);
lintwo = vpa(lintwo,8);
quatwo = vpa(quatwo,8);
linerr
quaerr
linterr
quaterr
return