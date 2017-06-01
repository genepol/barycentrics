%stargen generates star polgons
np = input('Number of star points');
vert = zeros(2,2*np);
sino  = zeros(2,2*np);
rout = input('Outer radius');
rin = input('Inner radius');
ang = 2*pi/np;
cosang = cos(ang);
sinang = sin(ang);
for k = 1:np
    cout = exp(2*i*pi*k/np);
    cin = exp(2*i*pi*(k+.5)/np);
    vout = rout*[real(cout);imag(cout)];
    vin = rin*[real(cin);imag(cin)];
    x1 = vout(1); y1 - vout(2);
    x2 = vin(1); y2 = vin(2);
    nodes(1,2*k-1) = 1;
    nodes(2,2*k-1) = vout(1);
    nodes(3,2*k-1) = vout(2);
    nodes(1,2*k) = 0;
    nodes(2,2*k) = vin(1);
    nodes(3,2*k) = vin(2);    
end
nodes = nodes(:,1:2*np);
return