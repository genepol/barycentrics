%polprint displays key results for poly2014.
display('vertices 1:n')
vert
display('sides')
s6 = vpa(sixy,6);
pretty(s6)
display('Exterior intersections, up to 3 per vertex')
eip
display('Side node coordinates')
sino
display('GADJ vertex values 1:n plus n+1 = 1')
rek
display('GADJ side values 1:n')
reno
display('Adjacent factors')
adj6 = vpa(adjac,6);
pretty(adj6)
Q = simple(Q1);Q = simplify(Q);
Q6 = vpa(Q,6);
display('Adjoint')
pretty(Q6)
display('Numerators')
nu6 = vpa(nu,6);
for k = 1:2*n    
    pretty(nu6(k))
end
return

