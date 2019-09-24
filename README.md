# barycentrics
MATLAB programs for polycons and starcons
The November 2018 description of barycentric construction for polygons and polycons is in nov25.pdf. Figures for this pdf file are in barfigs1 and barfigs2.   The general MATLAB program is poly2018  and the program for star polygons is starcon.  After running poly2018 one can run polprint17 to display results.
Coldcase6 tex, pdf, and IMG pdf 1-8 are a paper added in May 2019 in which polycon adjoint positivity was proved.
deg219 was added on July 31.  It computes degree two bases after a poly2018 run.  Updated poly2018, and deg219 were added on Sept 29.
degtest was also added.  One may verify approximation accuracy for specified functions linf and quaf.  Both degree-one and degree=two bases are used.  A quadratic is approximated well be deg219 but not by poly2018 which is only degree one.
