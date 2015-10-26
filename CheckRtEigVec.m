syms g u H
A = [0                  1           0
     (g-3)*u^2/2        (3-g)*u     (g-1)
     -H*u+(g-1)*u^3/2   H-(g-1)*u   g*u];
[P,L] = eig(A) 