function yprpr = hessianofbranin(x)
syms x1 x2
a=1;
b=5.1/(4*pi^2);
c=5/pi;
d=6;
g=10;
h=1/(8*pi);
F=a*(x2-(b*x1.^2)+(c*x1)-d).^2+g*(1-h)*cos(x1)+g;
H=hessian(F);
yprpr=double(subs(H,[x1 x2],[x(1) x(2)]));
end

