function y = func(x)
   a=1;
   b=5.1/(4*pi^2);
   c=5/pi;
   d=6;
   g=10;
   h=1/(8*pi);
    x1=x(1);
    x2=x(2);
   y=a*(x2-(b*x1.^2)+(c*x1)-d).^2+g*(1-h)*cos(x1)+g;
end