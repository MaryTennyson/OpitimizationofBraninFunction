clear all
close all
clc

X1=-5:0.01:10;
X2= 0:0.01:15;
[x1,x2]=meshgrid(X1,X2);
a=1;
b=5.1/(4*pi^2);
c=5/pi;
d=6;
g=10;
h=1/(8*pi);
F=a*(x2-(b*x1.^2)+(c*x1)-d).^2+g*(1-h)*cos(x1)+g;
realFMin = min(min(F))

mesh(x1,x2,F)
figure
contourf(x1,x2,F)
hold on

%% Newton-Raphson
fprintf('Newton-Raphson Algorithm\n');
x=[-1,1];
%%x=[9,3];
%%x=[-4;10];
epsilon=10^(-3);

tic
%%fprintf('k=1, x1=%f, x2=%f, f(x)=%f\n',x(1),x(2),func(x))
%%plot(x(1),x(2),'r.')
x_next = x-inv(hessianofbranin(x))*gradofbranin(x);
%%fprintf('k=2, x1=%f, x2=%f, f(x)=%f, abs. error=%f\n',x_next(1),x_next(2),func(x_next),abs(func(x_next)-func(x)));
%%plot(x_next(1),x_next(2),'r*');
k=1;
%while(abs(func(x_next)-func(x))>epsilon)
while(norm(gradofbranin(x_next))>epsilon)
   x=x_next;
   x_next=x-inv(hessianofbranin(x))*gradofbranin(x);
   fprintf('k=%d,x1=%f,x2=%f,f(x)=%f, error=%f\n', k , x_next(1),x_next(2),func(x_next),norm(gradofbranin(x_next)))
   plot(x_next(1),x_next(2),'r*')
   k=k+1;
end
toc
title('Newton-Raphson Algorithm')
set(gca,'fontsize',35)
set(findobj(gca, 'Type', 'Line', 'Linestyle', '--'), 'LineWidth', 2);

%% Hestenes-Stiefel Algorithm
fprintf('Hestenes-Stiefel Algorithm\n');
x=[-1,1];
tic
fprintf('k=1, x1=%f, x2=%f,f(x)=%f\n',x(1),x(2),func(x))
g= gradfunc(x);
d=-g;

% alpha argmin procedure
alpha=0:0.01:1;
funcalpha=zeros(length(alpha),1);
for i=1:length(alpha)
    funcalpha(i)=func(x+alpha(i)*d);
end
[val, ind]= min(funcalpha);
alpha = alpha(ind);
% end of alpha procedure

x_next = x+alpha*d;
g_next = gradofbranin(x_next);
beta = (g_next'*(g_next-g))/(d'*(g_next-g));
d_next = -g_next+beta*d;
fprintf('k=2, x1=%f, x2=%f,f(x)=%f, abs. error=%f\n',x_next(1),x_next(2),func(x_next),abs(func(x_next)-func(x)))

k=3;
while(norm(gradofbranin(x_next))>epsilon)
    
    x = x_next;
    g = g_next;
    d = d_next;

    % alpha argmin procedure
    alpha=0:0.01:1;
    funcalpha=zeros(length(alpha),1);
    for i=1:length(alpha)
        funcalpha(i)=func(x+alpha(i)*d);
    end
    [val, ind]= min(funcalpha);
    alpha = alpha(ind);
    % end of alpha procedure

    x_next = x+alpha*d;
    g_next = gradofbranin(x_next);
    beta = (g_next'*(g_next-g))/(d'*(g_next-g));
    d_next = -g_next+beta*d;

    fprintf('k=%d, x1=%f, x2=%f,f(x)=%f, abs. error=%f\n',k,x_next(1),x_next(2),func(x_next),abs(func(x_next)-func(x)))
    k=k+1;

    if x==x_next
        break;
    end
end
toc

%% Polak–Ribière Algorithm
fprintf('Polak–Ribière Algorithm\n');
x=[-1,1];
tic
fprintf('k=1, x1=%f, x2=%f,f(x)=%f\n',x(1),x(2),func(x))
g=gradofbranin(x);
d=-g;

% alpha argmin procedure
alpha=0:0.01:1;
funcalpha=zeros(length(alpha),1);
for i=1:length(alpha)
    funcalpha(i)=func(x+alpha(i)*d);
end
[val, ind]= min(funcalpha);
alpha = alpha(ind);
% end of alpha procedure

x_next = x+alpha*d;
g_next = gradofbranin(x_next);
beta = (g_next'*(g_next-g))/((g'*g));
d_next = -g_next+beta*d;
fprintf('k=2, x1=%f, x2=%f,f(x)=%f, abs. error=%f\n',x_next(1),x_next(2),func(x_next),abs(func(x_next)-func(x)))

k=3;
while(norm(gradofbranin(x_next))>epsilon)
    % if (norm(powellSingularGradient(x_next))>epsilon)
    %     break;
    % end
    x = x_next;
    g = g_next;
    d = d_next;

    % alpha argmin procedure
    alpha=0:0.01:1;
    funcalpha=zeros(length(alpha),1);
    for i=1:length(alpha)
        funcalpha(i)=func(x+alpha(i)*d);
    end
    [val, ind]= min(funcalpha);
    alpha = alpha(ind);
    % end of alpha procedure

    x_next = x+alpha*d;
    g_next = gradofbranin(x_next);
    beta = (g_next'*(g_next-g))/((g'*g));
    d_next = -g_next+beta*d;

    fprintf('k=%d, x1=%f, x2=%f,f(x)=%f, abs. error=%f\n',k,x_next(1),x_next(2),func(x_next),abs(func(x_next)-func(x)))
    k=k+1;
    if x==x_next
        break;
    end
end
toc

%% Flecther-Reeves Algorithm
fprintf('Flecther-Reeves Algorithm\n');
x=[-1,1];
tic
fprintf('k=1, x1=%f, x2=%f, f(x)=%f\n',x(1),x(2),func(x))
g=gradofbranin(x);
d=-g;

% alpha argmin procedure
alpha=0:0.01:1;
funcalpha=zeros(length(alpha),1);
for i=1:length(alpha)
    funcalpha(i)=func(x+alpha(i)*d);
end
[val, ind]= min(funcalpha);
alpha = alpha(ind);
% end of alpha procedure

x_next = x+alpha*d;
g_next = gradofbranin(x_next);
beta = (g_next'*g_next)/(g'*g);
d_next = -g_next+beta*d;
fprintf('k=2, x1=%f, x2=%f,f(x)=%f, abs. error=%f\n',x_next(1),x_next(2),func(x_next),abs(func(x_next)-func(x)))

k=3;
while(norm(gradofbranin(x_next))>epsilon)
  
    x = x_next;
    g = g_next;
    d = d_next;

    % alpha argmin procedure
    alpha=0:0.01:1;
    funcalpha=zeros(length(alpha),1);
    for i=1:length(alpha)
        funcalpha(i)=func(x+alpha(i)*d);
    end
    [val, ind]= min(funcalpha);
    alpha = alpha(ind);
    % end of alpha procedure

    x_next = x+alpha*d;
    g_next = gradofbranin(x_next);
    beta = (g_next'*g_next)/(g'*g);
    d_next = -g_next+beta*d;

    fprintf('k=%d, x1=%f, x2=%f,f(x)=%f, abs. error=%f\n',k,x_next(1),x_next(2),func(x_next),abs(func(x_next)-func(x)))
    k=k+1;
    if x==x_next
        break;
    end
end
toc


