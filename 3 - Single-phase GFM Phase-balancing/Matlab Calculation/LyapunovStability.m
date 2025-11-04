%%
syms x;
syms y;
syms a;
syms b;
f = (y-a/3*y^3-x+a/3*x^3)    * ( sqrt(3)/2 * (1-0.5*(y-x)^2+b/4*(y-x)^4) - 0.5 * (y-x-a/3*(y-x)^3) )...
  + (y-a/3*y^3+2*x-2*a/3*x^3)* ( sqrt(3)/2 * (1-0.5*x^2+b/4*x^4)         - 0.5 * (x-a/3*x^3)       )...
  - (2*y-2*a/3*y^3+x-a/3*x^3)* ( sqrt(3)/2 * (1-0.5*y^2+b/4*y^4)         + 0.5 * (y-a/3*y^3)       );
expand(f)

%%
syms x y 
dx = 2*sin(x+2*pi/3) + sin(y+4*pi/3) + sin(x-y-2*pi/3);
dy = 2*sin(y+4*pi/3) + sin(x+2*pi/3) + sin(y-x+2*pi/3);
T1 = dx * sin(x) + dy * sin(y);
dT1 = taylor(T1,[x,y],'Order',4)

%%
syms x y
dx2 = cos(y+4*pi/3) - cos(x-y-2*pi/3);
dy2 = cos(x+2*pi/3) - cos(y-x+2*pi/3);
T2 = dx2 * sin(x) + dy2 * sin(y);
dT2 = taylor(T2,[x,y],'Order',4)