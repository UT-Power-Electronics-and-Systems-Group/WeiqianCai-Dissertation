clear all; 
clc; close all;

%% plot the vector field
% x:theta_ba, y:theta_ca
[x,y] = meshgrid(-2*pi : 0.05 : 2*pi, -2*pi : 0.05 : 2*pi);

z1 = 2*sin(x) + sin(y) + sin(x-y);
z2 = 2*sin(y) + sin(x) + sin(y-x);

U1  = 1 - cos(x-2*pi/3) + 1 - cos(y-4*pi/3);
dU1 = sin(x-2*pi/3).*z1 + sin(y-4*pi/3).*z2;

U2  = 1 - cos(x-4*pi/3) + 1 - cos(y-2*pi/3);
dU2 = sin(x-4*pi/3).*z1 + sin(y-2*pi/3).*z2;

U3  = 1 - cos(x-pi) + 1 - cos(y-pi);
dU3 = sin(x-pi).*z1 + sin(y-pi).*z2;

U4  = 1 - cos(x-pi) + 1 - cos(y);
dU4 = sin(x-pi).*z1 + sin(y).*z2;

U5  = 1 - cos(x) + 1 - cos(y-pi);
dU5 = sin(x).*z1 + sin(y-pi).*z2;

U6  = 1 - cos(x) + 1 - cos(y);
dU6 = sin(x).*z1 + sin(y).*z2;
%%
% (0,0):  ; (0,pi):; (pi,0): (pi,pi);
U7= cos(x-pi)-cos(y);
dU7 = -sin(x-pi).*z1 +sin(y).*z2;
% h = mesh(x,y,dU7,'EdgeColor','r');
% hold on
mesh(x,y,U7,'EdgeColor','b');
text(pi,0,'$(\pi,0)$','Interpreter','latex');
hold on
mesh(x,y,0*U7,'EdgeColor','r')
% Set properties
% set(h,'linewidth',0.05);
cb = colorbar; 
gcaP=get(gca,'position');
cbP=get(cb,'Position');
cbP(3)=cbP(3)/2; % Change the colorbar width
set(cb,'Position',cbP)
set(gca,'position',gcaP)
set(gcf,'position',[300,300,600,500]); % Set the map size
caxis([-2*pi 2*pi]);
axis equal

