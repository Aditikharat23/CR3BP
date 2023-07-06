%Main file
%REFERENCE: Celestial Mechanics Notes : The Circular Restricted Three Body Problem J.D. Mireles James

clear all
close all
%% Parameters
%scaling d'=Ld,s'=Vs,t'=Tt/2pi;V = 2piL/T;

m1 = 5.97219e24;                %Mass of Earth [kg]
m2 = 7.35e22;                   %Mass of moon
m3 = 100;                       %Mass of satellite
d = 3.84e8;                     % distance between earth and moon in [m]
G  =6.6743e-11;                 %Gravitational constant in Nm^2/kg^2
T= 2*pi*sqrt(d^3/(G*(m1+m2)));  %Time period in seconds time taken for moon to orbit earth
n = 2*pi/T;                     %Orbital velocity of moon

u = m2/(m1+m2);                 %Mass parameter
a=  1;
% r = (1.743e3 + 150)/d;
%% % Finding co-ordinates of the lagrange points in rotating frame
%Location of masses in RF
x1 = -u;
x2 = 1-u;
%EQUILATERAL POINTS
L4_x = cosd(60)-u;
L4_y = sind(60);
L5_x = L4_x;
L5_y = -L4_y;
L5_z = 0;
L4_z =0;

%Collinear 
% figure(1)
% fplot(@fn.collinear) % Visualise the function to find the range for roots
L2 = fzero(@fn.collinear,[1,1.5]);
L1 = fzero(@fn.collinear,[0.4,0.9]);%Between the two bodies
L3 = fzero(@fn.collinear,[-1.5,-0.5]);

%%Plotting the lagrange points in rotating reference
figure(1)
plot3(x1,0,0,'Marker','o','MarkerFaceColor','b','MarkerSize',10)
hold on
plot3(x2,0,0,'Marker','o','MarkerFaceColor','#7E2F8E','MarkerSize',6)
hold on
plot3(L4_x,L4_y,0,'Marker','diamond','MarkerFaceColor','g','MarkerEdgeColor','g')
hold on
plot3(L5_x,L5_y,0,'Marker','diamond','MarkerFaceColor','red','MarkerEdgeColor','red')
hold on
plot3(L3,0,0,'Marker','diamond','MarkerFaceColor','cyan','MarkerEdgeColor','cyan')
hold on
plot3(L2,0,0,'Marker','diamond','MarkerFaceColor','magenta','MarkerEdgeColor','magenta')
hold on
plot3(L1,0,0,'Marker','diamond','MarkerFaceColor','yellow','MarkerEdgeColor','yellow')
grid on
hold on
title('Lagrange points in rotating reference frame')
legend('Earth','Moon','L4','L5','L3','L2','L1')

%% Jacobi constant and zero velocity 
% basically once an energy level is fixed the magnitude  velocity is uniquely determined
% As an orbit reaches the zero velocity curve its velocity goes to 0 preventing it from crossing and constraining the motion
%These regions are allowable regions called hill's  regions

%Finding Jacobi constant and zero velocity curves for each lagrange point
%Trajectory at each lagrange point is a fixed point
C1 = fn.jacobiconst([L1,zeros(1,5)],u);
C2 = fn.jacobiconst([L2,zeros(1,5)],u);
C3  = fn.jacobiconst([L3,zeros(1,5)],u);
C4 = fn.jacobiconst([L4_x,L4_y,zeros(1,4)],u);
C5 = fn.jacobiconst([L5_x,L5_y,zeros(1,4)],u);
%% 
% %% Initial Condition
% 
% r_in = [(1-u) 0.0455 0];
% v_in =  [-0.5 0.0012 0];
% tend=10;
% C = fn.jacobiconst([r_in v_in],u);
% %% Plotting the trajectory in initial
% 
% [t,state] = ode45(@(t,state) fn.cr3bp(t,state,u), [0 20],[r_in v_in]);
% 
% figure(1)
% plot3(-u,0,0,'ro',1-u,0,0,'bo',state(:,1),state(:,2),state(:,3))

%% Plotting Zero Velocity surface for the given Initial Conditions
x = linspace(-2,2,100);
y = x;
z = x;
[X,Y,Z] = meshgrid(x,y,z);
rho1 = sqrt((X+u).^2 +Y.^2 + Z.^2);
rho2 = sqrt((X-(1-u)).^2 + Y.^2 + Z.^2);

V = X.^2 + Y.^2  + 2*(1-u)/rho1 + 2*u/rho2;

hold on
patch(isosurface(X,Y,Z,V,C2),'FaceColor','magenta','EdgeColor','none')
hold on

camlight;
alpha(0.2)
lighting  gouraud
%% %Zero velocity 2d

for i = 1:length(x)
    for j = 1:length(y)
        rho12d = sqrt((x(i)+u).^2 +y(j).^2);
rho22d = sqrt((x(i)-(1-u)).^2 + y(j).^2 );
        V2d(j,i) = x(i)^2 + y(j)^2  + 2*(1-u)/rho12d + 2*u/rho22d;
    end
end
figure(2)
subplot(1,3,1)
contourf(x,y,-V2d,[-C1 -C1])
hold on
plot(x1,0,'o',x2,0,'o',L1,0,'*',L2,0,'*',L3,0,'*',L4_x,L4_y,'*',L5_x,L5_y,"*")
title('C1 = 3.1884')
subplot(1,3,2)
contourf(x,y,-V2d,[-C2 -C2])
title('C2 = 3.1722')
hold on
plot(x1,0,'o',x2,0,'o',L1,0,'*',L2,0,'*',L3,0,'*',L4_x,L4_y,'*',L5_x,L5_y,"*")
subplot(1,3,3)
contourf(x,y,-V2d,[-C3 -C3])
hold on
plot(x1,0,'o',x2,0,'o',L1,0,'*',L2,0,'*',L3,0,'*',L4_x,L4_y,'*',L5_x,L5_y,"*")
title('C3 = 3.0122')

%% TRAJECTORY SIMULATION

 [t,stateb] = ode45(@(t,state) fn.cr3bp(t,state,u), [0 48],[-0.5,0,0,0,1,0]);
%[t,stateb] = RK4([-0.75,0,0,0,1,0],0.01,30,u);
%Body fixed frame
% plot3(stateb(:,1),stateb(:,2),stateb(:,3))


statei = zeros(length(t),3);
for i = 1:length(t)
    A = [cos(t(i)),sin(t(i)),0;...
        -sin(t(i)),cos(t(i)),0;...
        0,0,1];
    temp = A'*stateb(i,1:3)';
    statei(i,:) = temp'; %Inertial position vector for SC
end
statem = zeros(length(t),3);
for i = 1:length(t)
    m = t(i);
    A = [cos(m),sin(m),0;...
        -sin(m),cos(m),0;...
        0,0,1];
    temp = A'*[1-u;0;0];
    statem(i,:) = temp'; %Inertial position vector for SC
end
% figure(3)
% plot3(stateb(:,1),stateb(:,2),stateb(:,3))
% hold on
% plot3(x1,0,0,'Marker','o','MarkerFaceColor','b','MarkerSize',10)
% hold on
% plot3(x2,0,0,'Marker','o','MarkerFaceColor','#7E2F8E','MarkerSize',6)
% title('Body frame')
% 
% figure(4)
% plot3(statei(:,1),statei(:,2),statei(:,3))
% hold on
% plot3(statem(:,1),statem(:,2),statem(:,3),'Marker','o','MarkerFaceColor','#7E2F8E','MarkerSize',3)
% hold on
% plot3(x1,0,0,'Marker','o','MarkerFaceColor','b','MarkerSize',10)
% title('inertial')
close all
%% ANIMATING
f1 = figure(4); %body frame
f1.Position = [0,f1.Position(2:4)];%%Not sure??
clf
hold on
plot(x1,0,'Marker','o','MarkerFaceColor','b','MarkerSize',10)
hold on
plot(x2,0,'Marker','o','MarkerFaceColor','g','MarkerSize',6)
hold on 
contourf(x,y,-V2d,[-C1 -C1])
title('Body frame')
hold on 
grid on

satb = plot(stateb(1,1),stateb(1,1),'r.','MarkerSize',10);%For the marker
satb1 = plot(stateb(1,1),stateb(1,1),'r');%For the trajectory
hold off
set(gca,'Xlim',[-2,2],'YLim',[-2,2],'Xtick',[],'YTick',[]);
f2 = figure(5); %body frame
f2.Position = [0,f2.Position(2:4)];%%Not sure??
clf
hold on
plot3(x1,0,0,'Marker','o','MarkerFaceColor','b','MarkerSize',20)
title('Inertial frame')
hold on 


moon = plot(statem(1,1),statem(1,1),'g.','MarkerSize',10);%For the marker
moon1 = plot(statem(1,1),statem(1,1),'r');%For the trajectory

sat_i = plot(statei(1,1),statei(1,1),'r.','MarkerSize',10);%For the marker
sat_i1 = plot(statei(1,1),statei(1,1),'r');
hold off
set(gca,'Xlim',[-max(statem(:,1))-1,max(abs(statem(:,1)))+1],'YLim',[-max(abs(statem(:,2)))-1,max(abs(statem(:,2)))+1],'Xtick',[],'YTick',[]);

%% 

%final animation loop
for i  = 2:length(t)
    if i==2 
        pause(5);
    end
    %BODY FRAME ANIMATION OF SATELLITE
 set(satb,'XData',stateb(i,1),'YData',stateb(i,2));
 set(satb1,'XData',stateb(1:i,1),'YData',stateb(1:i,2));
%INITIAL FRAME ANIMATION
%MOON
 set(moon,'XData',statem(i,1),'YData',statem(i,2));
 set(moon1,'XData',statem(1:i,1),'YData',statem(1:i,2));
%SATELLITE
 set(sat_i,'XData',statei(i,1),'YData',statei(i,2));
 set(sat_i1,'XData',statei(1:i,1),'YData',statei(1:i,2));
    pause(0.0000001)
end

