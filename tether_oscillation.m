function tether_oscillation()
global c n K rho L A Q 
 
n=3; %numebr of mode shapes
K=1e+2; % N/m
rho= 1; % kg/m
L= 1;   % m
A = [0 0 0]';
phi_0 = [pi/2 pi/2 pi/2]';
wave_speed = sqrt(K/rho);
for i=1:n
    lambda(:,i) = 2*L/i;
    omega(:,i)= wave_speed/lambda(i);
end
%initial_conditions = [sin(phi_0(1)) sin(phi_0(2)) sin(phi_0(3)) -omega(1)*cos(phi_0(1)) -omega(2)*cos(phi_0(2)) -omega(3)*cos(phi_0(3))]';
initial_conditions = [0 0 0 0 0 0]';
t0=0;
tf=2;
c= zeros(n,1);
Q= zeros(n,1);

options = odeset('RelTol',1e-12,'AbsTol',1e-12*ones(2*n,1));
[T,Z] = ode45(@do_physics,[t0 tf],initial_conditions,options);
%=============================== Plotting ================================%
figure(1)
plot(T,Z(:,1:3));

figure(2)
plot(T,Z(:,4:6));

x=(0:0.01:L)';
dt = (tf-t0)/100;
t= (t0:dt:tf)';
z = interp1(T,Z,t);

%figure(3)

for j=1:length(t)
    zz =zeros (length(x),1);
    for i=1:n
        tau(j,i) = z(j,i);
        phi(:,i) = sin((i*pi/L)*x);
        etha(:,i) = tau(j,i)*phi(:,i);
        element(:,i) = etha(:,i);
    end
    %plot(x, etha(:,1:3));
    for i=1:n
        zz = zz+ element(:,i);
    end
    figure(4)
    xlim([0 L]);
    plot(x, zz);
    ylim([-2 2]);
    pause(2*dt)
               
         
            
end
   
%===============================          ================================%
end
function dy = do_physics(t,y)
global c n K rho L Q
if t>=0.1 && t<=0.2
    F=20;
else
    F=0;
end
x_i = L/2;
for i=1:n
    c(i,1)= -(K/rho)*(i*pi/L)^2;  %*(i*pi/L)^2;
    Q(i,1)= F*sin(i*pi*x_i/L) ;
end
dy(4,1) = c(1)*y(1)+Q(1);
dy(5,1) = c(2)*y(2)+Q(2);
dy(6,1) = c(3)*y(3)+Q(3);

dy(1,1) = y(4);
dy(2,1) = y(5);
dy(3,1) = y(6);

end