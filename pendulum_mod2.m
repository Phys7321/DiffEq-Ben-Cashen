function [period,sol] = pendulum_mod2(omega,theta0,thetad0,grph) 
% Finds the period of a nonlinear pendulum given the length of the pendulum
% arm and initial conditions. All angles in radians.

%Setting initial conditions
if nargin == 0
    error('Must input length and initial conditions')
end
if nargin == 1
   theta0 = pi/2;
   thetad0 = 0;
   grph = 0;
end
if nargin == 2
    thetad0 = 0;
    grph = 1;
end
if nargin == 3
    grph = 1;
end

m = 1;
g = 9.81;
R = g/omega^2;
T = 2*pi/omega;

% number of oscillations to graph
N = 10;


tspan = [0 N*T];
%opts = odeset('events',@events,'refine',6); %Here for future event finder
opts = odeset('refine',6);
r0 = [theta0 thetad0];
[t,w] = ode45(@proj,tspan,r0,opts,g,R);
sol = [t,w];

ind = w(:,2).*circshift(w(:,2), [-1 0]) <= 0;     
%ind = chop(ind,4);
period = 2*mean(diff(t(ind)));                    

if nargin == 3
    % Calculating total energy
    E_k = (1/2)*m*(R*w(:,2)).^2;     % Kinetic energy
    E_p = m*g*R*(1-cos(w(:,1)));     % Potential energy
    E = E_k + E_p;
    delta = (E(:) - E(1))/E(1);

    cyc_start = find(ind);
    cyc_start = cyc_start(1:2:end);     % Finding index for start of each cycle
    
    % Average energy in every cycle
    Ek_avg = [];
    Ep_avg = [];
    for i = 1:N
        Ek_avg(i) = mean(E_k(cyc_start(i):cyc_start(i+1)));
        Ep_avg(i) = mean(E_p(cyc_start(i):cyc_start(i+1)));
    end
end

Etot_avg = Ek_avg + Ep_avg;

figure(1)
plot(t,delta,'o')
title('Energy change during one period')
xlim([0,T])
xlabel('time')
ylabel('Relative energy change (\Delta)')
    
figure(2)
subplot(2,1,1)
plot(t,w(:,1),'r')
title('Position vs time')
xlabel('time')
ylabel('\theta')
subplot(2,1,2)
plot(t,w(:,2),'g')
title('Velocity vs time')
xlabel('time')
ylabel('d\theta/dt')

figure(3)
plot(w(:,1),w(:,2),'-')
title('Phase diagram of simple armonic motion')
xlabel('\theta')
ylabel('d\theta/dt')

figure(4)
plot(1:N,Ek_avg,'o',1:N,Ep_avg,'o',1:N,Etot_avg,'o')
title('Average kinetic and potential energy during full cycle')
legend('Avg. Kinetic Energy','Avg. Potential Energy','Avg. Total Energy','Location','best')
xlabel('Cycles')
ylabel('Energy')

plot(t,E,'-')
xlim([0,T])
title('Total energy during one period')
xlabel('time')
ylabel('Total Energy')
    
end
%-------------------------------------------
%
function rdot = proj(t,r,g,R)
    rdot = [r(2); -g/R*sin(r(1))];     % sin(r1) not expanded to account for large amplitude oscillations
end