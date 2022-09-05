function mainHH_Euler_IniguezDeOnzono_Martin
%==============================
% %% HODKING-HUXLEY MODEL NEURAL CELL MODEL
%==============================
clc
clear

close all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Constants set
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Cm=0.01; % Membrane Capacitance mF/cm^2
dt=0.05; % Time Step ms
t=0:dt:100; %Time Array ms
Iap=0.1; %External Current Applied
ENa=55.17; % mV Na+ reversal potential
EK=-72.14; % mV K+ reversal potential
El=-49.42; % mV Leakage reversal potential
gbarNa=1.2; % mS/cm^2 Na conductance
gbarK=0.36; % mS/cm^2 K conductance
gbarl=0.003; % mS/cm^2 Leakage conductance

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Definition of variables
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% We proceed to the initialization of the time-depending variables

m=zeros(1,length(t));
n=zeros(1,length(t));
h=zeros(1,length(t));
gNa=zeros(1,length(t));
gK=zeros(1,length(t));
gl=zeros(1,length(t));
INa=zeros(1,length(t));
IK=zeros(1,length(t));
Il=zeros(1,length(t));
V=zeros(1,length(t));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculation of initial values
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% We assume that the inital membrane voltage is -60 mV

V(1)=-60; % Initial Membrane voltage
v=V(1); % Refreshed membrane voltage

% Initial M-Calculation
am=0.1*(v+35)/(1-exp(-(v+35)/10));
bm=4.0*exp(-0.0556*(v+60));
m(1)=am/(am+bm); % Initial m-value

% Inital N-Calculation
an=0.01*(v+50)/(1-exp(-(v+50)/10));
bn=0.125*exp(-(v+50)/80);
n(1)=an/(an+bn); % Initial n-value

% Inital H-Calculation
ah=0.07*exp(-0.05*(v+60));
bh=1/(1+exp(-(0.1)*(v+30)));
h(1)=ah/(ah+bh); % Initial h-value

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MAIN LOOP CALCULATING VM FOR EACH INSTANT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i=1:length(t)-1
    
    %In this part you have to calculate the value of the conductances
    %(i.e. gna etc) from the value of the variables of the gates (i.e. m,n %and h)
    
    gNa(i) = gbarNa*(m(i)^3)*h(i);
    gK(i) = gbarK*(n(i)^4);
    
    % 2) Calculate the current by using the value of the conductance
    
    INa(i) = -gNa(i)*(V(i)-ENa);
    IK(i) = -gK(i)*(V(i)-EK);
    Il(i)= -gbarl*(V(i)-El);
    
    %INCLUDE YOUR CODE
    
    % 3) Use the Forward Euler Method to calculate the current step of VM
    %    i+1
    
    V(i+1) = V(i)+dt*(1/Cm)*(INa(i)+IK(i)+Il(i)+Iap);
    v=V(i);
    
    %INCLUDE YOUR CODE
    
    % 4) Use the forward Euler Method to calculate the values of the gates
    % (i.e. m(i+1), n(i+1) and h(i+1)
    
    am=0.1*(v+35)/(1-exp(-(v+35)/10));
    bm=4.0*exp(-0.0556*(v+60));
    an=0.01*(v+50)/(1-exp(-(v+50)/10));
    bn=0.125*exp(-(v+50)/80);
    ah=0.07*exp(-0.05*(v+60));
    bh=1/(1+exp(-(0.1)*(v+30)));
    
    m(i+1) = m(i) + dt*(am*(1-m(i))-bm*m(i));
    n(i+1)=n(i)+ dt*(an*(1-n(i))-bn*n(i));
    h(i+1)=h(i)+dt*(ah*(1-h(i))-bh*h(i));
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PLOT ACTION POTENTIAL
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure,
plot(t,V,'LineWidth',2),hold on
legend('Action Potential');
xlabel('Time (ms)','FontWeight','bold')
ylabel('Voltage (mV)','FontWeight','bold')
title('Voltage Change for Hodgkin-Huxley Model','FontWeight','bold')
set(gca,'FontSize',8)
set(gca,'FontWeight','bold')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PLOT ION CHANNEL CURRENTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure
plot(t(1:length(INa)),INa,'g',t(1:length(INa)),IK,'r',t(1:length(INa)),Il,'k','LineWidth',2),hold on
ylabel('Ion channel currents (uA)','FontWeight','bold')
xlabel('Time (ms)','FontWeight','bold')
legend('INa','IK','Ir');
set(gca,'FontSize',8)
set(gca,'FontWeight','bold')  

end