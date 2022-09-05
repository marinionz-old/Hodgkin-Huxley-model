function mainHH_Euler_120bpm_IniguezDeOnzono_Martin
% ==============================
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


%% IaP calculation for the seeked Frequency (120 Hz)

seeked_Freq = 120;

% Number of iterations 
iter=60/dt; % 60ms / 0.05 ms time step
% 60 ms okay for 2 peaks at low frequencies of 50 Hz, which need 20 ms per
% peak (at least two peaks for finding the differece in time and thus, 
% calculating the  frequency from that)

% We are going to check which is the period per Iap from 0.1 uA (initial one,
% yielding a Frequency of around 80 Hz) to 1 uA, as we consider that value to
% be very high and I have seen in a preliminary study inserting random
% values on Iap that the frequency increases as Iap increases, and at 1uA it
% is very high.

possible_Iap = Iap:0.01:1; 

for j=1:length(possible_Iap)
    
    % We initialize the variables per iteration
    
    m=zeros(1,iter);
    n=zeros(1,iter);
    h=zeros(1,iter);
    gNa=zeros(1,iter);
    gK=zeros(1,iter);
    gl=zeros(1,iter);
    INa=zeros(1,iter);
    IK=zeros(1,iter);
    Il=zeros(1,iter);
    V=zeros(1,iter);
    Iap=possible_Iap(j);
    
    % Inital values
    V(1)=-60; % Initial Membrane voltage
    v=V(1);
    am=0.1*(v+35)/(1-exp(-(v+35)/10));
    bm=4.0*exp(-0.0556*(v+60));
    m(1)=am/(am+bm); % Initial m-value
    an=0.01*(v+50)/(1-exp(-(v+50)/10));
    bn=0.125*exp(-(v+50)/80);
    n(1)=an/(an+bn); % Initial n-value
    ah=0.07*exp(-0.05*(v+60));
    bh=1/(1+exp(-(0.1)*(v+30)));
    h(1)=ah/(ah+bh); % Initial h-value

    for i=1:iter-1
        
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
        v=V(i+1);
        
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
    
    % The peaks of the signal are found to calculate its period (time
    % between max. potential peaks)
    [pks,locations]=findpeaks(V);
    
    % The difference in index in peaks is the period when multiplied by the
    % time step.
    for k=1:length(locations)-1
        locs_dif(k) = locations(k+1)-locations(k);
    end
    
    % The mean of that differences of peaks is used
    avg_locs_dif=mean(locs_dif);
    
    % That mean is multiplied by the time step to transform it into the ms
    % domain and we divide 1000 by that number to go from KHz to Hz (we
    % were using 3 ms)
    freq_Iap(j)=1000/(avg_locs_dif*dt); % Frequency in Hz
   
end

% The frequency with the least difference with the seeked frequency (in
% this case 120 Hz) is used as Iap for the HH model (from now on, the same
% code as mainHH_Euler is used)

[val,idx]=min(abs(freq_Iap-seeked_Freq));
Iap = possible_Iap(idx);

% Plot of how the frequency depends on the Iap introduced
figure
scatter(possible_Iap,freq_Iap),hold on

xlabel('Iap (uA)','FontWeight','bold') 
ylabel('Action Potential Frequency (Hz)','FontWeight','bold') 
title('Action Potential Frequency (Hz) depending on the stimuli introduced','FontWeight','bold') 
set(gca,'FontSize',9)
set(gca,'FontWeight','bold')  

% Calculations of the regressed line and the regressed value for the Seeked
% Frequency (120 Hz in this case).
fitted_eq=polyfit(possible_Iap,freq_Iap,1);
fitted_line=polyval(fitted_eq,0:1);
plot(0:1,fitted_line,'--');
regressed_Iap= (seeked_Freq-fitted_eq(2))/fitted_eq(1);
fprintf('The regressed Iap value is %.3f uA\nwhile the actual Iap that showed the closest value to 120 Hz was %.3f uA.\n',regressed_Iap,Iap);

%% Performance of the model with the calculated Iap

% We erase the variables, so there are no problems with overwriting their
% values

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
    
%%%% Initial Conditions %%%%

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
    v=V(i+1);
    
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
title('Voltage Change for Hodgkin-Huxley Model @ 120 Hz','FontWeight','bold') 
set(gca,'FontSize',9)
set(gca,'FontWeight','bold')  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PLOT ION CHANNEL CURRENTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure
plot(t(1:length(INa)),INa,'g',t(1:length(INa)),IK,'r',t(1:length(INa)),Il,'k','LineWidth',2),hold on
ylabel('Ion channel currents','FontWeight','bold') 
xlabel('Time (ms)','FontWeight','bold') 
legend('INa','IK','Ir');
set(gca,'FontSize',9)
set(gca,'FontWeight','bold')  

end