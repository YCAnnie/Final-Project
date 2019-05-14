%sigma38 activator-repressor
clear;
clc;
clf;
%from Karzbrun et. al.
W=20; %um 
D=33.5; %um^2/s 
R=50; %um 
FL=1.23e-2; %1/uM, fig.S4

alpha=146.7/60; %nM/s from Garamella J. et. al. table.S6
beta=146.7/60*(0.4/7.65); %nM/s from Garamella J. et. al. table.S5 and S6
gamma=9.3/60/5000; %1/s from Garamella J. et. al. FigS13

K=5; 
n=2;  
A(1)=0.01; %nM
B(1)=0.01; %nM
C(1)=0.01; %nM
t=[0:0.0001:4]; %hr
l=[50;100;150;200;250;300]; %um

for j=1:6
    
L=l(j);
tauR=pi*R^2*L/D/W; 
for i=1:40000
    Tdiff=(L^2)/D/3600/9;%diffusion time from Karzbrun
    if i<Tdiff*10000 %before main flow diffusion into compartment
        C_onset=1e-04;
    else %Extract arrive compartment, conc. inside increasing
        C_onset=1-exp(-i*0.0001*3600/tauR); %Extract conc.
    end
    dAdt(i)=(alpha/(1+(B(i)/K)^n)*C_onset^2-A(i)*(gamma+1/tauR)); 
    A(i+1)=A(i)+1/10000*3600*dAdt(i);
    dBdt(i)=(beta*A(i)^n/(K+A(i)^n)*C_onset^2-B(i)*(gamma+1/tauR));
    B(i+1)=B(i)+1/10000*3600*dBdt(i);
    dCdt(i)=(alpha/(1+(B(i)/K)^n)*C_onset^2-C_onset(i)*(gamma+1/tauR));
    C_onset(i+1)=C_onset(i)+1/10000*3600*dCdt(i);

end
EE(j)=C_onset(40000)*FL;

figure(1)
hold on
plot(t,C_onset*FL)
end
xlabel('T(hours)')
ylabel('Fluorescence[AU]')
legend('50um','100um','150um','200um','250um','300um')
hold off
