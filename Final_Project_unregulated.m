%unregulated
clear;
clc;
clf;
%from Karzbrun et. al.
W=20; %um 
D=33.5; %um^2/s 
R=50; %um 
FL=1.23e-2; %1/uM, fig.S4

alpha=146.7/60; %nM/s from Garamella J. et. al. table.S6
gamma=5.78e-04; %1/s from Bionumber ID 105190 

K=2; 
n=1.5;
P(1)=0.5; %nM
t=[0:0.0001:4]; %hr
l=[50;100;150;200;250;300]; %um

for j=1:6
    
L=l(j);
tauR=pi*R^2*L/D/W; 
for i=1:40000
    Tdiff=(L^2)/D/3600/9; %diffusion time from Karzbrun 
    if i<Tdiff*10000 %before main flow diffusion into compartment
        C=1e-04;
    else %Extract arrive compartment, conc. inside increasing
        C=1-exp(-i*0.0001*3600/tauR); %Extract conc.
    end
    dPdt(i)=alpha*C^2-P(i)*(gamma+1/tauR); 
    P(i+1)=P(i)+1/10000*3600*dPdt(i);


end
kk(j)=P(40000)*FL;

figure(1)
hold on
plot(t,P*FL)
end
xlabel('T(hours)')
ylabel('Fluorescence[AU]')
legend('50um','100um','150um','200um','250um','300um')
hold off
