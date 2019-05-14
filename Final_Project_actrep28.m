%sigma28 activator-repressor
clear;
clc;
%from Karzbrun et. al.
W=20; %um 
D=33.5; %um^2/s 
R=50; %um 
FL=1.23e-2; %1/uM, fig.S4

alpha_1=146.7/60; %nM/s from Garamella J. et. al. table.S6
beta=146.7/60*(5.69/7.65); %nM/s from Garamella J. et. al. table.S5 and S6
gamma=5.78e-04; %Bionumber ID 105190

K=2; 
n=3;  
A(1)=0.5; %nM
B(1)=0.2; %nM
C(1)=0.9; %nM
t=[0:0.0001:10]; %hr
l=[50;100;150;200;250;300]; %um

for j=1:6
    
L=l(j);
tauR=pi*R^2*L/D/W; 
for i=1:100000
    Tdiff=(L^2)/D/3600/9; %diffusion time from Karzbrun
    if i<Tdiff*10000 %before main flow diffusion into compartment
        C_onset=1e-04;
    else %Extract arrive compartment, conc. inside increasing
        C_onset=1-exp(-i*0.0001*3600/tauR); %Extract conc.
    end
    dAdt(i)=(alpha_1/(K+B(i)^n)*C_onset^2-A(i)*(gamma+1/tauR)); 
    A(i+1)=A(i)+1/10000*3600*dAdt(i);
    dBdt(i)=(beta*A(i)^n/(K+A(i)^n)*C_onset^2-B(i)*(gamma+1/tauR));
    B(i+1)=B(i)+1/10000*3600*dBdt(i);
    dCdt(i)=(beta*A(i)^n/(K+A(i)^n)*C_onset^2-C_onset(i)*(gamma+1/tauR));
    C_onset(i+1)=C_onset(i)+1/10000*3600*dCdt(i);

end
figure(1)
DD(j)=C_onset(100000)*FL;
hold on
plot(t,C_onset*FL)
end
xlabel('T(hours)')
ylabel('Fluorescence[AU]')
legend('50um','100um','150um','200um','250um','300um')
hold off
