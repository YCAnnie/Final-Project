%seperated compartments
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

K=10;
n=2;  
A(1)=1; %nM
B(1)=1; %nM
C(1)=1; %nM
t=[0:0.0001:5]; %hr
l=200; %um
d=[50;100;150;200;250;300]; %um

for j=1:6
L=d(j);
tauRA=pi*R^2*l/D/W; 
tauRB=pi*R^2*(l+L)/D/W; 
for i=1:50000
    TdiffA=(l^2)/D/3600/9; %diffusion time from Karzbrun
    if i<TdiffA*10000 %before main flow diffusion into compartment
        C_onsetA=1e-04;
    else %Extract arrive compartment, conc. inside increasing
        C_onsetA=1-exp(-i*0.0001*3600/tauRA); %Extract conc.
    end
    
    TdiffB=((l+L)^2)/D/3600/9; %diffusion into the repressor compartment
    if i<TdiffB*10000 %before main flow diffusion into compartment
        C_onsetB=1e-04;
    else %Extract arrive compartment, conc. inside increasing
        C_onsetB=1-exp(-i*0.0001*3600/tauRB); %Extract conc.
    end
    
    dAdt(i)=(alpha/(1+(B(i)/K)^n)*C_onsetA^2-A(i)*(gamma+0.5/tauRA)); 
    A(i+1)=A(i)+1/10000*3600*dAdt(i);
    dBdt(i)=(beta*A(i)^n/(K+A(i)^n)*C_onsetB^2-B(i)*(gamma+0.5/tauRB));
    B(i+1)=B(i)+1/10000*3600*dBdt(i);
    dCdt(i)=(alpha/(1+(B(i)/K)^n)*C_onsetA^2-C(i)*(gamma+0.5/tauRA));
    C(i+1)=C(i)+1/10000*3600*dCdt(i);

end
EE(j)=C(50000)*FL;

figure(1)
hold on
plot(t,C*FL)
end
xlabel('T(hours)')
ylabel('Fluorescence[AU]')
legend('50um','100um','150um','200um','250um','300um')
hold off
