%% %% Kundur SRM
%%%% Please set the parameters based on the guidance provided in the Readme.md file 
%both torque and AVR
aaa=0;
tsimave=0;
while aaa<1
    aaa=aaa+1;
clearvars -except aaa tsimave
clc

tic
% This parameter chooses the predictor method. It can take values 4 and 5. 4 is for the first-order Adams-Bashford, and 5 is for the second-order Adams-Bashford. 
cp=5; %Predictor Choice
% This parameter chooses the corrector method. It can take values 2 and 3. 2 is for the second-order BDF, and 3 is for the second-order Adams-Moulton. 
cc=3; %Corrector Choice
%This parameters chooses the quantization method
%1 uniform quantization, 2 Nonuniform quantization 
cq=1; %Quantization method

t=0; %Starting point of the integration
h=0.001; %The first time-step size
holdd=h;
hmin=0.001; %The minimum acceptable time-step size
hmax=1; %The maximum acceptable time-step size

tsim=21; %Time of the simulation
%Initial state variables values
[ssi]=ini;
xx=ssi;
dae=[16 76];
cn=8;
T=[0.021 0.022 0.023 0.024 0.0041 0.0042 0.0043 0.0044]; %Time of sampling by the controller
%Input
Vf1=0.000869750976562500;
Vf2=0.000534057617187500;
Vf3=0.000549316406250000;
Vf4=0.000381469726562500;
Tm1=0.954788208007813;
Tm2=0.546951293945313;
Tm3=0.433273315429688;
Tm4=0.224349975585938;
avrov1=[1.03;0;0;Vf1;Vf1;Vf1;0;0];
avrov2=[1.01;0;0;Vf2;Vf2;Vf2;0;0];
avrov3=[1.03;0;0;Vf3;Vf3;Vf3;0;0];
avrov4=[1.01;0;0;Vf4;Vf4;Vf4;0;0];
avrov=[avrov1 avrov2 avrov3 avrov4];

s=size(xx);
s=s(1);
xpre1=xx;
xpre2=xx;
times(1)=t;
xhist(1:s,1)=xx; %Allocating a space to save the results

iterationsc=0; %Allocating a space to save the number of Newoton iterations
tsc=0; %Number of time steps
mtsc=0; %Number of accepted time steps with the manimim value
%Allocating the space to save the old controllers' signal
etpre2j=zeros(1,cn);
etprej=zeros(1,cn);
etj=zeros(1,cn);
tg=T; %Time of the events
k=1;
ETOL=3e-4; %The tolerance to compare with the error estimate
fitertt=0;
dhiterationsc=0;

while t<tsim
        %Setting old times
    if k>2
        tpre3=tpre2;
    else
        tpre3=0;
    end
    if k>1
        tpre2=t-holdd;
    else
        tpre2=0;    
    end
    flagg=0;
    tpre1=t;
    t=t+h;
    
    %manual time control
    if t>1 && tpre1<1
        t=1;
        h=t-tpre1;
    end
    if t>12 && tpre1<12
        t=12;
        h=t-tpre1;
    end
%     if t>3.1 && tpre1<3.1
%         t=3.1;
%         h=t-tpre1;
%     end
%     if t>13 && tpre1<13
%         t=13;
%         h=t-tpre1;
%     end
    %Setting old solutions
    if k>2
        xpre3=xhist(1:s,k-2);
    else
        xpre3=0;
    end
    if k>1
    xpre2=xhist(1:s,k-1);
    else
        xpre2=xpre1;
    end
    xpre1=xhist(1:s,k);
    k=k+1;
    %Setting new and old controllers's signal
    if t<min(T)
        et=[Tm1 Tm2 Tm3 Tm4,Vf1 Vf2 Vf3 Vf4];
        etpre=[Tm1 Tm2 Tm3 Tm4,Vf1 Vf2 Vf3 Vf4];
        etpre2=[Tm1 Tm2 Tm3 Tm4,Vf1 Vf2 Vf3 Vf4];
        avrov1=[1.03;0;0;Vf1;Vf1;Vf1;0;0];
        avrov2=[1.01;0;0;Vf2;Vf2;Vf2;0;0];
        avrov3=[1.03;0;0;Vf3;Vf3;Vf3;0;0];
        avrov4=[1.01;0;0;Vf4;Vf4;Vf4;0;0];
    else
        tgpre=tg;
        [t,conti,tg,mm]=detert(tpre1,t,tg,T,cn,tsim); %Finding the event if there is any
        if mm==0
            etpre2=etpre;
            etpre=et;
            et=et;
            tg=tg;
        else
            etpre2=etpre;
            etpre=et;
            [et,avrov]=contex(conti,et,T,cn,avrov,xx); %Calculating the controllers' signal in case of any event
            et=quz(et,cq);
            h=t-tpre1; %Zero-crossing
            flagg=1;
        end
    end
    %Summation of the controllers' signal
    times(k)=t;

    [alpha0,alpha1,alpha2,betha1,betha2]=coefs(cp,h,holdd); %calculating the coefficients of the predictor
    xx0=predictor(xx,xpre1,et,etpre,t,tpre1,h,betha1,betha2,dae); %predictor
    [alpha0,alpha1,alpha2,betha1,betha2]=coefs(cc,h,holdd); %calculating the coefficients of the corrector
    [iterations,xx,fitert,dhiterations]=NewtonNL(h,xx0,xpre1,xpre2,alpha1,alpha2,alpha0,betha1,betha2,tpre1,et,etpre,t,dae); %The main function (Newton solver)
    [dn,r]=lote(cp,cc,xx0,xx,h,holdd); %Error estimate
    iterationsc=iterations+iterationsc; %Counting the number of Newton iterations
    fitertt=fitertt+fitert;
    dhiterationsc=dhiterationsc+dhiterations;


   %Calculating the new time step and saving the results for the current time step         
   hpre=h;
    if dn>ETOL
        hc=[hmin,0.5*h;];
        h=max(hc);
        if h==hmin
            mtsc=mtsc+1;
            tsc=tsc+1;
            disp('The minimum time step is accepted')
            xhist(1:s,k)=xx;
            x2=xx(2,1);
            if k==2
                hhist(1,1)=h;
                hhist(2,1)=dn;
                hhist(3,1)=r;
                ethist(:,1)=et;
            end
            hhist(1,k)=h;
            hhist(2,k)=dn;
            hhist(3,k)=r;
            ethist(:,k)=et;
            holdd=hpre;
        else
                t=t-hpre;
                tpre1=tpre2;
                if k>2
                    tpre2=tpre3;
                else
                    tpre2=0;
                end
                xpre1=xpre2;
                if k>2
                    xpre2=xpre3;
                else
                    xpre2=0;
                end
                if flagg==1
%                     tg=tg-T;
                    tg=tgpre;
                    et=etpre;
                    etpre=etpre2;
                end
                k=k-1;

        end
    else
        hc=[1.25*h,hmax];
        h=min(hc);
        tsc=tsc+1;
        xhist(1:s,k)=xx;
        x2=xx(2,1);
        if k==2
            hhist(1,1)=h;
            hhist(2,1)=dn;
            hhist(3,1)=r;
            ethist(:,1)=et;
        end
        hhist(1,k)=h;
        hhist(2,k)=dn;
        hhist(3,k)=r;
        ethist(:,k)=et;
        holdd=hpre;
    end


end
tsimave=tsimave+toc;
end
endtsimave=tsimave/aaa;

% %%%Output
% %Generator1
% ssi1=[psif1 psiq11 omega1 deltas1 ]'; %4
% %Generator2
% ssi2=[psif2 psiq12 omega2 deltas2 ]'; %4
% %Generator3
% ssi3=[psif3 psiq13 omega3 deltas3 ]'; %4
% %Generator4
% ssi4=[psif4 psiq14 omega4 deltas4 ]'; %4
% %other Generator
% ssi5=[ifd1 iq11 vd1 vq1 id1 iq1 psid1 psiq1 ifd2 iq12 vd2 vq2 id2 iq2 psid2 psiq2 ifd3 iq13 vd3 vq3 id3 iq3 psid3 psiq3 ifd4 iq14 vd4 vq4 id4 iq4 psid4 psiq4]'; %32
% %Voltage
% ssi6=[vx1 vy1 vx2 vy2 vx3 vy3 vx4 vy4 vx5 vy5 vx6 vy6 vx7 vy7 vx8 vy8 vx9 vy9 vx10 vy10 vx11 vy11]'; %22
% %Current
% ssi7=[ix1 iy1 ix2 iy2 ix3 iy3 ix4 iy4 ix5 iy5 ix6 iy6 ix7 iy7 ix8 iy8 ix9 iy9]'; %18
% ssi=[ssi1;ssi2;ssi3;ssi4;ssi5;ssi6;ssi7];

%Plotting the results


% figure(21) %speed deviation 1
% plot(times,xhist(3,:))
% hold on
% xlabel('Time','FontSize', 24)
% ylabel('Speed deviation 1','FontSize', 24)
% set(gca,'FontSize',18)
% figure(22) %speed deviation 2
% plot(times,xhist(7,:))
% hold on
% xlabel('Time','FontSize', 24)
% ylabel('Speed deviation 2','FontSize', 24)
% set(gca,'FontSize',18)
% figure(23) %speed deviation 3
% plot(times,xhist(11,:))
% hold on
% xlabel('Time','FontSize', 24)
% ylabel('Speed deviation 3','FontSize', 24)
% set(gca,'FontSize',18)
% figure(24) %speed deviation 4
% plot(times,xhist(15,:))
% hold on
% xlabel('Time','FontSize', 24)
% ylabel('Speed deviation 4','FontSize', 24)
% set(gca,'FontSize',18)

figure (1) %voltage
plot(times,sqrt((xhist(49,:)).^2+(xhist(50,:)).^2))
hold on
% plot(times,sqrt((xhist(55,:)).^2+(xhist(56,:)).^2))
% plot(times,sqrt((xhist(65,:)).^2+(xhist(66,:)).^2))
xlabel('Time','FontSize', 24)
ylabel('Voltage','FontSize', 24)
set(gca,'FontSize',18)
% axes('position',[.3 .2 .25 .25])
% box on % put box around new pair of axes
% signal=sqrt((xhist(49,:)).^2+(xhist(50,:)).^2);
% indexOfInterest = (times > 1.5) & (times < 3.5); % range of t near perturbation
% plot(times(indexOfInterest),signal(indexOfInterest)) % plot on new axes
% hold on


figure(33)
% plot(times,xhist(3,:))
% hold on
% plot(times,xhist(7,:))
plot(times,xhist(11,:))
hold on
% plot(times,xhist(15,:))
xlabel('Time','FontSize', 24)
ylabel('Speed deviation','FontSize', 24)
set(gca,'FontSize',18)

figure (31) %Governor
% plot(times,ethist(1,:))
% hold on
% plot(times,ethist(2,:))
stairs(times,ethist(3,:))
% plot(times,ethist(4,:))
hold on
xlabel('Time','FontSize', 24)
ylabel('Governor signal','FontSize', 24)
set(gca,'FontSize',18)
figure (32) %Exciter
% plot(times,ethist(5,:).*1985.93939393939,2621)
hold on
% plot(times,ethist(6,:).*2621.44000000000)
stairs(times,ethist(7,:).*2032.12403100775)
% plot(times,ethist(8,:).*2449.94392523364)
hold on
xlabel('Time','FontSize', 24)
ylabel('Exciter signal','FontSize', 24)
set(gca,'FontSize',18)
% [1985.93939393939,2621.44000000000,2032.12403100775,2449.94392523364]

figure(34)
plot(times,hhist(1,:))
hold on
xlabel('Time','FontSize', 24)
ylabel('Step size','FontSize', 24)
set(gca,'FontSize',18)

iterationsc=iterationsc
fitertt=fitertt
dhiterationsc=dhiterationsc
endtsimave=endtsimave
SRMt=times;
SRMx=xhist;
SRMh=ethist;

function [t,minl,tg,mm]=detert(tpre1,t,tg,T,cn,tsim)
temp=zeros(1,cn);
temp=temp+tsim*2;
cni=0;
mm=0;
while cni<cn
    cni=cni+1;
    if tpre1<tg(1,cni) && tg(1,cni)<=t+0.0000001
        temp(1,cni)=tg(1,cni);
        mm(1,cni)=1;
    end
end
if mm==0
    t=t;
    tg=tg;
    minl=0;
else
    temp=round(temp,6); 
    minl=find(temp==min(temp));
    sminl=size(minl);
    sminl=sminl(1,2);
    j=0;
    while j<sminl
        j=j+1;
        t=min(temp);
        tg(1,minl(1,j))=tg(1,minl(1,j))+T(1,minl(1,j));
    end
    tg=round(tg,6);
end
end

%This function updates the controllers' signal in case of the event 
function [etj,avrov]=contex(conti,etj,T,cn,avrov,xx)
jj=size(conti);
jj=jj(1,2);
j=0;
while j<jj
    j=j+1;
    contj=conti(1,j);
    if contj==1 || contj==2 || contj==3 || contj==4
        etj(1,contj)=torque(etj(1,contj),T(1,contj),contj,xx);
        etj(1,contj)=max(etj(1,contj),0); %upper limnit anti wind up
        etj(1,contj)=min(etj(1,contj),2); %lower limnit anti wind up
%         etj(1,contj)=min(etj(1,contj),1); %lower limnit anti wind up
    elseif contj==5 || contj==6 || contj==7 || contj==8
        [etj(1,contj),avrov(:,contj-4)]=avr(avrov(:,contj-4),T(1,contj),xx,contj); %Controller's equation
        if contj==7
        etj(1,contj)=max(etj(1,contj),0); %upper limnit anti wind up
        etj(1,contj)=min(etj(1,contj),2/2032.12403100775); %lower limnit anti wind up
%         etj(1,contj)=min(etj(1,contj),0.01); %lower limnit anti wind up
        end
    end
end
end

function et2=torque(et2,T,contj,xx)
%Data
Tsm=0.9; %0.3
sigma=0.04;
if contj==1
    zo=0.954788208007813;
    wS=1+xx(3,1);
elseif contj==2
    zo=0.546951293945313;
    wS=1+xx(7,1);
elseif contj==3
    zo=0.433273315429688;
    wS=1+xx(11,1);
else
    zo=0.224349975585938;
    wS=1+xx(15,1);
end
%initialization (old values)
T2=T;
x9=et2;
%equations
% x9=(x9*(Tsm*sigma*wN-sigma*wN*T2)+T2*(zo*sigma*wN-wS+wN))/(Tsm*sigma*wN);
x9=(x9*(Tsm*sigma-sigma*T2)+T2*(zo*sigma-wS+1))/(Tsm*sigma);
%assigning old values for the next sample
et2=x9;
end

function [et1,avrov]=avr(avrov,T1,xx,contj)
%Data
Tmm=0.2; %0.2
Te=0.12; %0.12
if contj==5
    vo=1.03;
    Vter=sqrt(xx(49,1)^2+xx(50,1)^2);
elseif contj==6
    vo=1.01;
    Vter=sqrt(xx(51,1)^2+xx(52,1)^2);
elseif contj==7
    vo=1.03;
    Vter=sqrt(xx(53,1)^2+xx(54,1)^2);
else
    vo=1.01;
    Vter=sqrt(xx(55,1)^2+xx(56,1)^2);
end
kp=0.003; %0.01
ki=0.001; %0.01
Ge=1;
kf=0.8; %0.1
Tf=0.9; %0.9

%initialization (old values)
x1=avrov(1,1);
x2=avrov(2,1);
x3=avrov(3,1);
x4=avrov(4,1);
x5=avrov(5,1);
x61=avrov(6,1);
x7=avrov(7,1);
x8=avrov(8,1);
%equations
x1=((Tmm-T1)*x1+T1*Vter)/Tmm;
x2=vo-x1-x8;
x3=kp*x2;
x4=x4+T1*ki*x2;
x5=x4+x3;
x62=(x61*(Te-T1)+T1*Ge*x5)/Te;
x7=(kf/T1)*(x62-x61);
% x7=(kf/T1)*(x62);
x8=(Tf*x8-T1*x8+T1*x7)/Tf;
et1=x62;
%assigning old values for the next sample
avrov=[x1 x2 x3 x4 x5 x62 x7 x8]';
end

function [iterations,xx,fitert,dhiterations]=NewtonNL(h,xx,xpre1,xpre2,alpha1,alpha2,alpha0,betha1,betha2,tpre1,et,etpre,t,dae)
NTOL=1e-4; %Tolerance for the Newton solver
iterations=0;
dhiterations=0;
a=1;
fitert=0;
    while (iterations<30) && (a>NTOL)
    xxpre=xx;
    iterations=iterations+1;
    %Dishonest Newton
%     if iterations==1 || iterations>5
%         dhiterations=dhiterations+1;
%         [J,fiter]=Jacob(h,xx,xpre1,xpre2,alpha1,alpha2,alpha0,betha1,betha2,tpre1,et,etpre,t,dae); %calculating the Jacobian
%         fitert=fitert+fiter;
%      end
    
    [J,fiter]=Jacob(h,xx,xpre1,xpre2,alpha1,alpha2,alpha0,betha1,betha2,tpre1,et,etpre,t,dae); %calculating the Jacobian
    fitert=fitert+fiter;
    bb=func(h,xx,xpre1,xpre2,alpha1,alpha2,alpha0,betha1,betha2,tpre1,et,etpre,t,dae);
    fitert=fitert+1;
    cc=J\bb;
    yy=xx-cc; %Newton
    xx=yy;
    a=abs(xx-xxpre);
%     a=abs((xx-xxpre)./xx);
    a=max(a);
    if iterations==30
        disp('Did not converged')
    end
    end
end

% This function computes the Jacobian
function [J,fiter]=Jacob(h,xx,xpre1,xpre2,alpha1,alpha2,alpha0,betha1,betha2,tpre1,et,etpre,t,dae)
s=size(xx);
s=s(1);
delta=1e-4;
i=0;
J=zeros(s,s);
fiter=0;
while i<s
    i=i+1;
    f0=func(h,xx,xpre1,xpre2,alpha1,alpha2,alpha0,betha1,betha2,tpre1,et,etpre,t,dae);
    fiter=fiter+1;
    yy=xx;
    yy(i)=xx(i)+delta;
    f1=func(h,yy,xpre1,xpre2,alpha1,alpha2,alpha0,betha1,betha2,tpre1,et,etpre,t,dae);
    fiter=fiter+1;
    J(:,i)=(f1-f0)/delta;
end
end

% This function determines the predictor formulation
function [xxp]=predictor(xx,xpre1,et,etpre,t,tpre1,h,betha1,betha2,dae)
    xxp1=xx+h*(((betha1)*evaleval(xx,t,et))+((betha2)*evaleval(xpre1,tpre1,etpre)));
    xxp2=xx+evaleval(xx,t,et);

if dae(1,1)==0
    xxp=xxp2;
elseif dae(1,2)==0
    xxp=xxp1;
else
    xxp=[xxp1(1:dae(1,1),1);xxp2(dae(1,1)+1:end,1)];
end
end


%This function determines the corrector formulation
function [f]=func(h,xx,xpre1,xpre2,alpha1,alpha2,alpha0,betha1,betha2,tpre1,et,etpre,t,dae) 
    f1=h*(betha1*evaleval(xx,t,et)+betha2*evaleval(xpre1,tpre1,etpre))+(alpha0*xx)+(alpha1*xpre1)+(alpha2*xpre2);
    f2=evaleval(xx,t,et);

if dae(1,1)==0
    f=f2;
elseif dae(1,2)==0
    f=f1;
else
    f=[f1(1:dae(1,1),1);f2(dae(1,1)+1:end,1)];
end
end


%This function computes the coefficients based on the step sizes
function [alpha0,alpha1,alpha2,betha1,betha2]=coefs(c,h,holdd)
if c==1
    % corrector coefficients
    % BDF-1
    alpha1=1;
    alpha2=0;
    alpha0=-1;
    betha1=1;
    betha2=0;

elseif c==2
    %corrector coefficients
    % BDF-2
    alpha1=1+h/holdd;
    alpha2=-(h^2)/((h+holdd)*holdd);
    alpha0=-(2*h+holdd)/(h+holdd);
    betha1=1;
    betha2=0;

elseif c==3
    %corrector coefficients
    % Adams-Moulton-2
    alpha1=1;
    alpha2=0;
    alpha0=-1;
    betha1=1/2;
    betha2=1/2; 

elseif c==4
    %predictor coefficients
    %AB-1
    alpha1=0;
    alpha2=0;
    alpha0=0;
    betha1=1;
    betha2=0;
    
elseif c==5
    %predictor coefficients
    %AB-1
    alpha1=0;
    alpha2=0;
    alpha0=0;
    betha1=1+(h/(holdd*2));
    betha2=-h/(holdd*2);
end
end



%This function computes the error estimate for predictor corrector methods (Milne's estimate)
function [dn,r]=lote(cp,cc,yp,yy,h,holdd)
global ETOL
ETOL=3e-6;
if cc==1 && cp==4
    dn=(yy-yp)*(-1/2)/((-1/2)-(1/2));
    dn=max(abs(dn(:)));
    r=0.9*ETOL/(dn*h);
    r=nthroot(abs(r),2);
% elseif cc==1 && cp==5
%     dn=(yy-yp)/3*(1+holdd/h);
%     r=0.9*ETOL/(dn*h);
%     r=nthroot(r,3);  
% elseif cc==2 && cp==4
%     dn=(yy-yp)/3*(1+holdd/h);
%     r=0.9*ETOL/(dn*h);
%     r=nthroot(r,3);
elseif cc==2 && cp==5
    dn=(yy-yp)*(-2/9)/((-2/9)-(5/12));
    dn=max(abs(dn(1:4,1)));
    r=0.9*ETOL/(dn*h);
    r=nthroot(r,3);    
% elseif cc==3 && cp==4
%     dn=(yy-yp)/3*(1+holdd/h);
%     r=0.9*ETOL/(dn*h);
%     r=nthroot(r,3);                            
elseif cc==3 && cp==5
    dn=(yy-yp)/(3*(1+holdd/h));
    dn=max(abs(dn(:,1)));
    r=0.9*ETOL/dn;
    r=nthroot(r,3);
end
end


%This function does the quantification
function [x]=quz(x,cq)
%Uniform Quantization
if cq==1 
    range=1;
    q=18;
    x=round(x*(2^q))*range/(2^q);
%NonUniform Quantization
elseif cq==2 
    xp=x;
    if xp>0.8
        x=2*x;
    else
        x=(2/3)*x+(4/3);
    end    
    range=1;
    q=10;
    x=round(x*(2^q))*range/(2^q);
    if xp>0.8
        x=x/2;
    else
        x=((3/2)*x)-2;
    end
end
end

%This function defines the system's formulation
function [f]=evaleval(xx,t,et)
j=sqrt(-1);
Pshsc=0;
Qshsc=0;
%Parameters
pp=para;
Xd=pp(1,1); Xq=pp(1,2); Xl=pp(1,3); Xpd=pp(1,4); Xpq=pp(1,5);
Tpdo_s=pp(1,6); Tpqo_s=pp(1,7); Ra=pp(1,8); fN=pp(1,9); HH12=pp(1,10);
HH34=pp(1,11); Xt=pp(1,12); Pload7=pp(1,13); Qload7=pp(1,14); Pload9=pp(1,15);
Qload9=pp(1,16); Re56=pp(1,17); Xe56=pp(1,18); Re67=pp(1,19); Xe67=pp(1,20);
Re1011=pp(1,21); Xe1011=pp(1,22); Re910=pp(1,23); Xe910=pp(1,24); Re78=pp(1,25);
Xe78=pp(1,26); Re89=pp(1,27); Xe89=pp(1,28); wN=pp(1,29); TB=pp(1,30);
Ldd=pp(1,31); Ldf=pp(1,32); Lff=pp(1,33); Tpdo=pp(1,34); Rf=pp(1,35);
Lqq=pp(1,36); Lqq1=pp(1,37); Lq1q1=pp(1,38); Tpqo=pp(1,39); Rq1=pp(1,40);
%Assiginng variables
psif1=xx(1,1); psiq11=xx(2,1); omega1=xx(3,1); deltas1=xx(4,1);
psif2=xx(5,1); psiq12=xx(6,1); omega2=xx(7,1); deltas2=xx(8,1); 
psif3=xx(9,1); psiq13=xx(10,1); omega3=xx(11,1); deltas3=xx(12,1);
psif4=xx(13,1); psiq14=xx(14,1); omega4=xx(15,1); deltas4=xx(16,1);
ifd1=xx(17,1); iq11=xx(18,1); vd1=xx(19,1); vq1=xx(20,1); id1=xx(21,1);
iq1=xx(22,1); psid1=xx(23,1); psiq1=xx(24,1);
ifd2=xx(25,1); iq12=xx(26,1); vd2=xx(27,1); vq2=xx(28,1); id2=xx(29,1); iq2=xx(30,1);
psid2=xx(31,1); psiq2=xx(32,1);
ifd3=xx(33,1); iq13=xx(34,1); vd3=xx(35,1); vq3=xx(36,1); id3=xx(37,1); iq3=xx(38,1);
psid3=xx(39,1); psiq3=xx(40,1);
ifd4=xx(41,1); iq14=xx(42,1); vd4=xx(43,1); vq4=xx(44,1); id4=xx(45,1); iq4=xx(46,1);
psid4=xx(47,1); psiq4=xx(48,1);
vx1=xx(49,1); vy1=xx(50,1); vx2=xx(51,1); vy2=xx(52,1); vx3=xx(53,1); vy3=xx(54,1);
vx4=xx(55,1); vy4=xx(56,1); vx5=xx(57,1); vy5=xx(58,1); vx6=xx(59,1); vy6=xx(60,1);
vx7=xx(61,1); vy7=xx(62,1); vx8=xx(63,1); vy8=xx(64,1); vx9=xx(65,1); vy9=xx(66,1);
vx10=xx(67,1); vy10=xx(68,1); vx11=xx(69,1); vy11=xx(70,1);
ix1=xx(71,1); iy1=xx(72,1); ix2=xx(73,1); iy2=xx(74,1); ix3=xx(75,1); iy3=xx(76,1);
ix4=xx(77,1); iy4=xx(78,1); ix5=xx(79,1); iy5=xx(80,1); ix6=xx(81,1); iy6=xx(82,1);
ix7=xx(83,1); iy7=xx(84,1); ix8=xx(85,1); iy8=xx(86,1); ix9=xx(87,1); iy9=xx(88,1);
ixg3=xx(89,1); iyg3=xx(90,1); ixshsc=xx(91,1); iyshsc=xx(92,1);

%Input
Vf1=et(1,5);
Vf2=et(1,6);
Vf3=et(1,7);
Vf4=et(1,8);
Tm1=et(1,1);
Tm2=et(1,2);
Tm3=et(1,3);
Tm4=et(1,4);

%Extra
if t>1 && t<1.2
     Pshsc=-1000/900;
     Qshsc=900/900;
end
% if t>1 && t<12
%     Pload7=500/900;
% end
% if t>3 && t<13
%     Ppv=0.2;
% end
% %Generator1
% phi1=0.752279097307610;
% theta01=phi1+pi/2;
% Te1=psid1*iq1-psiq1*id1;
% %Generator2
% phi2=0.562954707645971;
% theta02=phi2+pi/2;
% Te2=psid2*iq2-psiq2*id2;
% %Generator3
% phi3=0.297897626122744;
% theta03=phi3+pi/2;
% Te3=psid3*iq3-psiq3*id3;
% %Generator4
% phi4=0.103540782371745;
% theta04=phi4+pi/2;
% Te4=psid4*iq4-psiq4*id4;

%Generator1
phasIgen1=((vx1*ix1+vy1*iy1)-j*(vy1*ix1-vx1*iy1))/(vx1+j*vy1);
phi1=angle(vx1+j*vy1+Ra*phasIgen1+j*Xq*phasIgen1);
phi1=0.752279097307610;
theta01=phi1+pi/2;
Te1=psid1*iq1-psiq1*id1;
%Generator2
phasIgen2=((vx2*ix2+vy2*iy2)-j*(vy2*ix2-vx2*iy2))/(vx2+j*vy2);
phi2=angle(vx2+j*vy2+Ra*phasIgen2+j*Xq*phasIgen2);
phi2=0.562954707645971;
theta02=phi2+pi/2;
Te2=psid2*iq2-psiq2*id2;
%Generator3
phasIgen3=((vx3*ix3+vy3*iy3)-j*(vy3*ix3-vx3*iy3))/(vx3+j*vy3);
phi3=angle(vx3+j*vy3+Ra*phasIgen3+j*Xq*phasIgen3);
phi3=0.297897626122744;
theta03=phi3+pi/2;
Te3=psid3*iq3-psiq3*id3;
%Generator4
phasIgen4=((vx4*ix4+vy4*iy4)-j*(vy4*ix4-vx4*iy4))/(vx4+j*vy4);
phi4=angle(vx4+j*vy4+Ra*phasIgen4+j*Xq*phasIgen4);
phi4=0.103540782371745;
theta04=phi4+pi/2;
Te4=psid4*iq4-psiq4*id4;


kdd=0.9;
%%%Differential Equations
%Generator1
f1=wN*(-Rf*ifd1+Vf1); %psif (INPUT=Vf1)
f2=wN*(-Rq1*iq11); %psiq1
f3=(Tm1-Te1-kdd*omega1)/(2*HH12); %omega
f4=omega1; %deltas
%Generator2
f5=wN*(-Rf*ifd2+Vf2); %psif (INPUT=Vf2)
f6=wN*(-Rq1*iq12); %psiq1
f7=(Tm2-Te2-kdd*omega2)/(2*HH12); %omega
f8=omega2; %deltas
%Generator3
f9=wN*(-Rf*ifd3+Vf3); %psif (INPUT=Vf3)
f10=wN*(-Rq1*iq13); %psiq1
f11=(Tm3-Te3-kdd*omega3)/(2*HH34); %omega
f12=omega3; %deltas
%Generator4
f13=wN*(-Rf*ifd4+Vf4); %psif (INPUT=Vf4)
f14=wN*(-Rq1*iq14); %psiq1
f15=(Tm4-Te4-kdd*omega4)/(2*HH34); %omega
f16=omega4; %deltas

%%%Algebraic equations

%%Generator1
%flux-current
f17=Ldd*id1+Ldf*ifd1-psid1; 
f18=Lqq*iq1+Lqq1*iq11-psiq1; 
f19=Lff*ifd1+Ldf*id1-psif1; 
f20=Lq1q1*iq11+Lqq1*iq1-psiq11;
%stator park
f21=vd1+Ra*id1+psiq1; 
f22=vq1+Ra*iq1-psid1;
%dqxy
f23=vd1-vx1*cos(theta01)-vy1*sin(theta01); 
f24=vq1-vx1*sin(theta01)+vy1*cos(theta01); 
f25=id1-ix1*cos(theta01)-iy1*sin(theta01); 
f26=iq1-ix1*sin(theta01)+iy1*cos(theta01);

%%Generator2
%flux-current
f27=Ldd*id2+Ldf*ifd2-psid2; 
f28=Lqq*iq2+Lqq1*iq12-psiq2; 
f29=Lff*ifd2+Ldf*id2-psif2; 
f30=Lq1q1*iq12+Lqq1*iq2-psiq12;
%stator park
f31=vd2+Ra*id2+psiq2; 
f32=vq2+Ra*iq2-psid2;
%dqxy
f33=vd2-vx2*cos(theta02)-vy2*sin(theta02); 
f34=vq2-vx2*sin(theta02)+vy2*cos(theta02); 
f35=id2-ix2*cos(theta02)-iy2*sin(theta02); 
f36=iq2-ix2*sin(theta02)+iy2*cos(theta02);

%%Generator3
%flux-current
f37=Ldd*id3+Ldf*ifd3-psid3; 
f38=Lqq*iq3+Lqq1*iq13-psiq3; 
f39=Lff*ifd3+Ldf*id3-psif3; 
f40=Lq1q1*iq13+Lqq1*iq3-psiq13;
%stator park
f41=vd3+Ra*id3+psiq3; 
f42=vq3+Ra*iq3-psid3;
%dqxy
f43=vd3-vx3*cos(theta03)-vy3*sin(theta03); 
f44=vq3-vx3*sin(theta03)+vy3*cos(theta03); 
f45=id3-ixg3*cos(theta03)-iyg3*sin(theta03); 
f46=iq3-ixg3*sin(theta03)+iyg3*cos(theta03);

%%Generator4
%flux-current
f47=Ldd*id4+Ldf*ifd4-psid4; 
f48=Lqq*iq4+Lqq1*iq14-psiq4; 
f49=Lff*ifd4+Ldf*id4-psif4; 
f50=Lq1q1*iq14+Lqq1*iq4-psiq14;
%stator park
f51=vd4+Ra*id4+psiq4; 
f52=vq4+Ra*iq4-psid4;
%dqxy
f53=vd4-vx4*cos(theta04)-vy4*sin(theta04); 
f54=vq4-vx4*sin(theta04)+vy4*cos(theta04); 
f55=id4-ix4*cos(theta04)-iy4*sin(theta04); 
f56=iq4-ix4*sin(theta04)+iy4*cos(theta04);

%%Transformers
%1-5
f57=vx1-vx5+Xt*iy1; 
f58=vy1-vy5-Xt*ix1;
%2-6
f59=vx2-vx6+Xt*iy2; 
f60=vy2-vy6-Xt*ix2;
%3-11
f61=vx3-vx11+Xt*iy3; 
f62=vy3-vy11-Xt*ix3;
%4-10
f63=vx4-vx10+Xt*iy4; 
f64=vy4-vy10-Xt*ix4;

%%Transmission Lines
%5-6
f65=vx5-vx6-Re56*ix1+Xe56*iy1; 
f66=vy5-vy6-Re56*iy1-Xe56*ix1;
%6-7
f67=vx6-vx7-Re67*ix1+Xe67*iy1; 
f68=vy6-vy7-Re67*iy1-Xe67*ix1;
%7-8
f69=vx7-vx8-Re78*ix5+Xe78*iy5; 
f70=vy7-vy8-Re78*iy5-Xe78*ix5;
%8-9
f71=vx8-vx9-Re89*ix8+Xe89*iy8; 
f72=vy8-vy9-Re89*iy8-Xe89*ix8;
%10-9
f73=vx10-vx9-Re910*ix6+Xe910*iy6; 
f74=vy10-vy9-Re910*iy6-Xe910*ix6;
%11-10
f75=vx11-vx10-Re1011*ix3+Xe1011*iy3; 
f76=vy11-vy10-Re1011*iy3-Xe1011*ix3;

%%KCL
%6
f77=ix1+ix2-ix5;
f78=iy1+iy2-iy5;
%7
f79=ix7+ix8-ix5;
f80=iy7+iy8-iy5;
%9
f81=ix6+ix8-ix9;
f82=iy6+iy8-iy9;
%10
f83=ix3+ix4-ix6;
f84=iy3+iy4-iy6;

%%Loads
%7
f85=-Pload7+vx7*ix7+vy7*iy7;
f86=-Qload7-vx7*iy7+vy7*ix7;
%9
f87=-Pload9+vx9*ix9+vy9*iy9;
f88=-Qload9-vx9*iy9+vy9*ix9;
%3 Short-circuit
f89=-Pshsc+vx3*ixshsc+vy3*iyshsc;
f90=-Qshsc-vx3*iyshsc+vy3*ixshsc;
f91=-ixg3+ixshsc+ix3;
f92=-iyg3+iyshsc+iy3;


%sum
%Differential
ff1=[f1 f2 f3 f4 f5 f6 f7 f8 f9 f10 f11 f12 f13 f14 f15 f16]';
%g1
ff2=[f17 f18 f19 f20 f21 f22 f23 f24 f25 f26]';
%g1
ff3=[f27 f28 f29 f30 f31 f32 f33 f34 f35 f36]';
%g1
ff4=[f37 f38 f39 f40 f41 f42 f43 f44 f45 f46]';
%g1
ff5=[f47 f48 f49 f50 f51 f52 f53 f54 f55 f56]';
%Other
ff6=[f57 f58 f59 f60 f61 f62 f63 f64 f65 f66 f67 f68 f69 f70 f71 f72 f73 f74 f75 f76 f77 f78 f79 f80 f81 f82 f83 f84 f85 f86 f87 f88]';
%shsc
ff7=[f89 f90 f91 f92]';
f=[ff1;ff2;ff3;ff4;ff5;ff6;ff7];
end

%initial values
function [ssi]=ini
%Per Unit 900MVA 20KV primary 230KV secondary
j=sqrt(-1);
%parameters
pp=para;
Xd=pp(1,1); Xq=pp(1,2); Xl=pp(1,3); Xpd=pp(1,4); Xpq=pp(1,5);
Tpdo_s=pp(1,6); Tpqo_s=pp(1,7); Ra=pp(1,8); fN=pp(1,9); H12=pp(1,10);
H34=pp(1,11); Xt=pp(1,12); Pload7=pp(1,13); Qload7=pp(1,14); Pload9=pp(1,15);
Qload9=pp(1,16); Re56=pp(1,17); Xe56=pp(1,18); Re67=pp(1,19); Xe67=pp(1,20);
Re1011=pp(1,21); Xe1011=pp(1,22); Re910=pp(1,23); Xe910=pp(1,24); Re78=pp(1,25);
Xe78=pp(1,26); Re89=pp(1,27); Xe89=pp(1,28); wN=pp(1,29); TB=pp(1,30);
Ldd=pp(1,31); Ldf=pp(1,32); Lff=pp(1,33); Tpdo=pp(1,34); Rf=pp(1,35);
Lqq=pp(1,36); Lqq1=pp(1,37); Lq1q1=pp(1,38); Tpqo=pp(1,39); Rq1=pp(1,40);
%%%Generator Values
%1
Pgen1=700/900;
Qgen1=185/900;
Vgen1=1.03;
% deltagen1=20.2;
% phasVgen1=Vgen1*(cos(deg2rad(deltagen1))+j*sin(deg2rad(deltagen1)));
phasVgen1=Vgen1*(cos(0)+j*sin(0));
vx1=real(phasVgen1); 
vy1=imag(phasVgen1);
It1=(sqrt(Pgen1^2+Qgen1^2))/Vgen1; %Terminal Current
% phi1=acos(Pgen1/(Vgen1*It1)); %Power Factor angle
% deltas1=atan((Xq*It1*cos(phi1)-Ra*It1*sin(phi1))/(Vgen1+Ra*It1*cos(phi1)+Xq*It1*sin(phi1))); %Rotor Angle
phasIgen1=(Pgen1+j*Qgen1)/(vx1+j*vy1); %PV current (phasor)
phasIgen1=conj(phasIgen1);
deltas1=angle(phasVgen1+Ra*phasIgen1+j*Xq*phasIgen1);
ix1=real(phasIgen1);
iy1=imag(phasIgen1);
a=sqrt(ix1^2+iy1^2);
% ix11=It1*cos(phi1);
% iy11=It1*sin(phi1);
S1=phasVgen1*(ix1-j*iy1);
% S2=phasVgen1*(ix11-j*iy11);
theta01=deltas1+pi/2;
vd1=Vgen1*cos(theta01);
vq1=Vgen1*sin(theta01);
id1=abs(phasIgen1)*cos(theta01+angle(phasIgen1));
iq1=abs(phasIgen1)*sin(theta01+angle(phasIgen1));
psid1=vq1+Ra*iq1;
psiq1=Lqq*iq1;
psiq11=Lqq1*iq1;
ifd1=(psid1-Ldd*id1)/Ldf;
psif1=Lff*ifd1+Ldf*id1;
Vf1=Rf*ifd1;
Tm1=psid1*iq1-psiq1*id1;
Tm11=Pgen1+Ra*It1^2;
iq11=0;
omega1=0;
%2
Pgen2=700/900;
Qgen2=235/900;
Vgen2=1.01;
% deltagen2=10.5;
% phasVgen2=Vgen2*(cos(deg2rad(deltagen2))+j*sin(deg2rad(deltagen2)));
phasVgen2=Vgen2*(cos(-0.1707)+j*sin(-0.1707));
vx2=real(phasVgen2); 
vy2=imag(phasVgen2);
It2=(sqrt(Pgen2^2+Qgen2^2))/Vgen2; %Terminal Current
% phi2=acos(Pgen2/(Vgen2*It2)); %Power Factor angle
% deltas2=atan((Xq*It2*cos(phi2)-Ra*It2*sin(phi2))/(Vgen2+Ra*It2*cos(phi2)+Xq*It2*sin(phi2))); %Rotor Angle
phasIgen2=(Pgen2+j*Qgen2)/(vx2+j*vy2); %PV current (phasor)
phasIgen2=conj(phasIgen2);
deltas2=angle(phasVgen2+Ra*phasIgen2+j*Xq*phasIgen2);
ix2=real(phasIgen2);
iy2=imag(phasIgen2);
theta02=deltas2+pi/2;
vd2=Vgen2*cos(theta02);
vq2=Vgen2*sin(theta02);
id2=abs(phasIgen2)*cos(theta02-angle(phasIgen2));
iq2=abs(phasIgen2)*sin(theta02-angle(phasIgen2));
psid2=vq2+Ra*iq2;
psiq2=Lqq*iq2;
psiq12=Lqq1*iq2;
ifd2=(psid2-Ldd*id2)/Ldf;
psif2=Lff*ifd2+Ldf*id2;
Vf2=Rf*ifd2;
Tm2=psid2*iq2-psiq2*id2;
Tm22=Pgen2+Ra*It2^2;
iq12=0;
omega2=0;
%3
Pgen3=719/900;
Qgen3=176/900;
Vgen3=1.03;
% deltagen3=-6.8;
% phasVgen3=Vgen3*(cos(deg2rad(deltagen3))+j*sin(deg2rad(deltagen3)));
phasVgen3=Vgen3*(cos(-0.4738)+j*sin(-0.4738));
vx3=real(phasVgen3); 
vy3=imag(phasVgen3);
It3=(sqrt(Pgen3^2+Qgen3^2))/Vgen3; %Terminal Current
% phi3=acos(Pgen3/(Vgen3*It3)); %Power Factor angle
% deltas3=atan((Xq*It3*cos(phi3)-Ra*It3*sin(phi3))/(Vgen3+Ra*It3*cos(phi3)+Xq*It3*sin(phi3))); %Rotor Angle
phasIgen3=(Pgen3+j*Qgen3)/(vx3+j*vy3); %PV current (phasor)
phasIgen3=conj(phasIgen3);
deltas3=angle(phasVgen3+Ra*phasIgen3+j*Xq*phasIgen3);
ix3=real(phasIgen3);
iy3=imag(phasIgen3);
theta03=deltas3+pi/2;
vd3=Vgen3*cos(theta03);
vq3=Vgen3*sin(theta03);
id3=abs(phasIgen3)*cos(theta03-angle(phasIgen3));
iq3=abs(phasIgen3)*sin(theta03-angle(phasIgen3));
psid3=vq3+Ra*iq3;
psiq3=Lqq*iq3;
psiq13=Lqq1*iq3;
ifd3=(psid3-Ldd*id3)/Ldf;
psif3=Lff*ifd3+Ldf*id3;
Vf3=Rf*ifd3;
Tm3=psid3*iq3-psiq3*id3;
Tm33=Pgen3+Ra*It3^2;
iq13=0;
omega3=0;
%4
Pgen4=700/900;
Qgen4=202/900;
Vgen4=1.01;
% deltagen4=-17;
% phasVgen4=Vgen4*(cos(deg2rad(deltagen4))+j*sin(deg2rad(deltagen4)));
phasVgen4=Vgen4*(cos(-0.6518)+j*sin(-0.6518));
vx4=real(phasVgen4); 
vy4=imag(phasVgen4);
It4=(sqrt(Pgen4^2+Qgen4^2))/Vgen4; %Terminal Current
% phi4=acos(Pgen4/(Vgen4*It4)); %Power Factor angle
% deltas4=atan((Xq*It4*cos(phi4)-Ra*It4*sin(phi4))/(Vgen4+Ra*It4*cos(phi4)+Xq*It4*sin(phi4))); %Rotor Angle
phasIgen4=(Pgen4+j*Qgen4)/(vx4+j*vy4); %PV current (phasor)
phasIgen4=conj(phasIgen4);
deltas4=angle(phasVgen4+Ra*phasIgen4+j*Xq*phasIgen4);
ix4=real(phasIgen4);
iy4=imag(phasIgen4);
theta04=deltas4+pi/2;
vd4=Vgen4*cos(theta04);
vq4=Vgen4*sin(theta04);
id4=abs(phasIgen4)*cos(theta04-angle(phasIgen4));
iq4=abs(phasIgen4)*sin(theta04-angle(phasIgen4));
psid4=vq4+Ra*iq4;
psiq4=Lqq*iq4;
psiq14=Lqq1*iq4;
ifd4=(psid4-Ldd*id4)/Ldf;
psif4=Lff*ifd4+Ldf*id4;
Vf4=Rf*ifd4;
Tm4=psid4*iq4-psiq4*id4;
Tm44=Pgen4+Ra*It4^2;
iq14=0;
omega4=0;

%%%Other
%5
Vgen5=1.005921 ;
phasVgen5=Vgen5*(cos(-0.1128836)+j*sin(-0.1128836));
vx5=real(phasVgen5); 
vy5=imag(phasVgen5);
%6
Vgen6=0.9773796 ;
phasVgen6=Vgen6*(cos(-0.2891676)+j*sin(-0.2891676));
vx6=real(phasVgen6); 
vy6=imag(phasVgen6);
%7
Vgen7=0.9600698 ;
phasVgen7=Vgen7*(cos(-0.4361968)+j*sin(-0.4361968));
vx7=real(phasVgen7); 
vy7=imag(phasVgen7);
%8
Vgen8=0.9475377 ;
phasVgen8=Vgen8*(cos(-0.6788807)+j*sin(-0.6788807));
vx8=real(phasVgen8); 
vy8=imag(phasVgen8);
%9
Vgen9=0.9704016 ;
phasVgen9=Vgen9*(cos(-0.9168022)+j*sin(-0.9168022));
vx9=real(phasVgen9); 
vy9=imag(phasVgen9);
%10
Vgen10=0.9827057 ;
phasVgen10=Vgen10*(cos(-0.7697021)+j*sin(-0.7697021));
vx10=real(phasVgen10); 
vy10=imag(phasVgen10);
%11
Vgen11=1.007725 ;
phasVgen11=Vgen11*(cos(-0.5895673)+j*sin(-0.5895673));
vx11=real(phasVgen11); 
vy11=imag(phasVgen11);

% ix5=(ix1+ix2);
% iy5=(iy1+iy2);
% ix6=(ix3+ix4);
% iy6=(iy3+iy4);
% phasIgen5=ix5+j*iy5;
phasIgen5=phasIgen1+phasIgen2;
ix5=real(phasIgen5); 
iy5=imag(phasIgen5);
phasIgen6=phasIgen3+phasIgen4;
ix6=real(phasIgen6); 
iy6=imag(phasIgen6);
phasIgen7=(Pload7+j*Qload7)/(vx7+j*vy7); %PV current (phasor)
phasIgen7=conj(phasIgen7);
ix7=real(phasIgen7); 
iy7=imag(phasIgen7);
% phasIgen8=(Pload8-j*Qload8)/(vx8-j*vy8); %PV current (phasor)
phasIgen8=phasIgen5+phasIgen7;
ix8=real(phasIgen8); 
iy8=imag(phasIgen8);
phasIgen9=(Pload9+j*Qload9)/(vx9+j*vy9); %PV current (phasor)
phasIgen9=conj(phasIgen9);
ix9=real(phasIgen9); 
iy9=imag(phasIgen9);



%%%Output
%Generator1
ssi1=[psif1 psiq11 omega1 deltas1 ]'; %4
%Generator2
ssi2=[psif2 psiq12 omega2 deltas2 ]'; %4
%Generator3
ssi3=[psif3 psiq13 omega3 deltas3 ]'; %4
%Generator4
ssi4=[psif4 psiq14 omega4 deltas4 ]'; %4
%other Generator
ssi5=[ifd1 iq11 vd1 vq1 id1 iq1 psid1 psiq1 ifd2 iq12 vd2 vq2 id2 iq2 psid2 psiq2 ifd3 iq13 vd3 vq3 id3 iq3 psid3 psiq3 ifd4 iq14 vd4 vq4 id4 iq4 psid4 psiq4]'; %32
%Voltage
ssi6=[vx1 vy1 vx2 vy2 vx3 vy3 vx4 vy4 vx5 vy5 vx6 vy6 vx7 vy7 vx8 vy8 vx9 vy9 vx10 vy10 vx11 vy11]'; %22
%Current
ssi7=[ix1 iy1 ix2 iy2 ix3 iy3 ix4 iy4 ix5 iy5 ix6 iy6 ix7 iy7 ix8 iy8 ix9 iy9]'; %18
ssi=[ssi1;ssi2;ssi3;ssi4;ssi5;ssi6;ssi7];
% ssi=[0.905070902653913;0.733652719078255;0.00174560535242228;-53.0944791339272;0.837979776944666;0.655871641541966;0.0174835630424538;6.16708841297472;0.993110830879738;0.507359436587652;0.0218829252023258;4.94298230064510;0.972658264153199;0.376599710110103;0.0299573994018910;9.72908599535070;1.26948677680599;-1.48679723541597e-05;-0.829512622990016;0.610746193893419;-0.788449914434362;0.489121205726255;0.611968996907735;0.831483747776102;0.802557297281113;-1.44719520641293e-06;-0.742491422561544;0.684815789811108;-0.332323756509437;0.437249648673899;0.685908913932793;0.743322231952818;0.816505609209263;8.95625663163273e-06;-0.574375323707403;0.855122319936946;-0.250245047190045;0.338227942317901;0.855967889792741;0.575000936325379;0.556942469681813;1.98605870489200e-05;-0.426833355597590;0.915391355453436;0.0138394474351487;0.251040568293193;0.916018956874169;0.426798756979002;1.01273599805209;-0.188329241766957;0.975394546415127;-0.262434101015294;0.986044553353040;-0.298089472040079;0.954604658127539;-0.329936348059411;0.976517763650269;-0.322710520149494;0.968241703276923;-0.344503762639886;0.964931279127585;-0.353221059636043;0.960596867691536;-0.367361981786251;0.958914012158546;-0.373526431296326;0.960561435247690;-0.367176205660000;0.965052124863901;-0.357606564121284;0.895875189216908;-0.241454896012126;0.547131077497281;-0.0476856209213608;0.396780613874703;-0.139949523260932;0.248265717337259;0.0397118474676737;1.44300626671419;-0.289140516933486;0.645046331211963;-0.100237675793259;0.816002795093051;-0.183555310796605;0.627003471621139;-0.105585206136882;1.27204980283310;-0.205822881930140;0.396482927787140;-0.138983297959703;0;0];
ssi=[0.905043779757177;0.733722524115933;3.69308438977887e-06;1.02147299460645;0.837432667426397;0.655973770115201;1.81354191779677e-05;0.793705459568758;0.992949534451635;0.507509655236909;-6.78158237972163e-06;0.452913395858058;0.973321016900431;0.376774521323109;-1.23953602281832e-05;-0.0602348340154869;1.26961372500876;-3.99104505496737e-05;-0.829609319714081;0.610674859207719;-0.788602277661112;0.489200406520035;0.611897860224019;0.831580825408234;0.801325533053331;4.15858612157834e-07;-0.742608261585254;0.684594292016833;-0.331351818115388;0.437315304319191;0.685687580277631;0.743436641130542;0.816137985376778;5.25607105433953e-05;-0.574515018358379;0.855058392905746;-0.249953725369625;0.338271212709405;0.855904070937519;0.575139902671803;0.558592248651113;9.05107108709648e-05;-0.426977427328431;0.915608962052713;0.0124939036682487;0.251064956766444;0.916236624444629;0.426946192569261;1.01274998752626;-0.188448586596930;0.975269585235510;-0.262651120217023;0.986024443872969;-0.298241777225941;0.954835989981744;-0.330057157284667;0.976523184040961;-0.322854155625336;0.968245290642937;-0.344651303738649;0.964934133283728;-0.353370162983974;0.960608021728203;-0.367508872500123;0.958934406027901;-0.373671886914067;0.960592394561647;-0.367321514042501;0.965075694046697;-0.357752248066225;0.896037126856041;-0.241512023235349;0.546667890144173;-0.0468286306171527;0.396736472268559;-0.139658332175148;0.248429045052228;0.0383760305326882;1.44270501700021;-0.288340653852502;0.645165517320788;-0.101282301642460;0.815935117596906;-0.183656073762700;0.626769899403309;-0.104684580089802;1.27193541672410;-0.205966881732262;0.396736472268559;-0.139658332175148;0;0];
% ssi(3,1)=0;
% ssi(7,1)=0;
% ssi(11,1)=0;
% ssi(15,1)=0;
% ssi(4,1)=1.018643899376857;
% ssi(8,1)=0.788267725737706;
% ssi(12,1)=0.453544908445960;
% ssi(16,1)=-0.061024033441272;
end

function pp=para
%%%DATA%%%
%Per Unit 900MVA 20KV primary 230KV secondary
%Generator parameters
Xd=1.8;
Xq=1.7;
Xl=0.2;
Xpd=0.3;
Xpq=0.55;
Tpdo_s=8;
Tpqo_s=0.4;
Ra=0.0025;
fN=50;
H12=6.5;
H34=6.175;
%Transformers
Xt=0.15;
%Loads
Qc7=200/900;
Qc9=350/900;
Pload7=767/900;
Qload7=100/900-Qc7;
Pload9=1167/900;
Qload9=100/900-Qc9;
% Qc7=0/900;
% Qc9=0/900;
% Pload7=0/900;
% Qload7=0/900-Qc7;
% Pload9=0/900;
% Qload9=0/900-Qc9;
%Transmission Lines
Re56=0.0001*25;
Xe56=0.001*25;
Re67=0.0001*10;
Xe67=0.001*10;
Re1011=0.0001*25;
Xe1011=0.001*25;
Re910=0.0001*10;
Xe910=0.001*10;
Re78=0.0001*10;
Xe78=0.001*10;
Re89=0.0001*10;
Xe89=0.001*10;
%Voltage and Power
%%%Parameters%%%
wN=2*pi*fN;
TB=1/(2*pi*fN);
Ldd=Xd;
Ldf=Xd-Xl;
Lff=Ldf^2/(Xd-Xpd);
Tpdo=Tpdo_s/TB;
Rf=Lff/Tpdo;
Lqq=Xq;
Lqq1=Xq-Xl;
Lq1q1=Lqq1^2/(Xq-Xpq);
Tpqo=Tpqo_s/TB;
Rq1=Lq1q1/Tpqo;
%%% 1  2  3  4   5   6      7      8  9  10  11  12 13     14     15     16     17   18   19   20   21     22     23    24    25   26   27   28   29 30 31  32  33  34   35 36  37   38    39   40          
pp=[Xd Xq Xl Xpd Xpq Tpdo_s Tpqo_s Ra fN H12 H34 Xt Pload7 Qload7 Pload9 Qload9 Re56 Xe56 Re67 Xe67 Re1011 Xe1011 Re910 Xe910 Re78 Xe78 Re89 Xe89 wN TB Ldd Ldf Lff Tpdo Rf Lqq Lqq1 Lq1q1 Tpqo Rq1];
end