%% %% This file solves synchronous machine using the traditional method (reference)
%%%% Updated Controllers
%%%% Please set the parameters based on the guidance provided in the Readme.md file 
%both torque and AVR
%[0 1.1]     [0 0.0115]
clearvars -except net fuzzyAVR
global net
global fuzzyAVR;
clc
tic
% This parameter chooses the predictor method. It can take values 4 and 5. 4 is for the first-order Adams-Bashford, and 5 is for the second-order Adams-Bashford. 
cp=5; %Predictor Choice
% This parameter chooses the corrector method. It can take values 2 and 3. 2 is for the second-order BDF, and 3 is for the second-order Adams-Moulton. 
cc=3; %Corrector Choice
%This parameters chooses the quantization method
%1 uniform quantization, 2 Nonuniform quantization 
cq=1; %Quantization method
save=0;

t=0; %Starting point of the integration
h=0.001; %The first time-step size
holdd=h;
hmin=0.001; %The minimum acceptable time-step size
hmax=1; %The maximum acceptable time-step size

tsim=20; %Time of the simulation
%Initial state variables values
[ssi]=ini;
xx=ssi;
dae=[8 13];
cn=2;
T=[0.02 0.006]; %Time of sampling by the controller
Tm=0.5013;
Vf=8.638587815812726e-04;
G=0.05; %Gain of the controller
u=1; %Input 

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
avrov=[1;0;0;Vf;Vf;Vf;0;0];
gvrov=Tm;


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
    flag=0;
    
       %manual time control
    if t>=1 && tpre1<1
        t=1;
        h=t-tpre1;
        flag=0;
    end
    if t>=1.2 && tpre1<1.2
        t=1.2;
        h=t-tpre1;
        flag=0;
    end
    
    %Setting old solutions
    if k>2
        xpre3=xhist(1:s,k-2);
    else
        xpre3=0;
    end
    if k>1
        xpre2=xhist(1:s,k-1);
        avrovpre=avrov;
        gvrovpre=gvrov;
    else
        xpre2=xpre1;
        avrovpre=avrov;
        gvrovpre=gvrov;
    end
    xpre1=xhist(1:s,k);
    %Setting new and old controllers's signal
    if t<min(T)
        et=[Tm,Vf];
        etpre=[Tm,Vf];
        etpre2=[Tm,Vf];
        avrov=[1;0;0;Vf;Vf;Vf;0;0];
        gvrov=Tm;
    else
        tgpre=tg;
        [t,tg,mm,temp]=detert(tpre1,t,tg,T,cn,tsim); %Finding the event if there is any
        if sum(mm)==0
            etpre2=etpre;
            etpre=et;
            et=et;
            tg=tg;
        else
            etpre2=etpre;
            etpre=et;
            jjj=0;
            mmsum=sum(mm);
            while jjj<cn
                jjj=jjj+1;
                if mm(1,jjj)~=0
                    [et,avrov,gvrov]=contex(mm,et,T,cn,avrov,gvrov,xx,net,jjj); %Calculating the controllers' signal in case of any event
                end
            end
            et=quz(et,cq);
            flagg=1;
        end
    end
    avrovhist(:,k)=avrov;
    %Summation of the controllers' signal
    if flag==0
        k=k+1;
        times(k)=t;
    end

    [alpha0,alpha1,alpha2,betha1,betha2]=coefs(cp,h,holdd); %calculating the coefficients of the predictor
    xx0=predictor(xx,xpre1,et,etpre,t,tpre1,h,betha1,betha2,dae); %predictor
    [alpha0,alpha1,alpha2,betha1,betha2]=coefs(cc,h,holdd); %calculating the coefficients of the corrector
    [iterations,xx]=NewtonNL(h,xx0,xpre1,xpre2,alpha1,alpha2,alpha0,betha1,betha2,tpre1,et,etpre,t,dae); %The main function (Newton solver)
    [dn,r]=lote(cp,cc,xx0,xx,h,holdd); %Error estimate
    iterationsc=iterations+iterationsc; %Counting the number of Newton iterations

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
                 hhist(4,1)=et(1,1);
                 hhist(5,1)=et(1,2);
            end
            hhist(1,k)=h;
            hhist(2,k)=dn;
            hhist(3,k)=r;
             hhist(4,k)=et(1,1);
             hhist(5,k)=et(1,2);
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
                avrov=avrovpre;
                gvrov=gvrovpre;
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
             hhist(4,1)=et(1,1);
             hhist(5,1)=et(1,2);
        end
        hhist(1,k)=h;
        hhist(2,k)=dn;
        hhist(3,k)=r;
         hhist(4,k)=et(1,1);
         hhist(5,k)=et(1,2);
        holdd=hpre;
    end


end
toc
%Plotting the results
figure(1)
plot(times,sqrt(xhist(14,:).^2+xhist(15,:).^2),'--','color',[0.4660 0.6740 0.1880])
% plot(times,xhist(9,:))
hold on
% plot(times,xhist(10,:))
xlabel('Time','FontSize', 24)
ylabel('Voltage','FontSize', 24)
set(gca,'FontSize',18)

figure(2)
% plot(times,xhist(13,:))
plot(times,sqrt((xhist(18,:).^2)+(xhist(19,:).^2)),'--','color',[0.4660 0.6740 0.1880])
hold on
% plot(times,xhist(14,:))
xlabel('Time','FontSize', 24)
ylabel('Current','FontSize', 24)
set(gca,'FontSize',18)

figure(5)
plot(times,xhist(12,:),'--')
hold on
plot(times,xhist(13,:),'--')
xlabel('Time','FontSize', 24)
ylabel('vd & vq','FontSize', 24)
set(gca,'FontSize',18)

figure(3)
plot(times,xhist(3,:),'--','color',[0.4660 0.6740 0.1880])
hold on
xlabel('Time','FontSize', 24)
ylabel('Speed deviation','FontSize', 24)
set(gca,'FontSize',18)

figure(4)
plot(times,xhist(4,:),'--','color',[0.4660 0.6740 0.1880])
hold on
xlabel('Time','FontSize', 24)
ylabel('Rotor angle','FontSize', 24)
set(gca,'FontSize',18)

figure(6)
stairs(times,hhist(4,:),'--','color',[0.4660 0.6740 0.1880])
hold on
xlabel('Time (s)','FontSize', 24)
ylabel('Governor output (Per unit)','FontSize', 24)
set(gca,'FontSize',18)
legend('Without quantization (SRM)','With quantization (SRM)','Location','southeast')
set(gca, 'FontName', 'Times New Roman')
xlim([0 tsim])
figure_property.units = 'inches';
figure_property.format = 'pdf';
figure_property.Preview= 'none';
figure_property.Width= '16'; % Figure width on canvas
figure_property.Height= '9'; % Figure height on canvas
figure_property.Units= 'inches';
figure_property.Color= 'rgb';
figure_property.Background= 'w';
figure_property.FixedfontSize= '12';
figure_property.ScaledfontSize= 'auto';
figure_property.FontMode= 'scaled';
figure_property.FontSizeMin= '12';
figure_property.FixedLineWidth= '1';
figure_property.ScaledLineWidth= 'auto';
figure_property.LineMode= 'none';
figure_property.LineWidthMin= '0.1';
figure_property.FontName= 'Times New Roman';% Might want to change this to something that is available
figure_property.FontWeight= 'auto';
figure_property.FontAngle= 'auto';
figure_property.FontEncoding= 'latin1';
figure_property.PSLevel= '3';
figure_property.Renderer= 'painters';
figure_property.Resolution= '600';
figure_property.LineStyleMap= 'none';
figure_property.ApplyStyle= '0';
figure_property.Bounds= 'tight';
figure_property.LockAxes= 'off';
figure_property.LockAxesTicks= 'off';
figure_property.ShowUI= 'off';
figure_property.SeparateText= 'off';
chosen_figure=gcf;
set(chosen_figure,'PaperUnits','inches');
set(chosen_figure,'PaperPositionMode','auto');
set(chosen_figure,'PaperSize',[str2num(figure_property.Width) str2num(figure_property.Height)]); % Canvas Size
set(chosen_figure,'Units','inches');
if save==1
    hgexport(gcf,'SMIBSRMQuantizationTest.pdf',figure_property);
end


figure(7)
plot(times,hhist(5,:).*(1/0.000747681),'--','color',[0.4660 0.6740 0.1880])
hold on
xlabel('Time','FontSize', 24)
ylabel('AVR Signal','FontSize', 24)
set(gca,'FontSize',18)
% figure(77)
% plot(times,xhist(5,:).*(1/0.000747681))
% hold on
% xlabel('Time','FontSize', 24)
% ylabel('AVR Signal','FontSize', 24)
% set(gca,'FontSize',18)


figure (8)
plot(times,sqrt((xhist(14,:).^2)+(xhist(15,:).^2)),'--','color',[0.4660 0.6740 0.1880])
hold on
xlabel('Time','FontSize', 24)
ylabel('Vter','FontSize', 24)
set(gca,'FontSize',18)

figure(9)
plot(times,xhist(9,:),'--','color',[0.4660 0.6740 0.1880])
hold on
xlabel('Time','FontSize', 24)
ylabel('Frequency','FontSize', 24)
set(gca,'FontSize',18)

% figure(2)
% plot(xhist(1,:),xhist(3,:))

figure(11)
plot(times,hhist(1,:),'color',[0.4660 0.6740 0.1880])
hold on
xlabel('Time','FontSize', 24)
ylabel('Step Size','FontSize', 24)
set(gca,'FontSize',18)

figure(12)
plot(times,hhist(2,:),'color',[0.4660 0.6740 0.1880])
xlabel('Time','FontSize', 24)
ylabel('Error Estimate','FontSize', 24)
set(gca,'FontSize',18)
hold on

figure(20)
plot(times,xhist(14,:).*xhist(18,:)+xhist(15,:).*xhist(19,:),'--','color',[0.4660 0.6740 0.1880])
xlabel('Time','FontSize', 24)
ylabel('Active power (Per unit)','FontSize', 24)
set(gca,'FontSize',18)
hold on
figure(21)
plot(times,xhist(15,:).*xhist(18,:)-xhist(14,:).*xhist(19,:),'--','color',[0.4660 0.6740 0.1880])
xlabel('Time','FontSize', 24)
ylabel('Reactive power (Per unit)','FontSize', 24)
set(gca,'FontSize',18)
hold on

figure(112)
plot(times,xhist(6,:))
hold on
plot(times,xhist(20,:).*xhist(17,:)-xhist(21,:).*xhist(16,:),'--','color',[0.4660 0.6740 0.1880])
xlabel('Time','FontSize', 24)
ylabel('TM & Te','FontSize', 24)
set(gca,'FontSize',18)


function [t,tg,mm,temp]=detert(tpre1,t,tg,T,cn,tsim)
temp=zeros(1,cn);
temp=tg;
cni=0;
mm=zeros(1,cn);
while cni<cn
    cni=cni+1;
    while temp(1,cni)<=t+0.0000001
        if tpre1<tg(1,cni) && tg(1,cni)<=t+0.0000001
            temp(1,cni)=temp(1,cni)+T(1,cni);
            mm(1,cni)=mm(1,cni)+1;
        else
            mm(1,cni)=0;
        end
    end
end
if sum(mm)==0
    t=t;
    tg=tg;
    minl=0;
else
    temp=round(temp,6); 
    tg=temp;
end
end

%This function updates the controllers' signal in case of the event 
function [etj,avrov,gvrov]=contex(mm,etj,T,cn,avrov,gvrov,xx,net,jjj)
global fuzzyAVR
contj=jjj;
    if contj==1
        [etj(1,contj),gvrov]=torque(etj(1,contj),T(1,1),xx,gvrov);
        etj(1,contj)=max(etj(1,contj),0.1); %upper limnit anti wind up
        etj(1,contj)=min(etj(1,contj),0.9); %lower limnit anti wind up
    elseif contj==2
        Vter=sqrt(xx(14,1)^2+xx(15,1)^2);
        ntime=tic;
         etj(1,contj)=evalfis(fuzzyAVR,Vter);
%        [etj(1,contj)]=sim(net,Vter);
%         etj(1,contj)=net(Vter);
%         etj(1,contj)=ANN2(net,Vter);
%         system('custom.exe CustomI.txt')
        finalntime=toc(ntime)
%         [y]=sim(net,[Vter avrov']'); %Controller's equation
%         etj(1,contj)=y(1,1);
%         avrov=y(2:end,1);
        etj(1,contj)=max(etj(1,contj),0); %upper limnit anti wind up
        etj(1,contj)=min(etj(1,contj),1); %lower limnit anti wind up
    end
end

function Output=ANN(Input)
%#codegen
    Output=net(Input);
end

function [et2,gvrov]=torque(et2,T,xg,gvrov)
omega=xg(3,1);
% Tm=xg(6,1);
wS=1+omega;
% Pe=xg(14,1)*xg(18,1)+xg(15,1)*xg(19,1);
Te=xg(20,1)*xg(17,1)-xg(21,1)*xg(16,1);
%Data
kp=12;
ki=0.2;
Droop=1;
zo=0.5013;
%initialization (old values)
T2=T;
x5=gvrov;
%equations
x3=zo-wS+1-(Te*(Droop));
x4=x3*kp;
x5=x5+T2*x3*ki;
x6=x5+x4;
%assigning old values for the next sample
gvrov=x5;
et2=x6;
end



function [et1,avrov]=avr(avrov,T1,xx)
%Data
Tmm=0.2;
Te=0.12;
vo=1;
kp=0.003;
ki=0.004;
Ge=1;
kf=0.8;
Tf=0.9;
Vter=sqrt(xx(14,1)^2+xx(15,1)^2);
x62=xx(5,1);
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
x5=x4+x3; %Controller output
%x62=(x61*(Te-T1)+T1*Ge*x5)/Te;
x7=(kf/T1)*(x62-x61);
x8=(Tf*x8-T1*x8+T1*x7)/Tf;
et1=x5;
%assigning old values for the next sample
avrov=[x1 x2 x3 x4 x5 x62 x7 x8]';
end

function [iterations,xx]=NewtonNL(h,xx,xpre1,xpre2,alpha1,alpha2,alpha0,betha1,betha2,tpre1,et,etpre,t,dae)
NTOL=1e-4; %Tolerance for the Newton solver
iterations=0;
a=1;
    while (iterations<7) && (a>NTOL)
    xxpre=xx;
    iterations=iterations+1;
    [J]=Jacob(h,xx,xpre1,xpre2,alpha1,alpha2,alpha0,betha1,betha2,tpre1,et,etpre,t,dae); %calculating the Jacobian
    yy=xx-J\func(h,xx,xpre1,xpre2,alpha1,alpha2,alpha0,betha1,betha2,tpre1,et,etpre,t,dae); %Newton
    xx=yy;
%     a=abs(xx-xxpre);
    a=abs((xx-xxpre)./xx);
    a=max(a);
    if iterations==30
        disp('Did not converged')
        t=t
    end
    end
end

% This function computes the Jacobian
function [J]=Jacob(h,xx,xpre1,xpre2,alpha1,alpha2,alpha0,betha1,betha2,tpre1,et,etpre,t,dae)
s=size(xx);
s=s(1);
delta=1e-6;
i=0;
J=zeros(s,s);
while i<s
    i=i+1;
    f0=func(h,xx,xpre1,xpre2,alpha1,alpha2,alpha0,betha1,betha2,tpre1,et,etpre,t,dae);
    yy=xx;
    yy(i)=xx(i)+delta;
    f1=func(h,yy,xpre1,xpre2,alpha1,alpha2,alpha0,betha1,betha2,tpre1,et,etpre,t,dae);
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
    dn=max(abs(dn(1:4,1)));
    r=0.9*ETOL/dn;
    r=nthroot(r,3);
end
end


%This function does the quantification
function [x]=quz(x,cq)
%Uniform Quantization
if cq==1 
    range=1;
    q=22;
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
%global Tm omega Vf;
j=sqrt(-1);
TTT=0.05;
Ge=1;
Tee=0.12;
%Parameters
pp=para;
P=pp(1,1);
Q=pp(1,2);
V=pp(1,3);
Xe=pp(1,4);
Re=pp(1,5);
Ra=pp(1,6);
fN=pp(1,7);
Xd=pp(1,8);
Xq=pp(1,9);
Xl=pp(1,10);
Xpd=pp(1,11);
Xpq=pp(1,12);
Tpdo_s=pp(1,13);
Tpqo_s=pp(1,14);
wN=pp(1,15);
TB=pp(1,16);
Le=pp(1,17);
Ldd=pp(1,18);
Ldf=pp(1,19);
Lff=pp(1,20);
Tpdo=pp(1,21);
Rf=pp(1,22);
Lqq=pp(1,23);
Lqq1=pp(1,24);
Lq1q1=pp(1,25);
Tpqo=pp(1,26);
Rq1=pp(1,27);
HH=pp(1,28);
%Extra
if t>1.2
    Xe=Xe*2;
    Re=Re*2;    
end
phasI=(P-j*Q)/V;
phasE=V-(Re+j*Xe)*phasI;
phi=angle(V+Ra*phasI+j*Xq*phasI);
theta0=phi+pi/2;
theta=theta0+wN*t;
ex=real(phasE);
ey=imag(phasE);
%Input
if t>1 && t<1.2
    ex=0;
    ey=0;
end
% if t>1.2
%     Xe=Xe*2;
%     Re=Re*2;    
% end

Le=Xe;
x5=et(1,2); %8.638587815812726e-04
x6=et(1,1); %5013
%Assiginng variables
psif=xx(1,1);
psiq1=xx(2,1);
omega=xx(3,1);
deltas=xx(4,1);
Vf=xx(5,1);
Tm=xx(6,1);
vxm=xx(7,1);
vym=xx(8,1);
freq=xx(9,1);
ifd=xx(10,1);
iq1=xx(11,1);
vd=xx(12,1);
vq=xx(13,1);
vx=xx(14,1);
vy=xx(15,1);
id=xx(16,1);
iq=xx(17,1);
ix=xx(18,1);
iy=xx(19,1);
% ix=0;
% iy=0;
psid=xx(20,1);
psiq=xx(21,1);
Te=psid*iq-psiq*id;
%Differential Equations
f1=wN*(-Rf*ifd+Vf); %psif (INPUT=Vf)
f2=wN*(-Rq1*iq1); %psiq1
f3=(Tm-Te-0.9*omega)/(2*HH); %omega
f4=omega; %deltas
f5=(x5*Ge-Vf)/Tee;
f6=(x6-Tm)/0.3;
% f6=(x6)/0.3;




%Frequency g
f7=(vx-vxm)/TTT;
f8=(vy-vym)/TTT;
%frequency 7
f9=1+omega+(((vy-vym)*vxm)-((vx-vxm)*vym))/(2*pi*fN*TTT*(vxm^2+vym^2))-freq;


%Algebraic equations
f10=Ldd*id+Ldf*ifd-psid; %psid
f11=Lqq*iq+Lqq1*iq1-psiq; %psiq
f12=Lff*ifd+Ldf*id-psif; %psif
f13=Lq1q1*iq1+Lqq1*iq-psiq1; %psiq1
f14=vd+Ra*id+psiq; %vd
f15=vq+Ra*iq-psid; %vq
f16=vx-ex-Re*ix+Xe*iy; %vx
f17=vy-ey-Re*iy-Xe*ix; %vy
f18=vd-vx*cos(theta0)-vy*sin(theta0); %vd
f19=vq-vx*sin(theta0)+vy*cos(theta0); %vq
f20=id-ix*cos(theta0)-iy*sin(theta0); %id
f21=iq-ix*sin(theta0)+iy*cos(theta0); %iq
%sum
f=[f1 f2 f3 f4 f5 f6 f7 f8 f9 f10 f11 f12 f13 f14 f15 f16 f17 f18 f19 f20 f21]';
end

%initial values
function [ssi]=ini
%global Tm omega Vf;
j=sqrt(-1);
pp=para;
P=pp(1,1);
Q=pp(1,2);
V=pp(1,3);
Xe=pp(1,4);
Re=pp(1,5);
Ra=pp(1,6);
fN=pp(1,7);
Xd=pp(1,8);
Xq=pp(1,9);
Xl=pp(1,10);
Xpd=pp(1,11);
Xpq=pp(1,12);
Tpdo_s=pp(1,13);
Tpqo_s=pp(1,14);
wN=pp(1,15);
TB=pp(1,16);
Le=pp(1,17);
Ldd=pp(1,18);
Ldf=pp(1,19);
Lff=pp(1,20);
Tpdo=pp(1,21);
Rf=pp(1,22);
Lqq=pp(1,23);
Lqq1=pp(1,24);
Lq1q1=pp(1,25);
Tpqo=pp(1,26);
Rq1=pp(1,27);
HH=pp(1,28);
%Network
It=(sqrt(P^2+Q^2))/V; %Terminal Current
phi=acos(P/(V*It)); %Power Factor angle
phiR=phi;
deltai=atan((Xq*It*cos(phiR)-Ra*It*sin(phiR))/(V+Ra*It*cos(phiR)+Xq*It*sin(phiR))); %Rotor Angle
vx=V;                                
vy=0;
phasI=(P-j*Q)/V;
ix=real(phasI);
iy=imag(phasI);
phasE=V-(Re+j*Xe)*phasI;
ex=real(phasE);
ey=imag(phasE);
%Machine
phi=angle(V+Ra*phasI+j*Xq*phasI);     
theta0=phi+pi/2;
vd=V*cos(theta0);
vq=V*sin(theta0);
id=abs(phasI)*cos(theta0-angle(phasI));
iq=abs(phasI)*sin(theta0-angle(phasI));
psid=vq+Ra*iq;
psiq=Lqq*iq;
ifd=(psid-Ldd*id)/Ldf;
psif=Lff*ifd+Ldf*id;
Vf=Rf*ifd;
%Vf=8.638587815812726e-04;
Tm=psid*iq-psiq*id;
Tm2=P+Ra*It^2;
iq1=0;
psiq1=Lqq1*iq;
omega=0;
deltas=deltai;
%Output
ssi=[psif psiq1 omega deltas Vf Tm vx vy 1 ifd iq1 vd vq vx vy id iq ix iy psid psiq]';
end

function pp=para
%%%DATA%%%
%Network
P=0.5;
Q=0.1;
V=1;
Xe=0.2;
Re=0.01;
%Machine
Ra=0.005;
fN=50;
Xd=2.4;
Xq=2.4;
Xl=0.2;
Xpd=0.40;
Xpq=0.25;
Tpdo_s=7;
Tpqo_s=0.3;

%%%Parameters%%%
wN=2*pi*fN ;
TB=1/(2*pi*fN) ;
Le=Xe ;
Ldd=Xd ;
Ldf=Xd-Xl ;
Lff=Ldf^2/(Xd-Xpd) ;
Tpdo=Tpdo_s/TB ;
Rf=Lff/Tpdo ;
Lqq=Xq ;
Lqq1=Xq-Xl ;
Lq1q1=Lqq1^2/(Xq-Xpq) ;
Tpqo=Tpqo_s/TB;
Rq1=Lq1q1/Tpqo;
HH=3.5;
pp=[P Q V Xe Re Ra fN Xd Xq Xl Xpd Xpq Tpdo_s Tpqo_s wN TB Le Ldd Ldf Lff Tpdo Rf Lqq Lqq1 Lq1q1 Tpqo Rq1 HH];
end