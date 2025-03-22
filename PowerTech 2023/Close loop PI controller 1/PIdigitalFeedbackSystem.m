%% %% This file solves the anti-windup PI controller using SRM (Digital)
clear all
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
hmax=0.001; %The maximum acceptable time-step size
global T;
T=0.2; %Sampling time of controller
tsim=6.5; %Time of the simulation
%Initial state variables values
ki=3;
kp=1;
kf=0.1;
wmax=2.2;
wmin=-2.2;

dae=[3 0];
x=0.5;
u=0;
y=0.5;
w=0.5;
xx=[u;0;0];
s=size(xx);
s=s(1);
xpre1=xx;
xpre2=xx;
times(1)=t;
xhist(1:s,1)=xx; %Allocating a space to save the results


iterationsc=0; %Allocating a space to save the number of Newoton iterations
tsc=0; %Number of time steps
mtsc=0; %Number of accepted time steps with the manimim value
etpre=0.5;
et=0.5;
et1pre=0.5;
et1=0.5;
tg=T; %Time of the events
k=1;
ETOL=3e-4; %The tolerance to compare with the error estimate
hhist(1,1)=h;
hhist(2,1)=et;
hhist(3,1)=w;
save=0;
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
    
    %Setting new and old controllers's signal
    if t<T
        et=0.5;
        etpre=0.5;
        etpre2=0.5;
        et1=0.5;
        et1pre=0.5;
        et1pre2=0.5;
    elseif tg>tpre1 && tg<=t+0.00000001
        etpre2=etpre;
        etpre=et;
        et1pre2=et1pre;
        et1pre=et1;
        if et>wmin && et<wmax
            x=ki*(xx(1,1)-kf*xx(3,1))*T+x;
        else
            x=x;
        end
        et=kp*(xx(1,1)-kf*xx(3,1))+x;
        et1=min(et,wmax); %upper limnit anti wind up
        et1=max(et1,wmin); %lower limnit anti wind up
        w=et1;
%         t=tg;
%         h=t-tpre1; %Zero-crossing
        tg=tg+T;
    else
        etpre2=etpre;
        etpre=et;
        et=et;
        et1pre2=et1pre;
        et1pre=et;
        et1=et1;
        tg=tg;
    end
    
    k=k+1;
    %Summation of the controllers' signal
    times(k)=t;
    [alpha0,alpha1,alpha2,betha1,betha2]=coefs(cp,h,holdd); %calculating the coefficients of the predictor
    xx0=predictor(xx,xpre1,t,tpre1,h,betha1,betha2,dae,et1,et1pre); %predictor
    [alpha0,alpha1,alpha2,betha1,betha2]=coefs(cc,h,holdd); %calculating the coefficients of the corrector
    [iterations,xx]=NewtonNL(h,xx0,xpre1,xpre2,alpha1,alpha2,alpha0,betha1,betha2,tpre1,t,dae,et1,et1pre); %The main function (Newton solver)
    [dn,r]=lote(cp,cc,xx0,xx,h,holdd); %Error estimate
    iterationsc=iterations+iterationsc; %Counting the number of Newton iterations
   %Calculating the new time step and saving the results for the current time step         
    hpre=h;
    tsc=tsc+1;
    xhist(1:s,k)=xx;
    hhist(1,k)=h;
    hhist(2,k)=et;
    hhist(3,k)=w;
    hhist(4,k)=x;
    hhist(5,k)=iterations;
    holdd=hpre;
end
toc
%Plotting the results
figure(1)
plot(times,hhist(2,:))
hold on
% plot(times,hhist(3,:))
% plot(times,xhist(1,:))
% plot(times,xhist(2,:))
xlabel('Time (s)','FontSize', 24)
ylabel('Controller Output','FontSize', 24)
set(gca,'FontSize',18)
legend('Analog','Digital (T=0.001)','Digital (T=0.002)','Digital (T=0.02)','Digital (T=0.2)','Location','northeast')
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
    hgexport(gcf,'C:\Users\mjafa\Desktop\Coding tutorial\MATLAB\PI controller paper\PIAnalogFeedback1.pdf',figure_property);
end

figure(2)
plot(times,xhist(3,:))
hold on
xlabel('Time (s)','FontSize', 24)
ylabel('System Output (z_{2})','FontSize', 24)
set(gca,'FontSize',18)
legend('Analog','Digital (T=0.001)','Digital (T=0.002)','Digital (T=0.02)','Digital (T=0.2)','Location','northwest')
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
% figure_property.FixedfontSize= '12';
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
    hgexport(gcf,'C:\Users\mjafa\Desktop\Coding tutorial\MATLAB\PI controller paper\PIAnalogFeedback2.pdf',figure_property);
end

figure(5)
bar(times,hhist(5,:))
xlim([4.8 5.6])
hold on
xlabel('Time','FontSize', 24)
ylabel('Number of iterations','FontSize', 24)
set(gca,'FontSize',18)

% figure(3)
% plot(times,xhist(1,:)-xhist(3,:))
% hold on
% xlabel('Time','FontSize', 24)
% ylabel('Output','FontSize', 24)
% set(gca,'FontSize',18)

figure(11)
plot(times,hhist(1,:))
hold on
xlabel('Time','FontSize', 24)
ylabel('Step Size','FontSize', 24)
set(gca,'FontSize',18)

figure(12)
plot(times,hhist(2,:))
xlabel('Time','FontSize', 24)
ylabel('Error Estimate','FontSize', 24)
set(gca,'FontSize',18)
hold on

function [iterations,xx]=NewtonNL(h,xx,xpre1,xpre2,alpha1,alpha2,alpha0,betha1,betha2,tpre1,t,dae,et1,et1pre)
NTOL=1e-4; %Tolerance for the Newton solver
iterations=0;
a=1;
    while (iterations<10) && (a>NTOL)
    xxpre=xx;
    iterations=iterations+1;
    [J]=Jacob(h,xx,xpre1,xpre2,alpha1,alpha2,alpha0,betha1,betha2,tpre1,t,dae,et1,et1pre); %calculating the Jacobian
    yy=xx-J\func(h,xx,xpre1,xpre2,alpha1,alpha2,alpha0,betha1,betha2,tpre1,t,dae,et1,et1pre); %Newton
    xx=yy;
    a=abs(xx-xxpre);
    a=max(a);
    if iterations==10
        disp('Did not converged')
        t=t
    end
    end
end

% This function computes the Jacobian
function [J]=Jacob(h,xx,xpre1,xpre2,alpha1,alpha2,alpha0,betha1,betha2,tpre1,t,dae,et1,et1pre)
s=size(xx);
s=s(1);
delta=1e-6;
i=0;
J=zeros(s,s);
while i<s
    i=i+1;
    f0=func(h,xx,xpre1,xpre2,alpha1,alpha2,alpha0,betha1,betha2,tpre1,t,dae,et1,et1pre);
    yy=xx;
    yy(i)=xx(i)+delta;
    f1=func(h,yy,xpre1,xpre2,alpha1,alpha2,alpha0,betha1,betha2,tpre1,t,dae,et1,et1pre);
    J(:,i)=(f1-f0)/delta;
end
end

% This function determines the predictor formulation
function [xxp]=predictor(xx,xpre1,t,tpre1,h,betha1,betha2,dae,et1,et1pre)
    xxp1=xx+h*(((betha1)*evaleval(t,xx,et1))+((betha2)*evaleval(tpre1,xpre1,et1pre)));
    xxp2=xx+evaleval(t,xx,et1);

if dae(1,1)==0
    xxp=xxp2;
elseif dae(1,2)==0
    xxp=xxp1;
else
    xxp=[xxp1(1:dae(1,1),1);xxp2(dae(1,1)+1:end,1)];
end
end


%This function determines the corrector formulation
function [f]=func(h,xx,xpre1,xpre2,alpha1,alpha2,alpha0,betha1,betha2,tpre1,t,dae,et1,et1pre) 
    f1=h*(betha1*evaleval(t,xx,et1)+betha2*evaleval(tpre1,xpre1,et1pre))+(alpha0*xx)+(alpha1*xpre1)+(alpha2*xpre2);
    f2=evaleval(t,xx,et1);

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
    q=16;
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
    q=12;
    x=round(x*(2^q))*range/(2^q);
    if xp>0.8
        x=x/2;
    else
        x=((3/2)*x)-2;
    end
end
end

function [f]=evaleval(t,xx,et1)
%Variable assignment
xxx=[xx(2,1);xx(3,1)];
%Data
%y calculation
if t<3
   f1=1; 
else
   f1=-1; 
end
[A,B]=basedmatrix;
f2=A*xxx+B*et1;
f=[f1;f2];
end

%System's eigenvalues
function [A,B]=basedmatrix
a=-0.1;
b=0.9;
A=[2*a sqrt(a^2+b^2);-sqrt(a^2+b^2) 0];
B=[-sqrt(a^2+b^2) 0]';
end