%% %% This file solves the RTE example (single PI controller) using the interpolation based method (IBM)
%%%% Please set the parameters based on the guidance provided in the Readme.md file 

clear all
clc
% This parameter chooses the predictor method. It can take values 4 and 5. 4 is for the first-order Adams-Bashford, and 5 is for the second-order Adams-Bashford. 
cp=5; %Predictor Choice
% This parameter chooses the corrector method. It can take values 2 and 3. 2 is for the second-order BDF, and 3 is for the second-order Adams-Moulton. 
cc=3; %Corrector Choice
%This parameters chooses the quantization method
%1 uniform quantization, 2 Nonuniform quantization
cq=1; %Quantization method

t=0; %Starting point of the integration
h=0.05; %The first time-step size
holdd=h;
hmin=0.05; %The minimum acceptable time-step size
hmax=1; %The maximum acceptable time-step size

tsim=65; %Time of simulation
%Initial state variables values
x=0;
y=0;
xx=[x;y];
T=0.1; %Time of sampling by the controller
G=0.07; %Gain of the controller
u=1; %Input

s=size(xx);
s=s(1);
xpre1=xx;
xpre2=xx;
times(1)=t;
xhist(1:s,1)=xx; %Allocating a space to save the results


eeghist=0;
feeghist=0;
iterationsc=0; %Allocating a space to save the number of Newoton iterations
tsc=0; %Number of time steps
mtsc=0; %number of accepted time steps with the manimim value
egpre=0;
eg=0;
tg=T; %Time of the events
k=1;
ETOL=3e-4; %The tolerance to compare with the error estimate

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
    if k>1
        xpre2=xhist(1:s,k-1);
    else
        xpre2=xpre1;
    end
    xpre1=xhist(1:s,k);
       
    %Setting old controllers's signal
    if k>1
        egpre=eg;
    else
        egpre=0;
    end
    k=k+1;
    times(k)=t;
    

    [alpha0,alpha1,alpha2,betha1,betha2]=coefs(cp,h,holdd); %calculating the coefficients of the predictor   
    [xx0g,ee0g,eg,tg,sg]=predictorg(xx,xpre1,t,tpre1,tpre2,holdd,T,eg,egpre,tg,cp,G,u,cq); %predictor of the intermediate points (explicit predictor)
    xx0=predictor(xx,xpre1,eg,egpre,t,tpre1,tpre2,h,betha1,betha2); %predictor
    [alpha0,alpha1,alpha2,betha1,betha2]=coefs(cc,h,holdd); %calculating the coefficients of the corrector
    [iterations,xx,eg,tg,feeg,eeg]=NewtonNL(h,holdd,xx0,xpre1,xpre2,alpha1,alpha2,alpha0,betha1,betha2,tpre1,t,sg,eg,T,tg,G,u,ee0g,egpre,cq); %The main function (Newton solver)
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
            if k==2
                hhist(1,1)=h;
                hhist(2,1)=dn;
                hhist(3,1)=r;
                hhist(4,1)=eg;
                eeghist=[eeghist eeg'];
            end
            hhist(1,k)=h;
            hhist(2,k)=dn;
            hhist(3,k)=r;
            hhist(4,k)=eg;
            holdd=hpre;
            eeghist=[eeghist eeg'];
        else

                t=t-hpre;
                tpre1=tpre2;
                if k>2
                    tpre2=tpre3;
                else
                    tpre2=0;
                end
                tg=tg-sg*T;
                eg=egpre;
                k=k-1;

        end
    else
        hc=[1.25*h,hmax];
        h=min(hc);
        tsc=tsc+1;
        xhist(1:s,k)=xx;
        if k==2
            hhist(1,1)=h;
            hhist(2,1)=dn;
            hhist(3,1)=r;
            hhist(4,1)=eg;
            eeghist=[eeghist eeg'];
        end
        hhist(1,k)=h;
        hhist(2,k)=dn;
        hhist(3,k)=r;
        hhist(4,k)=eg;
        holdd=hpre;
        eeghist=[eeghist eeg'];
    end
end

% sgsize=size(eeghist);
% sgsize=sgsize(1,2);
% timessg=0.1:0.1:sgsize*0.1;


%Plotting the results
figure(1)
plot(times,xhist(2,:),'+')
hold on
%,'+'


% figure(2)
% plot(xhist(1,:),xhist(3,:))

figure(11)
plot(times,hhist(1,:))
hold on

figure(12)
plot(times,hhist(2,:))
hold on
% figure(13)
% plot(times,hhist(4,:))
% figure(13)
% plot(times,hhist(3,:))
% figure(14)
% plot(timessg,feeghist)
% hold on
% figure(15)
% plot(timessg,eeghist)
% hold on
% figure(16)
% plot(times,hhist(4,:))
hold on

function [iterations,xx,eg,tg,feeg,eeg]=NewtonNL(h,holdd,xx,xpre1,xpre2,alpha1,alpha2,alpha0,betha1,betha2,tpre1,t,sg,eg,T,tg,G,u,eeg,egpre,cq)
NTOL=1e-4; %Tolerance for the Newton solver
iterations=0;
a=1;
    while (iterations<30) && (a>NTOL)
        xxpre=xx;
        iterations=iterations+1;
        s=size(xx);
        s=s(1);
        if sg==0
            xxx=xx;
        else
            xxx=[xx;eeg]; %extending the state variables for controller intermediate signals
        end
        [J]=Jacob(h,xx,xpre1,xpre2,alpha1,alpha2,alpha0,betha1,betha2,tpre1,eg,egpre,t); %calculating the Jacobian
        JJ=Jadd(J,sg); %Extending Jacobian for intermediate points
        [feeg,tg]=funcg(xx,xpre1,xpre2,eg,t,tpre1,T,eeg,egpre,tg,G,u,h,holdd,cq); %Corrector of intermediate points (Implicit corrector)
        fy=func(h,xx,xpre1,xpre2,alpha1,alpha2,alpha0,betha1,betha2,tpre1,eg,egpre,t); %Corrector
        if sg==0
            f=fy;
        else
            f=[fy;feeg']; %extending mismatch for intermediate points
        end
        xxx=xxx-JJ\f; %Newton
        %Updating the state variables
        xx=xxx(1:s,1);
        if sg~=0
            eeg=xxx(s+1:end,1);
            eeg=quz(eeg,cq);
        end
        if sg==0
            eg=eg;
        else
            eg=eeg(end,1);
        end
        a=abs(xx-xxpre);
        a=max(a);
        if iterations==30
            disp('Did not converged')
        end
    end
    tg=tg+sg*T; %Updating the  events time
end

% This function computes the Jacobian
function [J]=Jacob(h,xx,xpre1,xpre2,alpha1,alpha2,alpha0,betha1,betha2,tpre1,eg,egpre,t)
s=size(xx);
s=s(1);
delta=1e-6;
i=0;
J=zeros(s,s);
while i<s
    i=i+1;
    f0=func(h,xx,xpre1,xpre2,alpha1,alpha2,alpha0,betha1,betha2,tpre1,eg,egpre,t);
    yy=xx;
    yy(i)=xx(i)+delta;
    f1=func(h,yy,xpre1,xpre2,alpha1,alpha2,alpha0,betha1,betha2,tpre1,eg,egpre,t);
    J(:,i)=(f1-f0)/delta;
end
end

% This function extends the Jacobian for intermediate points
function JJ=Jadd(J,sg)
sj=size(J);
sj=sj(1);
JJ=[J zeros(sj(1),sg);zeros(sg,sj(1)) eye(sg)];
end

% This function determines the predictor formulation
function [xxp]=predictor(xx,xpre1,eg,egpre,t,tpre1,tpre2,h,betha1,betha2)
    xxp=xx+h*(((betha1)*evaleval(xx,tpre1,eg))+((betha2)*evaleval(xpre1,tpre2,egpre)));
end


%This function determines the corrector formulation
function [f]=func(h,xx,xpre1,xpre2,alpha1,alpha2,alpha0,betha1,betha2,tpre1,eg,egpre,t) 
f=h*(betha1*evaleval(xx,t,eg)+betha2*evaleval(xpre1,tpre1,egpre))+(alpha0*xx)+(alpha1*xpre1)+(alpha2*xpre2);
end

%Implicit interpolator
function [eeeg,tg]=funcg(xx,xpre1,xpre2,eg,t,tpre1,T,eeg,egpre,tg,G,u,h,holdd,cq)
tgpre=tg;
if tg>tpre1 && tg<=t+0.000001
    m=0;
    while (tg<=t+0.000001)       
        m=m+1;
        hg=tg-tpre1; %Calculating the time step for the intermediate points
        
%%%chossing the interpolator (only the chosen one should be decomented)
%         xg=hg*(evaleval(xx,t,eg))+(xpre1); %First order AM implicit interpolator
        xg=(hg*(evaleval(xx,t,eg))+((1+hg/holdd)*xpre1)+((-(hg^2)/((hg+holdd)*holdd))*xpre2))/((2*hg+holdd)/(hg+holdd)); %Second order BDF implicit interpolator
%         xg=xpre1+hg*evaleval(xpre1,tpre1,egpre)+((hg^2/h)*(xx-xpre1-(h*evaleval(xpre1,tpre1,egpre)))); %Second order AM implicit interpolator

        xg=quz(xg,cq);
        if m==1
            aw=egpre+G*T*(u-xg(2,1)); %Controller's equation
%             aw=min(aw,0.9); %upper limit anti wind up
%             aw=max(aw,-0.9); %Downer limit anti wind up
            f2=eeg(m,1)-aw;
            eeeg(1,1)=f2;
        else
            aw=+eeg(m-1,1)+G*T*(u-xg(2,1)); %Controller's equation
%             aw=min(aw,0.9); %upper limnit anti wind up
%             aw=max(aw,-0.9); %Downer limit anti wind up
            f2=eeg(m,1)-aw;            
            eeeg(1,m)=f2;
        end
    tg=tg+T;
    end
    

else   
    tg=tg;
    eeg=eeg;
    eeeg=0;
end
tg=tgpre; %reset the counter of the intermediate points to be counted again in the next iterations
end

%Explicit interpolator
function [xx0g,ee0g,eg,tg,sg]=predictorg(xx,xpre1,t,tpre1,tpre2,holdd,T,eg,egpre,tg,cp,G,u,cq)
tgpre=tg;
if (tg > tpre1) && (tg <= t+0.000001)
    m=0;
    while tg <= t+0.000001
        m=m+1;
        hg=tg-tpre1;
        [alpha0,alpha1,alpha2,betha1,betha2]=coefs(cp,hg,holdd);

%%%chossing the interpolator (only the chosen one should be decomented)
%         xg=xpre1+hg*((evaleval(xpre1,tpre1,egpre))); %first order Adams Bashford explicit interpolator
        xg=xx+hg*(((betha1)*evaleval(xx,tpre1,eg))+((betha2)*evaleval(xpre1,tpre2,egpre))); %second order Adams Bashford explicit interpolator
        xx0g(m,:)=xg;
        eg=eg+G*T*(u-xg(2,1));
%         eg=min(eg,0.9); %upper limnit anti wind up
%         eg=max(eg,-0.9); %Downer limit anti wind up
        ee0g(m,:)=eg;
        tg=tg+T;
    end
    xx0g1=quz(xx0g(:,2),cq);
    xx0g(:,2)=xx0g1;
    ee0g=quz(ee0g,cq);
    xx0g=xx0g(:,2);
    sg=size(xx0g);
    sg=sg(1,1);
else
    xx0g=0;
    tg=tg;
    eg=eg;
    ee0g=0;
    sg=0;
end
tg=tgpre; %reset the counter of the intermediate points to be counted again in the next iterations
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
    dn=max(abs(dn(:)));
    r=0.9*ETOL/(dn*h);
    r=nthroot(r,3);    
% elseif cc==3 && cp==4
%     dn=(yy-yp)/3*(1+holdd/h);
%     r=0.9*ETOL/(dn*h);
%     r=nthroot(r,3);                            
elseif cc==3 && cp==5
    dn=(yy-yp)/(3*(1+holdd/h));
    dn=max(abs(dn(:)));
    r=0.9*ETOL/dn;
    r=nthroot(r,3);
end
end

%This function defines the system's formulation
function [f]=evaleval(xx,t,eg)
%example1 RTE PI Controller Close-loop
[A,B]=basedmatrix;
f=A*xx+B*eg;
%example1 RTE PI Controller open-loop
% [A,B]=basedmatrix;
% f=A*xx+B;
%example2 PI controller
% ki=3.1;
% kp=2;
% f(1)=ki*et;
% f(2)=-xx(1)+xx(2)+xx(3);
% f(3)=-xx(2)+kp*et;
% f=[f(1);f(2);f(3)];
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
%System's eigenvalues
function [A,B]=basedmatrix
a=-0.01;
b=0.9;
A=[2*a sqrt(a^2+b^2);-sqrt(a^2+b^2) 0];
B=[-sqrt(a^2+b^2) 0]';
end