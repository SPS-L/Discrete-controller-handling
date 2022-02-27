%% %% IBM method for RTE example
%%%%% anti wind up added

clear all
clc
%5-3 2nd order Adams Bashford Adams Moulton, 5-2 2nd order Adams Bashfrod BDF 
cp=5; %Predictor Choice
cc=3; %Corrector Choice

t=0;
h=0.05;
holdd=h;
hmin=0.05;
hmax=1;

% example #1
tsim=65; %Time of simulation
x=0;
y=0;
xx=[x;y];
T=0.1; %Time of sampling by the controller
n=1;
G=0.07; %gainb of the controller
u=1; %input

%example2 PI controller
% tsim=65;
% y=0;
% x1=0;
% x2=0;
% xx=[y;x1;x2];
% T=0.1;
% n=1;
% G=0.09;
% u=1;

s=size(xx);
s=s(1);
xpre1=xx;
xpre2=xx;
times(1)=t;
xhist(1:s,1)=xx;

timessg=0;
eeghist=0;
feeghist=0;
iterationsc=0;
tsc=0;
mtsc=0;
egpre=0;
eg=0;
tg=T;
k=1;
ETOL=3e-4;

while t<tsim
    %old times
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
    
    %old solutions
    if k>1
        xpre2=xhist(1:s,k-1);
    else
        xpre2=xpre1;
    end
    xpre1=xhist(1:s,k);
    
    
    
    
    %Last digital signal
    if k>1
        egpre=eg;
    else
        egpre=0;
    end
    k=k+1;
    times(k)=t;
    

    [alpha0,alpha1,alpha2,betha1,betha2]=coefs(cp,h,holdd);   
    [xx0g,ee0g,eg,tg,sg]=predictorg(xx,xpre1,t,tpre1,tpre2,holdd,T,eg,egpre,tg,cp,G,u); %predictor of intermediate points
    xx0=predictor(xx,xpre1,eg,egpre,t,tpre1,tpre2,h,betha1,betha2); %predictor
    [alpha0,alpha1,alpha2,betha1,betha2]=coefs(cc,h,holdd);
    [iterations,xx,eg,tg,feeg,eeg]=NewtonNL(h,holdd,xx0,xpre1,xpre2,alpha1,alpha2,alpha0,betha1,betha2,tpre1,t,sg,eg,T,tg,G,u,ee0g,egpre);
    [dn,r]=lote(cp,cc,xx0,xx,h,holdd); %error estimate
    iterationsc=iterations+iterationsc;
    

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

sgsize=size(eeghist);
sgsize=sgsize(1,2);
timessg=0.1:0.1:sgsize*0.1;



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
figure(16)
plot(times,hhist(4,:))
hold on

function [iterations,xx,eg,tg,feeg,eeg]=NewtonNL(h,holdd,xx,xpre1,xpre2,alpha1,alpha2,alpha0,betha1,betha2,tpre1,t,sg,eg,T,tg,G,u,eeg,egpre)
NTOL=1e-4; 
iterations=0;
a=1;
    while (iterations<5) && (a>NTOL)
        xxpre=xx;
        iterations=iterations+1;
        s=size(xx);
        s=s(1);
        if sg==0
            xxx=xx;
        else
            xxx=[xx;eeg]; %extending the state variables for controller intermediate signals
        end
        [J]=Jacob(h,xx,xpre1,xpre2,alpha1,alpha2,alpha0,betha1,betha2,tpre1,eg,egpre,t); %Jacobian
        JJ=Jadd(J,sg); %Extending Jacobian for intermediate points
        [feeg,tg]=funcg(xx,xpre1,xpre2,eg,t,tpre1,T,eeg,egpre,tg,G,u,h,holdd); %corrector of intermediate points
        fy=func(h,xx,xpre1,xpre2,alpha1,alpha2,alpha0,betha1,betha2,tpre1,eg,egpre,t); %corrector
        if sg==0
            f=fy;
        else
            f=[fy;feeg']; %extending mismatch for intermediate points
        end
        xxx=xxx-JJ\f;
        xx=xxx(1:s,1);
        if sg~=0
            eeg=xxx(s+1:end,1);
            eeg=quz(eeg);
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
    tg=tg+sg*T;
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


%This function is the chosen method's equation
function [f]=func(h,xx,xpre1,xpre2,alpha1,alpha2,alpha0,betha1,betha2,tpre1,eg,egpre,t) 
f=h*(betha1*evaleval(xx,t,eg)+betha2*evaleval(xpre1,tpre1,egpre))+(alpha0*xx)+(alpha1*xpre1)+(alpha2*xpre2);
end

%This is the equations of the controller for the intermediate points to be
%added to the newton solver
function [eeeg,tg]=funcg(xx,xpre1,xpre2,eg,t,tpre1,T,eeg,egpre,tg,G,u,h,holdd)
tgpre=tg;
if tg>tpre1 && tg<=t+0.000001
    m=0;
    while (tg<=t+0.000001)       
        m=m+1;
        hg=tg-tpre1;
%         xg=hg*(evaleval(xx,t,eg))+(xpre1); %First order AM implicit interpolator
        xg=(hg*(evaleval(xx,t,eg))+((1+hg/holdd)*xpre1)+((-(hg^2)/((hg+holdd)*holdd))*xpre2))/((2*hg+holdd)/(hg+holdd)); %Second order BDF implicit interpolator
%         xg=xpre1+hg*evaleval(xpre1,tpre1,egpre)+((hg^2/h)*(xx-xpre1-(h*evaleval(xpre1,tpre1,egpre)))); %Second order AM implicit interpolator
        xg=quz(xg);
        if m==1
            aw=egpre+G*T*(u-xg(2,1));
%             aw=min(aw,0.9); %upper limnit anti wind up
%             aw=max(aw,-0.9); %Downer limit anti wind up
            f2=eeg(m,1)-aw;
            eeeg(1,1)=f2;
        else
            aw=+eeg(m-1,1)+G*T*(u-xg(2,1));
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

function [xx0g,ee0g,eg,tg,sg]=predictorg(xx,xpre1,t,tpre1,tpre2,holdd,T,eg,egpre,tg,cp,G,u)
tgpre=tg;
if (tg > tpre1) && (tg <= t+0.000001)
    m=0;
    while tg <= t+0.000001
        m=m+1;
        hg=tg-tpre1;
        [alpha0,alpha1,alpha2,betha1,betha2]=coefs(cp,hg,holdd);
%         xg=xpre1+hg*((evaleval(xpre1,tpre1,egpre))); %first order Adams Bashford explicit interpolator
        xg=xx+hg*(((betha1)*evaleval(xx,tpre1,eg))+((betha2)*evaleval(xpre1,tpre2,egpre))); %second order Adams Bashford explicit interpolator
        xx0g(m,:)=xg;
        eg=eg+G*T*(u-xg(2,1));
%         eg=min(eg,0.9); %upper limnit anti wind up
%         eg=max(eg,-0.9); %Downer limit anti wind up
        ee0g(m,:)=eg;
        tg=tg+T;
    end
    xx0g1=quz(xx0g(:,2));
    xx0g(:,2)=xx0g1;
    ee0g=quz(ee0g);
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

%This function computes the coefficients based on the step sizes and the
%method chosen
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
    alpha1=0;
    alpha2=0;
    alpha0=0;
    betha1=1;
    betha2=0;
    
elseif c==5
    %predictor coefficients
    alpha1=0;
    alpha2=0;
    alpha0=0;
    betha1=1+(h/(holdd*2));
    betha2=-h/(holdd*2);
end
end



%This function compute the local truncation error for predictor corrector methods (Milne's estimate)
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
function [x]=quz(x)
range=1;
q=16;
x=round(x*(2^q))*range/(2^q);

end

function [A,B]=basedmatrix
a=-0.11;
b=0.9;
A=[2*a sqrt(a^2+b^2);-sqrt(a^2+b^2) 0];
B=[-sqrt(a^2+b^2) 0]';
end