%% %% This is for RTE example variable step landing on the events using zero crossing(True Solution)

clear all
clc
%5-3 2nd order Adams Bashford Adams Moulton, 5-2 2nd order Adams Bashfrod BDF 
cp=5;
cc=3;
%1 uniform quantization, 2 Nonuniform quantization 
cq=1; %Quantization method

t=0;
h=0.05;
holdd=h;
hmin=0.05;
hmax=1;

% example #5
tsim=65; %Time of simulation
x=0; 
y=0;
xx=[x;y];
T=0.1; %Time of sampling by the controller
n=1;
G=0.07; %gainb of the controller
u=1; %input

s=size(xx);
s=s(1);
xpre1=xx;
xpre2=xx;
times(1)=t;
xhist(1:s,1)=xx;

iterationsc=0;
tsc=0;
mtsc=0;
etpre=0;
et=0;
tg=T;
k=1;
ETOL=3e-4;
ek=0;

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
    
    

    
    %Discrete part
    if t<T
        et=0;
        etpre=0;
        etpre2=0;
    elseif tg>tpre1 && tg<=t+0.00000001
        etpre2=etpre;
        etpre=et;
        et=et+(G*T*(u-x2));
%         et=min(et,0.9); %upper limnit anti wind up
%         et=max(et,-0.9); %lower limnit anti wind up
        et=quz(et,cq);
        t=tg;
        h=t-tpre1;
        tg=tg+T;
        flagg=1;
    else
        etpre2=etpre;
        etpre=et;
        et=et;
        tg=tg;
    end
    k=k+1;
    times(k)=t;

    [alpha0,alpha1,alpha2,betha1,betha2]=coefs(cp,h,holdd);
    xx0=predictor(xx,xpre1,et,etpre,t,tpre1,h,betha1,betha2); %predictor
    [alpha0,alpha1,alpha2,betha1,betha2]=coefs(cc,h,holdd);
    [iterations,xx]=NewtonNL(h,xx0,xpre1,xpre2,alpha1,alpha2,alpha0,betha1,betha2,tpre1,et,etpre,t);
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
            x2=xx(2,1);
            if k==2
                hhist(1,1)=h;
                hhist(2,1)=dn;
                hhist(3,1)=r;
                hhist(4,1)=et;
            end
            hhist(1,k)=h;
            hhist(2,k)=dn;
            hhist(3,k)=r;
            hhist(4,k)=et;
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
                    tg=tg-T;
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
            hhist(4,1)=et;
        end
        hhist(1,k)=h;
        hhist(2,k)=dn;
        hhist(3,k)=r;
        hhist(4,k)=et;
        holdd=hpre;
    end


end

figure(1)
plot(times,xhist(2,:))
hold on
xlabel('Time','FontSize', 24)
ylabel('Solution','FontSize', 24)
set(gca,'FontSize',18)

% figure(2)
% plot(xhist(1,:),xhist(3,:))

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



function [iterations,xx]=NewtonNL(h,xx,xpre1,xpre2,alpha1,alpha2,alpha0,betha1,betha2,tpre1,et,etpre,t)
NTOL=1e-4;
iterations=0;
a=1;
    while (iterations<30) && (a>NTOL)
    xxpre=xx;
    iterations=iterations+1;
    [J]=Jacob(h,xx,xpre1,xpre2,alpha1,alpha2,alpha0,betha1,betha2,tpre1,et,etpre,t);
    yy=xx-J\func(h,xx,xpre1,xpre2,alpha1,alpha2,alpha0,betha1,betha2,tpre1,et,etpre,t);
    xx=yy;
    a=abs(xx-xxpre);
    a=max(a);
    if iterations==30
        disp('Did not converged')
    end
    end
end

% This function computes the Jacobian
function [J]=Jacob(h,xx,xpre1,xpre2,alpha1,alpha2,alpha0,betha1,betha2,tpre1,et,etpre,t)
s=size(xx);
s=s(1);
delta=1e-6;
i=0;
J=zeros(s,s);
while i<s
    i=i+1;
    f0=func(h,xx,xpre1,xpre2,alpha1,alpha2,alpha0,betha1,betha2,tpre1,et,etpre,t);
    yy=xx;
    yy(i)=xx(i)+delta;
    f1=func(h,yy,xpre1,xpre2,alpha1,alpha2,alpha0,betha1,betha2,tpre1,et,etpre,t);
    J(:,i)=(f1-f0)/delta;
end
end

% This function determines the predictor formulation
function [xxp]=predictor(xx,xpre1,et,etpre,t,tpre1,h,betha1,betha2)
    xxp=xx+h*(((betha1)*evaleval(xx,t,et))+((betha2)*evaleval(xpre1,tpre1,etpre)));
end


%This function is the chosen method's equation
function [f]=func(h,xx,xpre1,xpre2,alpha1,alpha2,alpha0,betha1,betha2,tpre1,et,etpre,t) 
f=h*(betha1*evaleval(xx,t,et)+betha2*evaleval(xpre1,tpre1,etpre))+(alpha0*xx)+(alpha1*xpre1)+(alpha2*xpre2);
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
    %Adams Bashford 1
    alpha1=0;
    alpha2=0;
    alpha0=0;
    betha1=1;
    betha2=0;
    
elseif c==5
    %Adams Bashford 2
    alpha1=0;
    alpha2=0;
    alpha0=0;
    betha1=1+(h/(holdd*2));
    betha2=-h/(holdd*2);
end
end



%This function is supposed to compute the local truncation error for
%predictor corrector methods (Milne's estimate)
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

function [f]=evaleval(xx,t,et)
% example #1
% sigma=10;
% r=28;
% b=8/3;
% f(1)=sigma*(xx(2)-xx(1));
% f(2)=(r*xx(1))-xx(2)-(xx(1)*xx(3));
% f(3)=(xx(1)*xx(2))-(b*xx(3));
% f=[f(1);f(2);f(3)];
% example #2
%f=(5/xx)-(5*t*(xx^2))-(1/xx^2);
% example #3
%f=-100*(xx-sin(t));
%example %4 (Opne-loop RTE Example)
% u=1;
% [A,B]=basedmatrix;
% f=A*xx+B*u;
%example %5 (close-loop RTE Example)
[A,B]=basedmatrix;
f=A*xx+B*et;
end

%This function does the quantification
function [x]=quz(x,cq)
if cq==1 %Uniform Quantization
    range=1;
    q=16;
    x=round(x*(2^q))*range/(2^q);
elseif cq==4 %NonUniform Quantization A Law
    xp=x;
    if xp>0.8
        x=2*x;
    else
        x=(2/3)*x+(4/3);
    end
    %Uniform Quantization    
    range=1;
    q=12;
    x=round(x*(2^q))*range/(2^q);
    %Expanding the signal
    if xp>0.8
%         x=((3/2)*x)-2;
        x=x/2;
    else
%         x=x/2;
        x=((3/2)*x)-2;
    end
end
end

function [A,B]=basedmatrix
a=-0.11;
b=0.9;
A=[2*a sqrt(a^2+b^2);-sqrt(a^2+b^2) 0];
B=[-sqrt(a^2+b^2) 0]';
end