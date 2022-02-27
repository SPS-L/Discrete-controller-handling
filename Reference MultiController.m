%% %% This is for RTE example variable step landing on the events using zero crossing (True Solution) for multiple controllers
%%%%% anti wind up added
%%%%% Multi controller

clear all
clc

cp=5;
cc=3;
cq=1;

t=0;
h=0.05;
holdd=h;
hmin=0.05;
hmax=1;

% example #5
tsim=65;
x=0;
y=0;
xx=[x;y];
cn=3;
T=[0.1 0.12 0.15];
n=1;
G=[0.05 0.06 0.07];
u=1;

s=size(xx);
s=s(1);
xpre1=xx;
xpre2=xx;
times(1)=t;
xhist(1:s,1)=xx;

x2=0;
iterationsc=0;
tsc=0;
mtsc=0;
etpre2j=zeros(1,cn);
etprej=zeros(1,cn);
etj=zeros(1,cn);
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

    
    %Discrete part
    if t<min(T)
        etj=zeros(1,cn);
        etprej=zeros(1,cn);
        etpre2j=zeros(1,cn);
    else
        [t,conti,tg,mm]=detert(tpre1,t,tg,T,cn,tsim);
        if mm==0
            etpre2j=etprej;
            etprej=etj;
            etj=etj;
            tg=tg;
        else
            etpre2j=etprej;
            etprej=etj;
            etj=contex(conti,etj,G,T,u,x2,cn);
            etj=quz(etj,cq);
            h=t-tpre1;
            flagg=1;
        end
    end
    et=sum(etj);
    etpre=sum(etprej);
    etpre2=sum(etpre2j);
    times(k)=t;

    [alpha0,alpha1,alpha2,betha1,betha2]=coefs(cp,h,holdd);
    xx0=predictor(xx,xpre1,et,etpre,t,tpre1,h,betha1,betha2);
    [alpha0,alpha1,alpha2,betha1,betha2]=coefs(cc,h,holdd);
    [iterations,xx]=NewtonNL(h,xx0,xpre1,xpre2,alpha1,alpha2,alpha0,betha1,betha2,tpre1,et,etpre,t);
    [dn,r]=lote(cp,cc,xx0,xx,h,holdd);
    iterationsc=iterations+iterationsc;

            
   hpre=h;
    if dn>ETOL
%       h=h/2;
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
%         if r*h>hmin
%             hc=hmin;
%         end
        hc=[1.25*h,hmax];
        h=min(hc);
        tsc=tsc+1;
%       h=r*h;
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

% figure(13)
% plot(times,hhist(3,:))
hold on

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
    minl=find(temp==min(temp));
    sminl=size(minl);
    sminl=sminl(1,2);
    j=0;
    while j<sminl
        j=j+1;
        t=min(temp);
        tg(1,minl(1,j))=tg(1,minl(1,j))+T(1,minl(1,j));
    end
end
end

function etj=contex(conti,etj,G,T,u,x2,cn)
jj=size(conti);
jj=jj(1,2);
j=0;
while j<jj
    j=j+1;
    contj=conti(1,j);
    if contj==1
        etj(1,contj)=etj(1,contj)+(G(1,1)*T(1,1)*(u-x2));
        etj(1,contj)=min(etj(1,contj),0.98); %upper limnit anti wind up
        etj(1,contj)=max(etj(1,contj),-0.98); %lower limnit anti wind up
    elseif contj==2
        etj(1,contj)=etj(1,contj)+(G(1,2)*T(1,2)*(u-x2));
        etj(1,contj)=min(etj(1,contj),0.98); %upper limnit anti wind up
        etj(1,contj)=max(etj(1,contj),-0.98); %lower limnit anti wind up
    elseif contj==3
        etj(1,contj)=etj(1,contj)+(G(1,3)*T(1,3)*(u-x2));
        etj(1,contj)=min(etj(1,contj),0.98); %upper limnit anti wind up
        etj(1,contj)=max(etj(1,contj),-0.98); %lower limnit anti wind up
    end
end
end


function [iterations,xx]=NewtonNL(h,xx,xpre1,xpre2,alpha1,alpha2,alpha0,betha1,betha2,tpre1,et,etpre,t)
NTOL=1e-6; %default=1e-8
iterations=0;
a=1;
    while (iterations<30) && (a>NTOL)
    xxpre=xx;
    iterations=iterations+1;
    [J]=Jacob(h,xx,xpre1,xpre2,alpha1,alpha2,alpha0,betha1,betha2,tpre1,et,etpre,t);
    yy=xx-J\func(h,xx,xpre1,xpre2,alpha1,alpha2,alpha0,betha1,betha2,tpre1,et,etpre,t);
    xx=yy;
%     xx=quz(xx,1);
    a=abs(xx-xxpre);
    a=max(a);
%     a=abs(func(h,xx,xpre1,xpre2,alpha1,alpha2,alpha0,betha1,betha2,tpre1,et,t));
%     a=a(1,1);
    if iterations==30
        disp('Did not converged')
    end
    end
end

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
elseif cq==2 %NonUniform Quantization Mu Law
    %Compressing the signal
    mu=255;
    x=sign(x)*((log(1+mu*abs(x)))/(log(1+mu)));
    %Uniform Quantization
    range=1;
    q=16;
    x=round(x*(2^q))*range/(2^q);
    %Expanding the signal
    x=sign(x)*(((1+mu)^abs(x))-1/mu);
elseif cq==3 %NonUniform Quantization A Law
    %Compressing the signal
    A=87.6;
    if abs(x)>=0 && abs(x)<=(1/A) 
        x=(sign(x)/(1+log(A)))*A*abs(x);
    elseif abs(x)>=(1/A) && abs(x)<=1
        x=(sign(x)/(1+log(A)))*(1+log(A*abs(x)));
    end
    %Uniform Quantization    
    range=1;
    q=16;
    x=round(x*(2^q))*range/(2^q);
    %Expanding the signal 
    if abs(x)>=0 && abs(x)<=(1/(1+log(A))) 
        x=sign(x)*((abs(x)*(1+log(A))/A));
    elseif abs(x)>=(1/(1+log(A))) && abs(x)<=1
        x=sign(x)*((exp(abs(x)*((1+log(A)))-1))/(A+(A*log(A))));
    end
end
end

function [A,B]=basedmatrix
a=-0.17;
b=0.9;
A=[2*a sqrt(a^2+b^2);-sqrt(a^2+b^2) 0];
B=[-sqrt(a^2+b^2) 0]';
end