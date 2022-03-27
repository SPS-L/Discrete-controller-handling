%% %% This file solves the RTE example (multiple PI controller) using the interpolation based method (IBM)
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
hmin=0.001; %The minimum acceptable time-step size
hmax=1; %The maximum acceptable time-step size

tsim=65; %Time of simulation
%Initial state variables values
x=0;
y=0;
xx=[x;y];
cn=3; %Number of controllers
T=[0.1 0.12 0.15]; %Time of sampling by the controllers
G=[0.05 0.06 0.07]; %Gain of the controllers
u=1; %Input

s=size(xx);
s=s(1);
xpre1=xx;
xpre2=xx;
times(1)=t;
xhist(1:s,1)=xx; %Allocating a space to save the results

timessg=0;
eeghist=0;
feeghist=0;
iterationsc=0; %Allocating a space to save the number of Newoton iterations
tsc=0; %Number of time steps
mtsc=0; %Number of accepted time steps with the manimim value
egpre=0;
eg=zeros(1,cn);
egt=0; %Summation of controllers' signal for the current time step
egpret=0; %Summation of controllers' signal for the previous time step
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
    tpre1=t;
    t=t+h;
    
    %Setting old solutions
    if k>1
        xpre2=xhist(1:s,k-1);
    else
        xpre2=xpre1;
    end
    xpre1=xhist(1:s,k);
    

    
    %Setting old controllers' signal
    if k>2
        egpre2=egpre;
    else
        egpre2=zeros(1,cn);
    end
    if k>1
        egpre=eg;
    else
        egpre=zeros(1,cn);
    end
    if k>1
        egpret=egt;
    else
        egpret=0;
    end
    
    k=k+1;
    times(k)=t;
    tghold=tg;

    [alpha0,alpha1,alpha2,betha1,betha2]=coefs(cp,h,holdd); %calculating the coefficients of the predictor  
    [xx0g,ee0gm,ee0gv,eg,egt,egpret,tg,sg]=predictorg(xx,xpre1,t,tpre1,tpre2,holdd,T,eg,egpre,tg,cp,G,u,cn,cq); %predictor of the intermediate points (explicit predictor)
    xx0=predictor(xx,xpre1,egt,egpret,t,tpre1,tpre2,h,betha1,betha2); %predictor
    [alpha0,alpha1,alpha2,betha1,betha2]=coefs(cc,h,holdd); %calculating the coefficients of the corrector
    [iterations,xx,eg,egt,tg,feeg,eeg]=NewtonNL(h,xx0,xpre1,xpre2,alpha1,alpha2,alpha0,betha1,betha2,tpre1,t,sg,eg,egt,T,tg,G,u,ee0gm,ee0gv,egpre,egpret,cn,cq); %The main function (Newton solver)
    iterationsc=iterations+iterationsc; %Counting the number of Newton iterations
    [dn,r]=lote(cp,cc,xx0,xx,h,holdd); %Error estimate
    
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
                hhist(4,1)=egt;
            end
            hhist(1,k)=h;
            hhist(2,k)=dn;
            hhist(3,k)=r;
            hhist(4,k)=egt;
            holdd=hpre;
        else
            t=t-hpre;
            tpre1=tpre2;
            if k>2
                tpre2=tpre3;
            else
                tpre2=0;
            end
            tg=tghold;
            eg=egpre;
            egpre=egpre2;
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
            hhist(4,1)=egt;
        end
        hhist(1,k)=h;
        hhist(2,k)=dn;
        hhist(3,k)=r;
        hhist(4,k)=egt;
        holdd=hpre;
    end

    

    feeghist=[feeghist feeg];
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

function [iterations,xx,eg,egt,tg,feeg,eegv]=NewtonNL(h,xx,xpre1,xpre2,alpha1,alpha2,alpha0,betha1,betha2,tpre1,t,sg,eg,egt,T,tg,G,u,eegm,eegv,egpre,egpret,cn,cq)
NTOL=1e-4; %Tolerance for the Newton solver
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
            xxx=[xx;eegv]; %extending the state variables for controller intermediate signals
        end
        [J]=Jacob(h,xx,xpre1,xpre2,alpha1,alpha2,alpha0,betha1,betha2,tpre1,egt,egpret,t); %calculating the Jacobian
        JJ=Jadd(J,sg); %Extending the Jacobian for intermediate points
        [feeg,tg,cnhist]=funcg(xx,xpre1,eg,egt,t,tpre1,T,eegm,eegv,egpre,egpret,tg,G,u,h,cn,cq); %Corrector of intermediate points (Implicit corrector)
        fy=func(h,xx,xpre1,xpre2,alpha1,alpha2,alpha0,betha1,betha2,tpre1,egt,egpret,t); %Corrector
        if sg==0
            f=fy;
        else
            f=[fy;feeg']; %extending mismatch for intermediate points
        end
        xxx=xxx-JJ\f; %Newton
        %Updating the state variables
        xx=xxx(1:s,1);
        if sg~=0
            eegv=xxx(s+1:end,1);
            eegv=quz(eegv,cq);
        end
        if sg==0
            egt=egt;
        else
            %Forming new eegm (Saving new values in matrix form)
            cni=0;
            nnzj=1;
            eegm=zeros(sg,cn);
            while cni<cn
                cni=cni+1;
                nnzcni=nnz(cnhist(:,cni));
                eegm(nnzj:nnzj+nnzcni-1,cni)=eegv(nnzj:nnzj+nnzcni-1,1);
                nnzj=nnzj+nnzcni;    
            end
            %Calculating egt
            cni=0;
            while cni<cn
                cni=cni+1;
                cc(1,cni)=nnz(eegm(:,cni));
                if cni==1
                    b(1,cni)=nnz(eegm(:,cni));
                else
                    b(1,cni)=nnz(eegm(:,cni))+b(1,cni-1);
                end
            end
            %calculating eg vector
            cni=0;
            while cni<cn
                cni=cni+1;
                if cc(1,cni)~=0
                    eg(1,cni)=eegm(b(1,cni),cni);
                end
            end
            
            %calculating egt scalar
            egt=sum(eg);
        end

        a=abs(xx-xxpre);
        a=max(a);
        if iterations==30
            disp('Did not converged')
        end
    end
%calculating new set of tg for the next time step
if sg~=0
cni=0;
    while cni<cn
        cni=cni+1;
        sgt=nnz(eegm(:,cni));
        if sgt>0
            tg(1,cni)=tg(1,cni)+(sgt)*T(1,cni);
        end
    end
else
    tg=tg;
end
end

% This function computes the Jacobian
function [J]=Jacob(h,xx,xpre1,xpre2,alpha1,alpha2,alpha0,betha1,betha2,tpre1,egt,egpret,t)
s=size(xx);
s=s(1);
delta=1e-6;
i=0;
J=zeros(s,s);
while i<s
    i=i+1;
    f0=func(h,xx,xpre1,xpre2,alpha1,alpha2,alpha0,betha1,betha2,tpre1,egt,egpret,t);
    yy=xx;
    yy(i)=xx(i)+delta;
    f1=func(h,yy,xpre1,xpre2,alpha1,alpha2,alpha0,betha1,betha2,tpre1,egt,egpret,t);
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
function [xxp]=predictor(xx,xpre1,egt,egpret,t,tpre1,tpre2,h,betha1,betha2)
    xxp=xx+h*(((betha1)*evaleval(xx,tpre1,egt))+((betha2)*evaleval(xpre1,tpre2,egpret)));
end


%This function determines the corrector formulation
function [f]=func(h,xx,xpre1,xpre2,alpha1,alpha2,alpha0,betha1,betha2,tpre1,egt,egpret,t) 
f=h*(betha1*evaleval(xx,t,egt)+betha2*evaleval(xpre1,tpre1,egpret))+(alpha0*xx)+(alpha1*xpre1)+(alpha2*xpre2);
end

%Implicit interpolator
function [eeeg,tg,cnhist]=funcg(xx,xpre1,eg,egt,t,tpre1,T,eegm,eegv,egpre,egpret,tg,G,u,h,cn,cq)
tgpre=tg;
cni=0;
mm=0;
aa=0;
cnhist=zeros(1,cn);
while cni<cn
    cni=cni+1;
    if tg(1,cni)>tpre1 && tg(1,cni)<=t+0.000001
        m=0;
        while (tg(1,cni)<=t+0.000001)       
            m=m+1;
            mm=mm+1;
            hg=tg(1,cni)-tpre1; %Calculating the time step for the intermediate points
            xg=xpre1+hg*evaleval(xpre1,tpre1,egpre(1,cni))+((hg^2/h)*(xx-xpre1-(h*evaleval(xpre1,tpre1,egpre(1,cni))))); 
            xg=quz(xg,cq);
            f2=contim(eegm,egpre,G,T,u,xg,m,mm,cni); %This function contains the controllers' equations
            eeeg(1,mm)=f2;
            temp(1,cni)=tg(1,cni);
            tg(1,cni)=tg(1,cni)+T(1,cni);
            cnhist(mm,cni)=1;
        end
    end
end
if mm==0
        eeeg=0;
        maxl=0;
        cnhist=cnhist;
end


tg=tgpre; %reset the counter of the intermediate points to be counted again in the next iterations
end

%This function computes the mismatch for intermediate points regarding the
%controllers' signals (cni is used to varry the controller equation)
function f2=contim(eegm,egpre,G,T,u,xg,m,mm,cni)
if cni==1
            if m==1
                aw=egpre(1,cni)+G(1,cni)*T(1,cni)*(u-xg(2,1));
                aw=min(aw,0.98); %upper limnit anti wind up
                aw=max(aw,-0.98); %Downer limit anti wind up
                f2=eegm(mm,cni)-aw;
            else
                aw=+eegm(mm-1,cni)+G(1,cni)*T(1,cni)*(u-xg(2,1));
                aw=min(aw,0.98); %upper limnit anti wind up
                aw=max(aw,-0.98); %Downer limit anti wind up
                f2=eegm(mm,cni)-aw;            
            end
elseif cni==2
            if m==1
                aw=egpre(1,cni)+G(1,cni)*T(1,cni)*(u-xg(2,1));
                aw=min(aw,0.98); %upper limnit anti wind up
                aw=max(aw,-0.98); %Downer limit anti wind up
                f2=eegm(mm,cni)-aw;
            else
                aw=+eegm(mm-1,cni)+G(1,cni)*T(1,cni)*(u-xg(2,1));
                aw=min(aw,0.98); %upper limnit anti wind up
                aw=max(aw,-0.98); %Downer limit anti wind up
                f2=eegm(mm,cni)-aw;            
            end
elseif cni==3
             if m==1
                aw=egpre(1,cni)+G(1,cni)*T(1,cni)*(u-xg(2,1));
                aw=min(aw,0.98); %upper limnit anti wind up
                aw=max(aw,-0.98); %Downer limit anti wind up
                f2=eegm(mm,cni)-aw;
            else
                aw=+eegm(mm-1,cni)+G(1,cni)*T(1,cni)*(u-xg(2,1));
                aw=min(aw,0.98); %upper limnit anti wind up
                aw=max(aw,-0.98); %Downer limit anti wind up
                f2=eegm(mm,cni)-aw;            
            end
end
end

%Explicit interpolator
function [xx0g,ee0gm,ee0gv,eg,egt,egpret,tg,sg]=predictorg(xx,xpre1,t,tpre1,tpre2,holdd,T,eg,egpre,tg,cp,G,u,cn,cq)
tgpre=tg;  
cni=0;
mm=0;
while cni<cn
    cni=cni+1;
    if tg(1,cni)>tpre1 && tg(1,cni)<=t+0.000001
        m=0;
        while (tg(1,cni)<=t+0.000001)
            mm=mm+1;
            m=m+1;
            hg=tg(1,cni)-tpre1;
            [alpha0,alpha1,alpha2,betha1,betha2]=coefs(cp,hg,holdd);
            xg=xx+hg*(((betha1)*evaleval(xx,tpre1,eg(1,cni)))+((betha2)*evaleval(xpre1,tpre2,egpre(1,cni)))); 
            xx0g(mm,:)=xg;
            
            eg(1,cni)=contex(eg,G,T,u,xg,cni); %This function contains the controllers' equations
            ee0gm(mm,cni)=eg(1,cni);
            ee0gv(mm,1)=eg(1,cni);
            temp(1,cni)=tg(1,cni);
            tg(1,cni)=tg(1,cni)+T(1,cni);
        end     
    end
end
if mm~=0
    xx0g=quz(xx0g,cq);
    ee0gm=quz(ee0gm,cq);
    ee0gv=quz(ee0gv,cq);
    xx0g=xx0g(:,2);
    sg=mm;
    egt=sum(eg);
    egpret=sum(egpre);
else
    xx0g=0;
    tg=tg;
    eg=eg;
    egt=sum(eg);
    egpret=sum(egpre);
    ee0gm=0;
    ee0gv=0;
    sg=0;
end
tg=tgpre; %reset the counter of the intermediate points to be counted again in the next iterations
end

%This function computes the intermediate points regarding the
%controllers' signals
function egtemp=contex(eg,G,T,u,xg,cni)
if cni==1
    egtemp=eg(1,cni)+G(1,cni)*T(1,cni)*(u-xg(2,1));
    egtemp=min(egtemp,0.98); %upper limnit anti wind up
    egtemp=max(egtemp,-0.98); %Downer limit anti wind up
elseif cni==2
    egtemp=eg(1,cni)+G(1,cni)*T(1,cni)*(u-xg(2,1));
    egtemp=min(egtemp,0.98); %upper limnit anti wind up
    egtemp=max(egtemp,-0.98); %Downer limit anti wind up
elseif cni==3
    egtemp=eg(1,cni)+G(1,cni)*T(1,cni)*(u-xg(2,1));
    egtemp=min(egtemp,0.98); %upper limnit anti wind up
    egtemp=max(egtemp,-0.98); %Downer limit anti wind up
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
    q=10;
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
a=-0.17;
b=0.9;
A=[2*a sqrt(a^2+b^2);-sqrt(a^2+b^2) 0];
B=[-sqrt(a^2+b^2) 0]';
end