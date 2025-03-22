%% %% This file solves PV using the traditional method (reference)
%%%% Please set the parameters based on the guidance provided in the Readme.md file 
%both torque and AVR
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
hmax=0.500; %The maximum acceptable time-step size

tsim=16; %Time of the simulation
%Initial state variables values
[ssi]=ini;
xx=ssi;
dae=[4 22];
cn=2;
T=[0.2 0.05]; %Time of sampling by the controller
%Input
Tm=0.606760138846862;
Vf=0.000996515736969354;
avrov=[1;0;0;Vf;Vf;Vf;0;0];
G=0.05;
u=1; %Input 

s=size(xx);
s=s(1);
xpre1=xx;
xpre2=xx;
times(1)=t;
xhist(1:s,1)=xx; %Allocating a space to save the results
hhist(1:5,1)=[h;0;0;Tm;Vf];

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
        et=[Tm,Vf];
        etpre=[Tm,Vf];
        etpre2=[Tm,Vf];

save=1;
while t<tsim
    t=t
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
%     if t>3 && tpre1<3
%         t=3;
%         h=t-tpre1;
%     end
% %     if t>6 && tpre1<6
% %         t=6;
% %         h=t-tpre1;
% %     end
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
        et=[Tm,Vf];
        etpre=[Tm,Vf];
        etpre2=[Tm,Vf];
        avrov=[1;0;0;Vf;Vf;Vf;0;0];
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
            xxq=quz(xx,cq);
            [et,avrov]=contex(conti,et,T,u,cn,avrov,xxq); %Calculating the controllers' signal in case of any event
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
plot(times,xhist(9,:),'--')
hold on
plot(times,xhist(10,:),'--')
% plot(times,sqrt((xhist(9,:).^2)+(xhist(10,:).^2)))
xlabel('Time','FontSize', 24)
ylabel('Voltage','FontSize', 24)
set(gca,'FontSize',18)

figure(2)
plot(times,xhist(13,:),'--')
hold on
plot(times,xhist(14,:),'--')
% plot(times,sqrt((xhist(13,:).^2)+(xhist(14,:).^2)))
xlabel('Time','FontSize', 24)
ylabel('Current','FontSize', 24)
set(gca,'FontSize',18)

figure(5)
plot(times,xhist(7,:),'--')
hold on
plot(times,xhist(8,:),'--')
xlabel('Time','FontSize', 24)
ylabel('vd & vq','FontSize', 24)
set(gca,'FontSize',18)

figure(3)
plot(times,xhist(3,:),'--')
hold on
xlabel('Time (s)','FontSize', 24)
ylabel('Speed deviation (Per unit)','FontSize', 24)
set(gca,'FontSize',18)
legend('Analog fixed time steps','Analog variable time steps','Digital fixed time steps','Digital variable time steps','Location','southeast')
set(gca, 'FontName', 'Times New Roman')
xlim([0 tsim])
box on
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
    hgexport(gcf,'C:\Users\mjafa\Desktop\Coding tutorial\MATLAB\PI controller paper\3BusSpeedDeviations.pdf',figure_property);
end

figure(4)
plot(times,xhist(4,:),'--')
hold on
xlabel('Time (s)','FontSize', 24)
ylabel('Rotor angle (Radian)','FontSize', 24)
set(gca,'FontSize',18)
legend('Analog fixed time steps','Analog variable time steps','Digital fixed time steps','Digital variable time steps','Location','northeast')
set(gca, 'FontName', 'Times New Roman')
xlim([0 tsim])
box on
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
    hgexport(gcf,'C:\Users\mjafa\Desktop\Coding tutorial\MATLAB\PI controller paper\3BusRotorAngle.pdf',figure_property);
end


figure(6)
plot(times,hhist(4,:),'--')
hold on
xlabel('Time (s)','FontSize', 24)
ylabel('Governor output','FontSize', 24)
set(gca,'FontSize',18)
legend('Analog fixed time steps','Analog variable time steps','Digital fixed time steps','Digital variable time steps','Location','northeast')
set(gca, 'FontName', 'Times New Roman')
xlim([0 tsim])
box on
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
    hgexport(gcf,'C:\Users\mjafa\Desktop\Coding tutorial\MATLAB\PI controller paper\3BusGovernor.pdf',figure_property);
end


figure(7)
plot(times,hhist(5,:),'--')
hold on
xlabel('Time (s)','FontSize', 24)
ylabel('Exciter output','FontSize', 24)
set(gca,'FontSize',18)
legend('Analog fixed time steps','Analog variable time steps','Digital fixed time steps','Digital variable time steps','Location','northeast')
set(gca, 'FontName', 'Times New Roman')
xlim([0 tsim])
box on
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
    hgexport(gcf,'C:\Users\mjafa\Desktop\Coding tutorial\MATLAB\PI controller paper\3BusExciter.pdf',figure_property);
end


figure (8)
% plot(times,sqrt((xhist(9,:).^2)+(xhist(10,:).^2))) $Bus 1 Generator
hold on
% plot(times,sqrt((xhist(17,:).^2)+(xhist(18,:).^2))) %Bus 2 PV
plot(times,sqrt((xhist(21,:).^2)+(xhist(22,:).^2)),'--')  %Bus 3 Load
xlabel('Time (s)','FontSize', 24)
ylabel('Voltage of bus 2 (per unit)','FontSize', 24)
set(gca,'FontSize',18)
legend('Analog fixed time steps','Analog variable time steps','Digital fixed time steps','Digital variable time steps','Location','northwest')
set(gca,'FontSize',18)
xlim([0 tsim])
box on
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
    hgexport(gcf,'C:\Users\mjafa\Desktop\Coding tutorial\MATLAB\PI controller paper\3BusVoltage.pdf',figure_property);
end


figure (9)
% plot(times,sqrt((xhist(13,:).^2)+(xhist(14,:).^2)))
hold on
% plot(times,sqrt((xhist(19,:).^2)+(xhist(20,:).^2)))
plot(times,sqrt((xhist(23,:).^2)+(xhist(24,:).^2)),'--') %Bus 3 load
xlabel('Time (s)','FontSize', 24)
ylabel('Load current (per unit)','FontSize', 24)
set(gca,'FontSize',18)
legend('Analog fixed time steps','Analog variable time steps','Digital fixed time steps','Digital variable time steps','Location','southeast')
set(gca, 'FontName', 'Times New Roman')
xlim([0 tsim])
box on
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
    hgexport(gcf,'C:\Users\mjafa\Desktop\Coding tutorial\MATLAB\PI controller paper\3BusCurrent.pdf',figure_property);
end


figure (101)
plot(times,xhist(9,:).*xhist(13,:)+xhist(10,:).*xhist(14,:),'--')
hold on
plot(times,xhist(17,:).*xhist(19,:)+xhist(18,:).*xhist(20,:),'--')
plot(times,xhist(21,:).*xhist(23,:)+xhist(24,:).*xhist(22,:),'--')
xlabel('Time','FontSize', 24)
ylabel('P','FontSize', 24)
set(gca,'FontSize',18)
figure (9)

figure (111)
plot(times,xhist(10,:).*xhist(13,:)-xhist(9,:).*xhist(14,:),'--')
hold on
plot(times,xhist(18,:).*xhist(19,:)-xhist(17,:).*xhist(20,:),'--')
plot(times,xhist(22,:).*xhist(23,:)-xhist(21,:).*xhist(24,:),'--')
xlabel('Time','FontSize', 24)
ylabel('Q','FontSize', 24)
set(gca,'FontSize',18)
figure (9)
% 
% figure(4)
% plot(times,xhist(9,:))
% hold on
% xlabel('Time','FontSize', 24)
% ylabel('deltas','FontSize', 24)
% set(gca,'FontSize',18)

% figure(2)
% plot(xhist(1,:),xhist(3,:))

figure(11)
plot(times,hhist(1,:))
hold on
xlabel('Time (s)','FontSize', 24)
ylabel('Step Size (s)','FontSize', 24)
set(gca,'FontSize',18)
legend('Analog variable time steps','Digital variable time steps','Location','southeast')
set(gca, 'FontName', 'Times New Roman')
xlim([0 tsim])
box on
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
    hgexport(gcf,'C:\Users\mjafa\Desktop\Coding tutorial\MATLAB\PI controller paper\3BusTimeStep.pdf',figure_property);
end


figure(12)
plot(times,hhist(2,:))
xlabel('Time','FontSize', 24)
ylabel('Error Estimate','FontSize', 24)
set(gca,'FontSize',18)
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
function [etj,avrov]=contex(conti,etj,T,u,cn,avrov,xx)
jj=size(conti);
jj=jj(1,2);
j=0;
while j<jj
    j=j+1;
    contj=conti(1,j);
    if contj==1
        etj(1,contj)=torque(etj(1,contj),T(1,1),xx);
%         etj(1,contj)=max(etj(1,contj),0); %upper limnit anti wind up
%         etj(1,contj)=min(etj(1,contj),1); %lower limnit anti wind up
%         etj(1,contj)=min(etj(1,contj),1); %lower limnit anti wind up
    elseif contj==2
        [etj(1,contj),avrov]=avr(avrov,T(1,2),xx); %Controller's equation
%         etj(1,contj)=max(etj(1,contj),0); %upper limnit anti wind up
%         etj(1,contj)=min(etj(1,contj),0.01); %lower limnit anti wind up
%         etj(1,contj)=min(etj(1,contj),0.01); %lower limnit anti wind up
    end
end
end

function et2=torque(et2,T,xx)
% global Tm wN omega;
%Data
Tsm=1; %0.3
sigma=0.04;
zo=0.606760138846862;
wS=1+xx(3,1);
%initialization (old values)
T2=T;
x9=et2;
%equations
% x9=(x9*(Tsm*sigma*wN-sigma*wN*T2)+T2*(zo*sigma*wN-wS+wN))/(Tsm*sigma*wN);
x9=(x9*(Tsm*sigma-sigma*T2)+T2*(zo*sigma-wS+1))/(Tsm*sigma);
%assigning old values for the next sample
et2=x9;
end

function [et1,avrov]=avr(avrov,T1,xx)
%Data
Tmm=0.2; %0.2
Te=0.12; %0.12
vo=1;
kp=0.005; %0.01
ki=0.01; %0.01
Ge=1;
kf=300; %0.1
Tf=0.9; %0.9
Vter=sqrt(xx(9,1)^2+xx(10,1)^2);
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

function [iterations,xx]=NewtonNL(h,xx,xpre1,xpre2,alpha1,alpha2,alpha0,betha1,betha2,tpre1,et,etpre,t,dae)
NTOL=1e-4; %Tolerance for the Newton solver
iterations=0;
a=1;
    while (iterations<10) && (a>NTOL)
    xxpre=xx;
    iterations=iterations+1;
    [J]=Jacob(h,xx,xpre1,xpre2,alpha1,alpha2,alpha0,betha1,betha2,tpre1,et,etpre,t,dae); %calculating the Jacobian
    bb=func(h,xx,xpre1,xpre2,alpha1,alpha2,alpha0,betha1,betha2,tpre1,et,etpre,t,dae);
    cc=J\bb;
    yy=xx-cc; %Newton
    xx=yy;
%     a=abs(xx-xxpre);
    a=abs((xx-xxpre)./xx);
    a=max(a);
    if iterations==10
        disp('Did not converged')
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
% global Tm omega Vf;
j=sqrt(-1);
%Parameters
pp=para;
Pgen=pp(1,1);
Qgen=pp(1,2);
Vgen=pp(1,3);
Xe1=pp(1,4);
Re1=pp(1,5);
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
Re2=pp(1,29);
Xe2=pp(1,30);
Rload=pp(1,31);
Xload=pp(1,32);
Ppv=pp(1,33);
Qpv=pp(1,34);
Vpv=pp(1,35);
deltapv=pp(1,36);
Vload=pp(1,37);
deltaload=pp(1,38);
Pload=pp(1,39);
Qload=pp(1,40);
G=0;
B=0;
% theta=theta0+wN*t;
%Input
% Vf=0.0010873;
Vf=et(1,2); %0.001196163166093
% Vf=0.0010873;
Tm=et(1,1); %0.709060632610000
% Tm=0.606766;
%Assiginng variables
psif=xx(1,1);
psiq1=xx(2,1);
omega=xx(3,1);
deltas=xx(4,1);
ifd=xx(5,1);
iq1=xx(6,1);
vd=xx(7,1);
vq=xx(8,1);
vx1=xx(9,1);
vy1=xx(10,1);
id=xx(11,1);
iq=xx(12,1);
ix1=xx(13,1);
iy1=xx(14,1);
psid=xx(15,1);
psiq=xx(16,1);
vx2=xx(17,1);
vy2=xx(18,1);
ix2=xx(19,1);
iy2=xx(20,1);
vx3=xx(21,1);
vy3=xx(22,1);
ix3=xx(23,1);
iy3=xx(24,1);
ix32=xx(25,1);
iy32=xx(26,1);
%Extra
% if t>=3 && t<=3.1
%     B=-0.8; 
%     G=0;
% end
if t>3 && t<6
    Ppv=0.2;
end
% if t>3 && t<13
%     Ppv=0.2;
% end
% phasIgen=(Pgen-j*Qgen)/Vgen;
% phi=angle(Vgen+Ra*phasIgen+j*Xq*phasIgen);
phasIgen=((vx1*ix1+vy1*iy1)-j*(vy1*ix1-vx1*iy1))/(vx1+j*vy1);
phi=angle(vx1+j*vy1+Ra*phasIgen+j*Xq*phasIgen);
theta0=phi+pi/2;
Te=psid*iq-psiq*id;
%%%Differential Equations
f1=wN*(-Rf*ifd+Vf); %psif (INPUT=Vf)
f2=wN*(-Rq1*iq1); %psiq1
f3=(Tm-Te-0.9*omega)/(2*HH); %omega
f4=omega; %deltas
%%%Algebraic equations
%flux-current
f5=Ldd*id+Ldf*ifd-psid; 
f6=Lqq*iq+Lqq1*iq1-psiq; 
f7=Lff*ifd+Ldf*id-psif; 
f8=Lq1q1*iq1+Lqq1*iq-psiq1;
%stator park
f9=vd+Ra*id+psiq; 
f10=vq+Ra*iq-psid;
%network line1
f11=vx1-vx3-Re1*ix1+Xe1*iy1; 
f12=vy1-vy3-Re1*iy1-Xe1*ix1;
%dqxy
f13=vd-vx1*cos(theta0)-vy1*sin(theta0); 
f14=vq-vx1*sin(theta0)+vy1*cos(theta0); 
f15=id-ix1*cos(theta0)-iy1*sin(theta0); 
f16=iq-ix1*sin(theta0)+iy1*cos(theta0);
%networkPV PV
f17=-Ppv+vx2*ix2+vy2*iy2;
f18=-Qpv-vx2*iy2+vy2*ix2;
%networkpv line2
f19=vx2-vx3-Re2*ix2+Xe2*iy2;
f20=vy2-vy3-Re2*iy2-Xe2*ix2;
%network load
f21=-Pload+vx3*ix3+vy3*iy3;
f22=-Qload-vx3*iy3+vy3*ix3;
% f25=-ix32+G*vx3-B*vy3;
% f26=-iy32+B*vx3+G*vy3;
f25=-ix32+G*vx3-B*vy3;
f26=-iy32+B*vx3+G*vy3;
%network KCL
f23=ix1+ix2-ix3-ix32;
f24=iy1+iy2-iy3-iy32;
%sum
f=[f1 f2 f3 f4 f5 f6 f7 f8 f9 f10 f11 f12 f13 f14 f15 f16 f17 f18 f19 f20 f21 f22 f23 f24 f25 f26]';
end

%initial values
function [ssi]=ini
global Tm omega Vf;
j=sqrt(-1);
%parameters
pp=para;
Pgen=pp(1,1);
Qgen=pp(1,2);
Vgen=pp(1,3);
Xe1=pp(1,4);
Re1=pp(1,5);
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
Re2=pp(1,29);
Xe2=pp(1,30);
Rload=pp(1,31);
Xload=pp(1,32);
Ppv=pp(1,33);
Qpv=pp(1,34);
Vpv=pp(1,35);
deltapv=pp(1,36);
Vload=pp(1,37);
deltaload=pp(1,38);
Pload=pp(1,39);
Qload=pp(1,40);
%%%Network
It=(sqrt(Pgen^2+Qgen^2))/Vgen; %Terminal Current
phi=acos(Pgen/(Vgen*It)); %Power Factor angle
deltai=atan((Xq*It*cos(phi)-Ra*It*sin(phi))/(Vgen+Ra*It*cos(phi)+Xq*It*sin(phi))); %Rotor Angle
vx1=Vgen;                                
vy1=0;
phasVpv=Vpv*(cos(deg2rad(deltapv))+j*sin(deg2rad(deltapv)));
vx2=real(phasVpv);                                
vy2=imag(phasVpv);
phasVload=Vload*(cos(deg2rad(deltaload))+j*sin(deg2rad(deltaload)));
vx3=real(phasVload);                                
vy3=imag(phasVload);
phasIgen=(Pgen-j*Qgen)/Vgen; %PV current (phasor)
ix1=real(phasIgen);
iy1=imag(phasIgen);
phasIpv=(Ppv-j*Qpv)/(conj(phasVpv)); %PV current (phasor)
ix2=real(phasIpv);
iy2=imag(phasIpv);
phasIload=(Pload-j*Qload)/(conj(phasVload)); %Load current (phasor)
ix3=real(phasIload);
iy3=imag(phasIload);
Zload=phasVload/-phasIload;
Rload=real(Zload);
Xload=imag(Zload);
%Machine
phi=angle(Vgen+Ra*phasIgen+j*Xq*phasIgen);     
theta0=phi+pi/2;
vd=Vgen*cos(theta0);
vq=Vgen*sin(theta0);
%id=abs(phasIgen)*cos(theta0-angle(phasIgen));
id=-0.5948;
%iq=abs(phasIgen)*sin(theta0-angle(phasIgen));
iq=0.2795;
%psid=vq+Ra*iq;
psid=0.7462;
%psiq=Lqq*iq;
psiq=0.6707;
ifd=(psid-Ldd*id)/Ldf;
psif=Lff*ifd+Ldf*id;
Vf=Rf*ifd;
%Tm=psid*iq-psiq*id;
Tm=0.7462*0.2795+0.6707*0.5948;
Tm2=Pgen+Ra*It^2;
iq1=0;
psiq1=Lqq1*iq;
omega=0;
deltas=deltai;
ix32=0;
iy32=0;
%Output
ssi=[0.994857646454475;0.670038962213594;2.34348578337354e-07;0.804636233538080;0.905559175933906;3.10080972807644e-09;-0.728232058986814;0.685330397459342;0.999999842708965;-2.60990946565615e-09;-0.543907072411626;0.304563161469621;0.604817054479558;-0.150963417271962;0.686853213266690;0.730951594348872;0.996759781735134;0.000401575503337228;0.300995386096771;-0.0500412722707348;0.993649745329625;0.000299997453851034;0.905812440576328;-0.201004689542696;0;0];
% ssi=[0.997143582747761 0.670393695105975 3.37790720235446e-07 0.804622750785925 0.904919168624512 -9.19460014048529e-05 -0.728651277537273 0.688104472710673 1.00220778462504 7.78061864817495e-05 -0.542164002419800 0.304818491146860 0.603474833688653 -0.150579690422287 0.689628565166408 0.731362097549372 0.998975009006049 0.000478241597140939 0.300331705144652 -0.0499075236475681 0.995871876907308 0.000376653423327316 0.903806538833305 -0.200487214069855 0 0]';
%ssi=[   1    2     3     4     5   6  7  8  9   10  11 12 13  14  15   16   17  18  19  20  21  22  23  24  25   26] 
% ssi=[psif psiq1 omega deltas ifd iq1 vd vq vx1 vy1 id iq ix1 iy1 psid psiq vx2 vy2 ix2 iy2 vx3 vy3 ix3 iy3 ix32 iy32]';
end

function pp=para
%%%DATA%%%
%Network
Pgen=6.0048/10;
Qgen=1.5958/10;
Vgen=1;
Vpv=0.9976;
deltapv=-0.3385;
Vload=0.9963;
deltaload=-0.6815;
Xe1=0.002;
Re1=0.01;
Xe2=0.002;
Re2=0.01;
Rload=1.051002730588235/10;
Xload=0.233556162352941/10;
Ppv=3/10;
Qpv=0.5/10;
Pload=9/10;
Qload=2/10;
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
Le=Xe1 ;
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
pp=[Pgen Qgen Vgen Xe1 Re1 Ra fN Xd Xq Xl Xpd Xpq Tpdo_s Tpqo_s wN TB Le Ldd Ldf Lff Tpdo Rf Lqq Lqq1 Lq1q1 Tpqo Rq1 HH Re2 Xe2 Rload Xload Ppv Qpv Vpv deltapv Vload deltaload Pload Qload];
end