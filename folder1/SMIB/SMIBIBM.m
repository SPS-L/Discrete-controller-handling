%% %% This file solves the synchronous machine using the interpolation based method (IBM)
%%%% Please set the parameters based on the guidance provided in the Readme.md file
%both torque and AVR (with new J)
tic
clear all
clc
% This parameter chooses the predictor method. It can take values 4 and 5. 4 is for the first-order Adams-Bashford, and 5 is for the second-order Adams-Bashford. 
cp=5; %Predictor Choice
% This parameter chooses the corrector method. It can take values 2 and 3. 2 is for the second-order BDF, and 3 is for the second-order Adams-Moulton. 
cc=3; %Corrector Choice
%This parameters chooses the quantization method
%1 uniform quantization, 2 Nonuniform quantization
cq=1; %Quantization method
save=0;
t=0; %Starting point of the integration
h=0.001; %The first time-step size 0.05
holdd=h;
hmin=0.001; %The minimum acceptable time-step size
hmax=1; %The maximum acceptable time-step size

tsim=18; %Time of simulation
%Initial state variables values
[ssi]=ini;
xx=ssi;
dae=[6 13];
cn=2; %Number of controllers
T=[0.2 0.04]; %Time of sampling by the controllers
Tm=0.5013;
Vf=8.638587815812726e-04;
avrov=[1;0;0;Vf;Vf;Vf;0;0];
% avrov=[0;0;0;0;0;Vf;0;0];

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
egpre=[Tm,Vf];
eg=[Tm,Vf];
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
    
    %manual time control
    if t>1 && tpre1<1
        t=1;
        h=t-tpre1;
    end
    if t>1.2 && tpre1<1.2
        t=1.2;
        h=t-tpre1;
    end
    
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
        egpre2=[Tm,Vf];
    end
    if k>1
        egpre=eg;
        avrovpre=avrov;
    else
        egpre=[Tm,Vf];
        avrovpre=[1;0;0;Vf;Vf;Vf;0;0];
    end
    
    k=k+1;
    times(k)=t;
    tghold=tg;

    [alpha0,alpha1,alpha2,betha1,betha2]=coefs(cp,h,holdd); %calculating the coefficients of the predictor  
    [ee0gm,ee0gv,eg,tg,sg,avrovmm]=predictorg(xx,xpre1,xpre2,t,tpre1,tpre2,holdd,T,eg,egpre,tg,cp,cn,cq,avrov,dae); %predictor of the intermediate points (explicit predictor)
    xx0=predictor(xx,xpre1,eg,egpre,t,tpre1,tpre2,h,betha1,betha2,dae); %predictor
    [alpha0,alpha1,alpha2,betha1,betha2]=coefs(cc,h,holdd); %calculating the coefficients of the corrector
    [iterations,xx,eg,tg,feeg,eeg,avrov]=NewtonNL(h,holdd,xx0,xpre1,xpre2,alpha1,alpha2,alpha0,betha1,betha2,tpre1,tpre2,t,sg,eg,T,tg,ee0gm,ee0gv,egpre,cn,cq,dae,avrovpre,avrovmm); %The main function (Newton solver)
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
        hhist(4,1)=eg(1,1);
        hhist(5,1)=eg(1,2);
        hhist(6,1)=iterations;
            end
            hhist(1,k)=h;
            hhist(2,k)=dn;
            hhist(3,k)=r;
        hhist(4,k)=eg(1,1);
        hhist(5,k)=eg(1,2);
        hhist(6,k)=iterations;
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
            avrov=avrovpre;
            eg=egpre;
            egpre=egpre2;
            xx=xpre1;
            xpre1=xpre2;
            if k>2
                xpre2=xhist(1:s,k-2);
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
        if k==2
            hhist(1,1)=h;
            hhist(2,1)=dn;
            hhist(3,1)=r;
            hhist(4,1)=eg(1,1);
            hhist(5,1)=eg(1,2);
            hhist(6,1)=iterations;
        end
        hhist(1,k)=h;
        hhist(2,k)=dn;
        hhist(3,k)=r;
        hhist(4,k)=eg(1,1);
        hhist(5,k)=eg(1,2);
        hhist(6,k)=iterations;
        holdd=hpre;
    end

    

    feeghist=[feeghist feeg];
end
toc
% sgsize=size(eeghist);
% sgsize=sgsize(1,2);
% timessg=0.1:0.1:sgsize*0.1;


%Plotting the results
figure(1)
plot(times,sqrt(xhist(12,:).^2+xhist(13,:).^2),'+')
hold on
% plot(times,xhist(10,:),'+')
xlabel('Time (s)','FontSize', 24)
ylabel('Voltage (Per unit)','FontSize', 24)
set(gca,'FontSize',18)
% legend('V_{g} (SRM)','V_{g} (IBM)','Location','southeast')
legend('V_{g} (SRM)','V_{g} (IBM)','V_{g} (IBM AM1)','V_{g} (IBM BDF2)','Location','southeast')
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
    hgexport(gcf,'SMIBvoltage.pdf',figure_property);
end




figure(2)
plot(times,sqrt((xhist(16,:).^2)+(xhist(17,:).^2)),'+')
% plot(times,xhist(13,:),'+')
hold on
% plot(times,xhist(14,:),'+')
% plot(times,sqrt((xhist(13,:).^2)+(xhist(14,:).^2)))
xlabel('Time (s)','FontSize', 24)
ylabel('Current (Per unit)','FontSize', 24)
set(gca,'FontSize',18)
legend('I_{g} (SRM)','I_{g} (IBM)','Location','northeast')
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
    hgexport(gcf,'SMIBcurrent.pdf',figure_property);
end

figure(5)
plot(times,xhist(10,:),'--*')
hold on
plot(times,xhist(11,:),'--*')
xlabel('Time','FontSize', 24)
ylabel('vd & vq','FontSize', 24)
set(gca,'FontSize',18)

figure(3)
plot(times,xhist(3,:),'+')
hold on
xlabel('Time (s)','FontSize', 24)
ylabel('Speed deviation (Per unit)','FontSize', 24)
set(gca,'FontSize',18)
legend('SRM','IBM','Location','southeast')
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
    hgexport(gcf,'SMIBspeed.pdf',figure_property);
end

figure(4)
plot(times,xhist(4,:),'+')
hold on
xlabel('Time (s)','FontSize', 24)
ylabel('Rotor angle (Radian)','FontSize', 24)
set(gca,'FontSize',18)
legend('SRM','IBM','Location','northeast')
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
    hgexport(gcf,'SMIBdelta.pdf',figure_property);
end

figure(6)
plot(times,hhist(4,:),'+')
% stairs(times,hhist(4,:))
hold on
xlabel('Time (s)','FontSize', 24)
ylabel('Governor output (Per unit)','FontSize', 24)
set(gca,'FontSize',18)
% legend('SRM','IBM','Location','northeast')
% legend('SRM','IBM','IBM AM1','IBM BDF2','Location','northeast')
legend('Without quantization (IBM)','With quantization (IBM)','Location','southeast')
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
%     hgexport(gcf,'SMIBtorquesignal.pdf',figure_property);
    hgexport(gcf,'SMIBIBMQuantizationTest.pdf.pdf',figure_property);
end

figure(7)
plot(times,hhist(5,:),'+')
hold on
xlabel('Time (s)','FontSize', 24)
ylabel('Exciter output (Per unit)','FontSize', 24)
set(gca,'FontSize',18)
legend('SRM','IBM','Location','northeast')
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
    hgexport(gcf,'SMIBavrsignal.pdf',figure_property);
end

figure (8)
plot(times,sqrt((xhist(12,:).^2)+(xhist(13,:).^2)),'--*')
hold on
xlabel('Time (s)','FontSize', 24)
ylabel('Vter','FontSize', 24)
set(gca,'FontSize',18)

figure(9)
plot(times,xhist(7,:),'+')
hold on
xlabel('Time (s)','FontSize', 24)
ylabel('Frequency (Per unit)','FontSize', 24)
set(gca,'FontSize',18)
legend('SRM','IBM','Location','southeast')
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
    hgexport(gcf,'SMIBfrequency.pdf',figure_property);
end

figure(11)
plot(times,hhist(1,:))
hold on
xlabel('Time (s)','FontSize', 24)
ylabel('Step size (s)','FontSize', 24)
set(gca,'FontSize',18)
% legend('SRM','IBM','Location','northwest')
legend('Without quantization (IBM)','With quantization (IBM)','Location','northwest')
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
%     hgexport(gcf,'SMIBstepsize.pdf',figure_property);
    hgexport(gcf,'SMIBIBMQuantizationstepsize.pdf',figure_property);
end

figure(13)
plot(times,hhist(6,:))
xlabel('Time','FontSize', 24)
ylabel('Iterations','FontSize', 24)
set(gca,'FontSize',18)
hold on

figure(20)
plot(times,xhist(12,:).*xhist(16,:)+xhist(13,:).*xhist(17,:),'+')
hold on
xlabel('Time','FontSize', 24)
ylabel('Active power (Per unit)','FontSize', 24)
set(gca,'FontSize',18)
legend('SRM','IBM','Location','northeast')
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
    hgexport(gcf,'SMIBP.pdf',figure_property);
end
hold on
figure(21)
plot(times,xhist(13,:).*xhist(16,:)-xhist(12,:).*xhist(17,:),'+')
hold on
xlabel('Time','FontSize', 24)
ylabel('Reactive power (Per unit)','FontSize', 24)
set(gca,'FontSize',18)
legend('SRM','IBM','Location','northeast')
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
    hgexport(gcf,'SMIBRP.pdf',figure_property);
end


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

function [iterations,xx,eg,tg,feeg,eegv,avrov]=NewtonNL(h,holdd,xx,xpre1,xpre2,alpha1,alpha2,alpha0,betha1,betha2,tpre1,tpre2,t,sg,eg,T,tg,eegm,eegv,egpre,cn,cq,dae,avrovpre,avrovmm)
NTOL=1e-4; %Tolerance for the Newton solver
iterations=0;
a=1;
    while (iterations<10) && (a>NTOL)
        xxpre=xx;
        iterations=iterations+1;
        s=size(xx);
        s=s(1);
        if sg==0
            xxx=xx;
        else
            xxx=[xx;eegv]; %extending the state variables for controller intermediate signals
        end
        xxxpre=xxx;
%         [JJ]=Jacob2(h,xx,xpre1,xpre2,alpha1,alpha2,alpha0,betha1,betha2,tpre1,eg,egpre,t,dae,sg,eegm,cn); %calculating the Jacobian
        [J]=Jacob(h,xx,xpre1,xpre2,alpha1,alpha2,alpha0,betha1,betha2,tpre1,eg,egpre,t,dae); %calculating the Jacobian
        JJ=Jadd(J,sg); %Extending the Jacobian for intermediate points
        [feeg,tg,cnhist,avrovmm]=funcg(xx,xpre1,xpre2,eg,t,tpre1,tpre2,T,eegm,eegv,egpre,tg,h,holdd,cn,cq,avrovpre,avrovmm,dae); %Corrector of intermediate points (Implicit corrector)
        fy=func(h,xx,xpre1,xpre2,alpha1,alpha2,alpha0,betha1,betha2,tpre1,eg,egpre,t,dae); %Corrector
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
        if sg~=0
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
%             
% %             %calculating egt scalar
% %             egt=sum(eg);
         end

%         a=abs((xxx-xxxpre)./xxx);
        a=abs(xxx-xxxpre);
%         a=abs((xx-xxpre)./xx);
        a=max(a);
        if iterations==10
            disp('Did not converged')
            disp(t)
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
avrov=avrovmm(:,end);
end

% This function computes the Jacobian
function [J]=Jacob(h,xx,xpre1,xpre2,alpha1,alpha2,alpha0,betha1,betha2,tpre1,eg,egpre,t,dae)
s=size(xx);
s=s(1);
delta=1e-6;
i=0;
J=zeros(s,s);
while i<s
    i=i+1;
    f0=func(h,xx,xpre1,xpre2,alpha1,alpha2,alpha0,betha1,betha2,tpre1,eg,egpre,t,dae);
    yy=xx;
    yy(i)=xx(i)+delta;
    f1=func(h,yy,xpre1,xpre2,alpha1,alpha2,alpha0,betha1,betha2,tpre1,eg,egpre,t,dae);
    J(:,i)=(f1-f0)/delta;
end
end

% This function computes the Jacobian (New Jacobian)
function [J]=Jacob2(h,xx,xpre1,xpre2,alpha1,alpha2,alpha0,betha1,betha2,tpre1,eg,egpre,t,dae,sg,eegm,cn)
s=size(xx);
s=s(1);
delta=1e-6;
i=0;
J1=zeros(s,s);
J2=zeros(s,sg);
fiter=0;
while i<s
    i=i+1;
    f0=func(h,xx,xpre1,xpre2,alpha1,alpha2,alpha0,betha1,betha2,tpre1,eg,egpre,t,dae);
    fiter=fiter+1;
    yy=xx;
    yy(i)=xx(i)+delta;
    f1=func(h,yy,xpre1,xpre2,alpha1,alpha2,alpha0,betha1,betha2,tpre1,eg,egpre,t,dae);
    fiter=fiter+1;
    J1(:,i)=(f1-f0)/delta;
end
if sg==0
    J=J1;
else
    cni=0;
    while cni<cn
        cni=cni+1;
        if size(eegm)~=0    
            cc(1,cni)=nnz(eegm(:,cni)); %Number of intermediate points for every controller
            if cni==1
                b(1,cni)=nnz(eegm(:,cni)); %Number of intermediate points for both controllers (first controller)
            else
                b(1,cni)=nnz(eegm(:,cni))+b(1,cni-1); %Number of intermediate points for both controllers (adding up)
            end
        end
    end
    i=1;
    j=0;
        while i<sg
        j=j+1;
            while i<b(1,j)+1
                eg1=eg;
                eg2=eg;
                eg1(1,j)=eegm(i,j);
                eg2(1,j)=eegm(i,j)+delta;            
                f0=func(h,xx,xpre1,xpre2,alpha1,alpha2,alpha0,betha1,betha2,tpre1,eg1,egpre,t,dae);
                fiter=fiter+1;
                f1=func(h,xx,xpre1,xpre2,alpha1,alpha2,alpha0,betha1,betha2,tpre1,eg2,egpre,t,dae);
                fiter=fiter+1;
                J2(:,i)=(f1-f0)/delta;
                i=i+1;
            end
        end
        if cc(1,1)>1
            ii=0;
            while ii<cc(1,1)-1
               ii=ii+1;
               J2(:,ii)=0;
            end
        end
        if cc(1,2)>1
            ii=0+cc(1,1);
            while ii<cc(1,2)+cc(1,1)-1
               ii=ii+1;
               J2(:,ii)=0;
            end
        end
    J=[J1 J2;zeros(sg,s(1)) eye(sg)];
end
end

% This function extends the Jacobian for intermediate points
function JJ=Jadd(J,sg)
sj=size(J);
sj=sj(1);
JJ=[J zeros(sj(1),sg);zeros(sg,sj(1)) eye(sg)];
end

% This function determines the predictor formulation
function [xxp]=predictor(xx,xpre1,eg,egpre,t,tpre1,tpre2,h,betha1,betha2,dae)
    xxp1=xx+h*(((betha1)*evaleval(xx,t,eg))+((betha2)*evaleval(xpre1,tpre1,egpre)));
    xxp2=xx+evaleval(xx,t,eg);

if dae(1,1)==0
    xxp=xxp2;
elseif dae(1,2)==0
    xxp=xxp1;
else
    xxp=[xxp1(1:dae(1,1),1);xxp2(dae(1,1)+1:end,1)];
end
end


%This function determines the corrector formulation
function [f]=func(h,xx,xpre1,xpre2,alpha1,alpha2,alpha0,betha1,betha2,tpre1,eg,egpre,t,dae) 
    f1=h*(betha1*evaleval(xx,t,eg)+betha2*evaleval(xpre1,tpre1,egpre))+(alpha0*xx)+(alpha1*xpre1)+(alpha2*xpre2);
    f2=evaleval(xx,t,eg);

if dae(1,1)==0
    f=f2;
elseif dae(1,2)==0
    f=f1;
else
    f=[f1(1:dae(1,1),1);f2(dae(1,1)+1:end,1)];
end
end

%Implicit interpolator
function [eeeg,tg,cnhist,avrovmm]=funcg(xx,xpre1,xpre2,eg,t,tpre1,tpre2,T,eegm,eegv,egpre,tg,h,holdd,cn,cq,avrovpre,avrovmm,dae)
tgpre=tg;
cni=0;
mm=0;
avrc=0;
cnhist=zeros(1,cn);
while cni<cn
    cni=cni+1;
    if tg(1,cni)>tpre1 && tg(1,cni)<=t+0.000001
        m=0;
        while (tg(1,cni)<=t+0.000001)       
            m=m+1;
            mm=mm+1;
            hg=tg(1,cni)-tpre1; %Calculating the time step for the intermediate points
            xg=xpre1+hg*evaleval(xpre1,tpre1,egpre)+((hg^2/h)*(xx-xpre1-(h*evaleval(xpre1,tpre1,egpre)))); %Interpolator Taylor expansion
%         xg=hg*(evaleval(xx,t,eg))+(xpre1); %First order AM implicit interpolator
%        xg=(hg*(evaleval(xx,t,eg))+((1+hg/holdd)*xpre1)+((-(hg^2)/((hg+holdd)*holdd))*xpre2))/((2*hg+holdd)/(hg+holdd)); %Second order BDF implicit interpolator
            xg=quz(xg,cq);
            [f2,avrov]=contim(eegm,egpre,T,xg,m,mm,cni,avrovpre,avrovmm); %This function contains the controllers' equations
            eeeg(1,mm)=f2;
            temp(1,cni)=tg(1,cni);
            tg(1,cni)=tg(1,cni)+T(1,cni);
            cnhist(mm,cni)=1;
            if cni==2
               avrovmm(:,m)=avrov;
               avrc=avrc+1;
            end
        end
    end
end
if mm==0
        eeeg=0;
        maxl=0;
        cnhist=cnhist;
end
% if avrc==0
%     avrovmm=avrov;
% end

tg=tgpre; %reset the counter of the intermediate points to be counted again in the next iterations
end

%This function computes the mismatch for intermediate points regarding the
%controllers' signals (cni is used to varry the controller equation)
function [f2,avrov]=contim(eegm,egpre,T,xg,m,mm,cni,avrovpre,avrovmm)
if cni==1
            if m==1
                aw=torque(egpre(1,cni),T(1,cni),xg);
                aw=max(aw,0); %upper limnit anti wind up
                aw=min(aw,1); %lower limnit anti wind up
                f2=eegm(mm,cni)-aw;
                avrov=0;
            else
                aw=torque(eegm(mm-1,cni),T(1,cni),xg);
                aw=max(aw,0); %upper limnit anti wind up
                aw=min(aw,1); %lower limnit anti wind up
                f2=eegm(mm,cni)-aw;
                avrov=0;
            end
elseif cni==2
            if m==1
                [aw,avrov]=avr(avrovpre,T(1,cni),xg); %Controller's equation
                aw=max(aw,0); %upper limnit anti wind up
                aw=min(aw,1); %lower limnit anti wind up
                f2=eegm(mm,cni)-aw;
            else
                [aw,avrov]=avr(avrovmm(:,m-1),T(1,cni),xg); %Controller's equation
                aw=max(aw,0); %upper limnit anti wind up
                aw=min(aw,1); %lower limnit anti wind up
                f2=eegm(mm,cni)-aw; 
            end
end
end

%Explicit interpolator
function [ee0gm,ee0gv,eg,tg,sg,avrovmm]=predictorg(xx,xpre1,xpre2,t,tpre1,tpre2,holdd,T,eg,egpre,tg,cp,cn,cq,avrov,dae)
tgpre=tg;  
cni=0;
mm=0;
avrc=0;
while cni<cn
    cni=cni+1;
    if tg(1,cni)>tpre1 && tg(1,cni)<=t+0.000001
        m=0;
        while (tg(1,cni)<=t+0.000001)
            mm=mm+1;
            m=m+1;
            hg=tg(1,cni)-tpre1;
            [alpha0,alpha1,alpha2,betha1,betha2]=coefs(cp,hg,holdd);
            xg=xx+hg*(((betha1)*evaleval(xx,tpre1,eg))+((betha2)*evaleval(xpre1,tpre2,egpre)));
            xx0g(mm,:)=xg;
            
            [eg(1,cni),avrov]=contex(eg,T,xg,cni,avrov); %This function contains the controllers' equations
            ee0gm(mm,cni)=eg(1,cni);
            ee0gv(mm,1)=eg(1,cni);
            temp(1,cni)=tg(1,cni);
            tg(1,cni)=tg(1,cni)+T(1,cni);
            if cni==2
                avrc=avrc+1;
                avrovmm(:,m)=avrov;
            end
        end     
    end
end
if avrc==0
    avrovmm=avrov;
end
if mm~=0
%     xx0g=quz(xx0g,cq);
    ee0gm=quz(ee0gm,cq);
    ee0gv=quz(ee0gv,cq);
%     xx0g=xx0g(:,2);
    sg=mm;
else
%     xx0g=0;
    tg=tg;
    eg=eg;
    ee0gm=0;
    ee0gv=0;
    sg=0;
end
tg=tgpre; %reset the counter of the intermediate points to be counted again in the next iterations
end

%This function computes the intermediate points regarding the
%controllers' signals
function [egtemp,avrov]=contex(eg,T,xg,cni,avrov)
if cni==1
    egtemp=torque(eg(1,cni),T(1,cni),xg);
    egtemp=max(egtemp,0); %upper limnit anti wind up
    egtemp=min(egtemp,1); %lower limnit anti wind up
elseif cni==2
    [egtemp,avrov]=avr(avrov,T(1,cni),xg); %Controller's equation
    egtemp=max(egtemp,0); %upper limnit anti wind up
    egtemp=min(egtemp,1); %lower limnit anti wind up
end
end

function et2=torque(et2,T,xg)
omega=xg(3,1);
%Data
Tsm=0.9;
sigma=0.04;
zo=0.5013;
wS=1+omega;
%initialization (old values)
T2=T;
x9=et2;
%equations
% x9=(x9*(Tsm*sigma*wN-sigma*wN*T2)+T2*(zo*sigma*wN-wS+wN))/(Tsm*sigma*wN);
x9=(x9*(Tsm*sigma-sigma*T2)+T2*(zo*sigma-wS+1))/(Tsm*sigma);
%assigning old values for the next sample
et2=x9;
end

function [et1,avrov]=avr(avrov,T1,xg)
%Data
Tmm=0.2;
Te=0.12;
vo=1;
kp=0.003;
ki=0.004;
Ge=1;
kf=0.8;
Tf=0.9;
Vter=sqrt(xg(12,1)^2+xg(13,1)^2);
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
x8=(Tf*x8-T1*x8+T1*x7)/Tf;
et1=x62;
%assigning old values for the next sample
avrov=[x1 x2 x3 x4 x5 x62 x7 x8]';
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

%This function defines the system's formulation
function [f]=evaleval(xx,t,et)
global Tm omega Vf;
j=sqrt(-1);
TTT=0.05;
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
Le=Xe;
Vf=et(1,2); %8.638587815812726e-04
Tm=et(1,1); %5013
%Assiginng variables
psif=xx(1,1);
psiq1=xx(2,1);
omega=xx(3,1);
deltas=xx(4,1);
vxm=xx(5,1);
vym=xx(6,1);
freq=xx(7,1);
ifd=xx(8,1);
iq1=xx(9,1);
vd=xx(10,1);
vq=xx(11,1);
vx=xx(12,1);
vy=xx(13,1);
id=xx(14,1);
iq=xx(15,1);
ix=xx(16,1);
iy=xx(17,1);
psid=xx(18,1);
psiq=xx(19,1);
Te=psid*iq-psiq*id;
%Differential Equations
f1=wN*(-Rf*ifd+Vf); %psif (INPUT=Vf)
f2=wN*(-Rq1*iq1); %psiq1
f3=(Tm-Te-0.9*omega)/(2*HH); %omega
f4=omega; %deltas

%Frequency g
f5=(vx-vxm)/TTT;
f6=(vy-vym)/TTT;
%frequency 7
f7=1+(((vy-vym)*vxm)-((vx-vxm)*vym))/(2*pi*fN*TTT*(vxm^2+vym^2))-freq;




%Algebraic equations
f8=Ldd*id+Ldf*ifd-psid; %psid
f9=Lqq*iq+Lqq1*iq1-psiq; %psiq
f10=Lff*ifd+Ldf*id-psif; %psif
f11=Lq1q1*iq1+Lqq1*iq-psiq1; %psiq1
f12=vd+Ra*id+psiq; %vd
f13=vq+Ra*iq-psid; %vq
f14=vx-ex-Re*ix+Xe*iy; %vx
f15=vy-ey-Re*iy-Xe*ix; %vy
f16=vd-vx*cos(theta0)-vy*sin(theta0); %vd
f17=vq-vx*sin(theta0)+vy*cos(theta0); %vq
f18=id-ix*cos(theta0)-iy*sin(theta0); %id
f19=iq-ix*sin(theta0)+iy*cos(theta0); %iq
%sum
f=[f1 f2 f3 f4 f5 f6 f7 f8 f9 f10 f11 f12 f13 f14 f15 f16 f17 f18 f19]';
end

%initial values
function [ssi]=ini
global Tm omega Vf;
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
Tm=psid*iq-psiq*id;
Tm2=P+Ra*It^2;
iq1=0;
psiq1=Lqq1*iq;
omega=0;
deltas=deltai;
%Output
ssi=[psif psiq1 omega deltas vx vy 1 ifd iq1 vd vq vx vy id iq ix iy psid psiq]';
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