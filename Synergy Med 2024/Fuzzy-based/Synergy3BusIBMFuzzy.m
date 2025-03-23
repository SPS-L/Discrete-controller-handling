%% %% This file solves the synchronous machine using the interpolation based method (IBM)
%%%% Please set the parameters based on the guidance provided in the Readme.md file
%both torque and AVR
aaa=0;
tsimave=0;
endtsimave=0;
while aaa<1
    aaa=aaa+1;
clearvars -except net fuzzyAVR3Bus tsimave aaa
clc
global fuzzyAVR3Bus;
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
h=0.001; %The first time-step size 0.05
holdd=h;
hmin=0.001; %The minimum acceptable time-step size
hmax=1; %The maximum acceptable time-step size

tsim=22; %Time of simulation
%Initial state variables values
[ssi]=ini;
xx=ssi;
dae=[6 22];
cn=2; %Number of controllers
T=[0.02 0.004]; %Time of sampling by the controllers
%Input
Tm=0.606781005859375; %0.607463
Vf=0.000991821289062500; %0.0010873
avrov=[1;0;0;Vf;Vf;Vf;0;0];
gvrov=Tm;


s=size(xx);
s=s(1);
xpre1=xx;
xpre2=xx;
times(1)=t;
xhist(1:s,1)=xx; %Allocating a space to save the results

timessg=0;
eeghist=0;
feeghist=0;
dhiterationsc=0;
iterationsc=0; %Allocating a space to save the number of Newoton iterations
tsc=0; %Number of time steps
mtsc=0; %Number of accepted time steps with the manimim value
egpre=[Tm,Vf];
eg=[Tm,Vf];
egt=0; %Summation of controllers' signal for the current time step
egpret=0; %Summation of controllers' signal for the previous time step
tg=T; %Time of the events
k=1;
ETOL=1e-4; %The tolerance to compare with the error estimate
fitertt=0;
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
    if t>3 && tpre1<3
        t=3;
        h=t-tpre1;
    end
%     if t>6 && tpre1<6
%         t=6;
%         h=t-tpre1;
%     end
    if t>3.1 && tpre1<3.1
        t=3.1;
        h=t-tpre1;
    end
%     if t>13 && tpre1<13
%         t=13;
%         h=t-tpre1;
%     end
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
        gvrovpre=gvrov;
    else
        egpre=[Tm,Vf];
        avrovpre=[1;0;0;Vf;Vf;Vf;0;0];
        gvrovpre=gvrov;
    end
    
    k=k+1;
    times(k)=t;
    tghold=tg;

    [alpha0,alpha1,alpha2,betha1,betha2]=coefs(cp,h,holdd); %calculating the coefficients of the predictor  
    [ee0gm,ee0gv,eg,tg,sg,avrovmm,gvrovmm]=predictorg(xx,xpre1,xpre2,t,tpre1,tpre2,holdd,T,eg,egpre,tg,cp,cn,cq,avrov,gvrov,dae); %predictor of the intermediate points (explicit predictor)
    xx0=predictor(xx,xpre1,eg,egpre,t,tpre1,tpre2,h,betha1,betha2,dae); %predictor
    [alpha0,alpha1,alpha2,betha1,betha2]=coefs(cc,h,holdd); %calculating the coefficients of the corrector
    [iterations,xx,eg,tg,feeg,eeg,avrov,gvrov,fitert,dhiterations]=NewtonNL(h,xx0,xpre1,xpre2,alpha1,alpha2,alpha0,betha1,betha2,tpre1,tpre2,t,sg,eg,T,tg,ee0gm,ee0gv,egpre,cn,cq,dae,avrovpre,avrovmm,gvrovpre,gvrovmm,holdd); %The main function (Newton solver)
    dhiterationsc=dhiterationsc+dhiterations;
    iterationsc=iterations+iterationsc; %Counting the number of Newton iterations
    fitertt=fitertt+fitert;
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
            end
            hhist(1,k)=h;
            hhist(2,k)=dn;
            hhist(3,k)=r;
        hhist(4,k)=eg(1,1);
        hhist(5,k)=eg(1,2);
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
            gvrov=gvrovpre;
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
        end
        hhist(1,k)=h;
        hhist(2,k)=dn;
        hhist(3,k)=r;
        hhist(4,k)=eg(1,1);
        hhist(5,k)=eg(1,2);
        holdd=hpre;
    end

    

    feeghist=[feeghist feeg];
end
tsimave=tsimave+toc;
end

endtsimave=tsimave/aaa;
% sgsize=size(eeghist);
% sgsize=sgsize(1,2);
% timessg=0.1:0.1:sgsize*0.1;


%Plotting the results
figure(1)
plot(times,xhist(11,:),'+')
hold on
plot(times,xhist(12,:),'+')
% plot(times,sqrt((xhist(9,:).^2)+(xhist(10,:).^2)),'+')
xlabel('Time (s)','FontSize', 24)
ylabel('Voltage (Per unit)','FontSize', 24)
set(gca,'FontSize',18)
legend('V_{x} (SRM)','V_{y} (SRM)','V_{x} (IBM)','V_{y} (IBM)','Location','southeast')
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
    hgexport(gcf,'voltagePV.pdf',figure_property);
end





figure(2)
plot(times,xhist(15,:),'+')
hold on
plot(times,xhist(16,:),'+')
% plot(times,sqrt((xhist(13,:).^2)+(xhist(14,:).^2)))
xlabel('Time (s)','FontSize', 24)
ylabel('Current (Per unit)','FontSize', 24)
set(gca,'FontSize',18)
legend('I_{x} (SRM)','I_{y} (SRM)','I_{x} (IBM)','I_{y} (IBM)','Location','southeast')
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
    hgexport(gcf,'currentPV.pdf',figure_property);
end

figure(5)
plot(times,xhist(9,:),'--*')
hold on
plot(times,xhist(10,:),'--*')
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
    hgexport(gcf,'speedPV.pdf',figure_property);
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
    hgexport(gcf,'deltaPV.pdf',figure_property);
end

figure(6)
plot(times,hhist(4,:),'*')
hold on
xlabel('Time (s)','FontSize', 24)
ylabel('Governor output','FontSize', 24)
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
    hgexport(gcf,'torquesignalPV.pdf',figure_property);
end

figure(7)
plot(times,hhist(5,:).*1387,'*')
hold on
xlabel('Time (s)','FontSize', 24)
ylabel('Excitation system output','FontSize', 24)
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
    hgexport(gcf,'avrsignalPV.pdf',figure_property);
end

figure (8)
plot(times,sqrt((xhist(11,:).^2)+(xhist(12,:).^2)),'*')
hold on
plot(times,sqrt((xhist(19,:).^2)+(xhist(20,:).^2)),'*')
plot(times,sqrt((xhist(23,:).^2)+(xhist(24,:).^2)),'*')
xlabel('Time (s)','FontSize', 24)
ylabel('Voltage (per unit)','FontSize', 24)
set(gca,'FontSize',18)
legend('Bus1 (SRM)','Bus3 (SRM)','Bus2 (SRM)','Bus1 (IBM)','Bus3 (IBM)','Bus2 (IBM)','Location','southeast')
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
    hgexport(gcf,'tvoltagePV.pdf',figure_property);
end

figure (9)
plot(times,sqrt((xhist(15,:).^2)+(xhist(16,:).^2)),'*')
hold on
plot(times,sqrt((xhist(21,:).^2)+(xhist(22,:).^2)),'*')
plot(times,sqrt((xhist(25,:).^2)+(xhist(26,:).^2)),'*')
xlabel('Time (s)','FontSize', 24)
ylabel('Current (per unit)','FontSize', 24)
set(gca,'FontSize',18)
legend('Generator (SRM)','PV (SRM)','Load (SRM)','Generator (IBM)','PV (IBM)','Load (IBM)','Location','southeast')
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
    hgexport(gcf,'tcurrentPV.pdf',figure_property);
end

figure(11)
plot(times,hhist(1,:))
hold on
xlabel('Time (s)','FontSize', 24)
ylabel('Step size (s)','FontSize', 24)
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
    hgexport(gcf,'stepsizePV.pdf',figure_property);
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

iterationsc=iterationsc
fitertt=fitertt
dhiterationsc=dhiterationsc
endtsimave=endtsimave


function [iterations,xx,eg,tg,feeg,eegv,avrov,gvrov,fitert,dhiterations]=NewtonNL(h,xx,xpre1,xpre2,alpha1,alpha2,alpha0,betha1,betha2,tpre1,tpre2,t,sg,eg,T,tg,eegm,eegv,egpre,cn,cq,dae,avrovpre,avrovmm,gvrovpre,gvrovmm,holdd)
NTOL=0.3e-4; %Tolerance for the Newton solver
iterations=0;
dhiterations=0;
a=1;
fitert=0;
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
        
        %Dishonest Jacobian
%         if iterations==1 || iterations>5
%             dhiterations=dhiterations+1;
%             %C
% %             [JJ,fiter]=Jacob4(h,xx,xpre1,xpre2,alpha1,alpha2,alpha0,betha1,betha2,tpre1,eg,egpre,t,dae,sg,eegm,eegv,cn,avrovpre,avrovmm,tg,tpre2,T,holdd,cq);
%             %B and C
%             [JJ,fiter]=Jacob3(h,xx,xpre1,xpre2,alpha1,alpha2,alpha0,betha1,betha2,tpre1,eg,egpre,t,dae,sg,eegm,eegv,cn,avrovpre,avrovmm,tg,tpre2,T,holdd,cq); %Exact Jacobian
%             %B
% %             [JJ,fiter]=Jacob2(h,xx,xpre1,xpre2,alpha1,alpha2,alpha0,betha1,betha2,tpre1,eg,egpre,t,dae,sg,eegm,cn); %Simplified Jacobian with B
%             %Simplified
% %             [JJ,fiter]=Jacob1(h,xx,xpre1,xpre2,alpha1,alpha2,alpha0,betha1,betha2,tpre1,eg,egpre,t,dae,sg); %Simplified Jacobian 
%             fitert=fitert+fiter;
%         end
        
        %Honest Jacobian
        %All
%       [JJ,fiter]=Jacob5(h,xx,xpre1,xpre2,alpha1,alpha2,alpha0,betha1,betha2,tpre1,eg,egpre,t,dae,sg,eegm,eegv,cn,avrovpre,avrovmm,tg,tpre2,T,holdd,cq);
        %C
%         [JJ,fiter]=Jacob4(h,xx,xpre1,xpre2,alpha1,alpha2,alpha0,betha1,betha2,tpre1,eg,egpre,t,dae,sg,eegm,eegv,cn,avrovpre,avrovmm,tg,tpre2,T,holdd,cq);
        %B&C  
%       [JJ,fiter]=Jacob3(h,xx,xpre1,xpre2,alpha1,alpha2,alpha0,betha1,betha2,tpre1,eg,egpre,t,dae,sg,eegm,eegv,cn,avrovpre,avrovmm,tg,tpre2,T,holdd,cq);
        %B
%         [JJ,fiter]=Jacob2(h,xx,xpre1,xpre2,alpha1,alpha2,alpha0,betha1,betha2,tpre1,eg,egpre,t,dae,sg,eegm,cn); %Simplified Jacobian with B
       [JJ,fiter]=Jacob1(h,xx,xpre1,xpre2,alpha1,alpha2,alpha0,betha1,betha2,tpre1,eg,egpre,t,dae,sg); %Simplified Jacobian
%         fitert=fitert+fiter;
        
        [feeg,tg,cnhist,avrovmm,gvrovmm]=funcg(xx,xpre1,xpre2,eg,t,tpre1,tpre2,T,eegm,eegv,egpre,tg,h,holdd,cn,cq,avrovpre,avrovmm,gvrovpre,gvrovmm,dae); %Corrector of intermediate points (Implicit corrector)
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
        a=abs((xxx-xxxpre));
%         a=abs((xxx-xxxpre)./xxx);
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
gvrov=gvrovmm(:,end);
end

% This function computes the Jacobian
function [JJ,fiter]=Jacob1(h,xx,xpre1,xpre2,alpha1,alpha2,alpha0,betha1,betha2,tpre1,eg,egpre,t,dae,sg)
s=size(xx);
s=s(1);
delta=1e-6;
i=0;
J=zeros(s,s);
fiter=0;
while i<s
    i=i+1;
    f0=func(h,xx,xpre1,xpre2,alpha1,alpha2,alpha0,betha1,betha2,tpre1,eg,egpre,t,dae);
    fiter=fiter+1;
    yy=xx;
    yy(i)=xx(i)+delta;
    f1=func(h,yy,xpre1,xpre2,alpha1,alpha2,alpha0,betha1,betha2,tpre1,eg,egpre,t,dae);
    fiter=fiter+1;
    J(:,i)=(f1-f0)/delta;
end
sj=size(J);
sj=sj(1);
JJ=[J zeros(sj(1),sg);zeros(sg,sj(1)) eye(sg)];
end

% This function computes the Jacobian (New Jacobian)
function [J,fiter]=Jacob2(h,xx,xpre1,xpre2,alpha1,alpha2,alpha0,betha1,betha2,tpre1,eg,egpre,t,dae,sg,eegm,cn)
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
        if size(eegm)~=1    
            cc(1,cni)=nnz(eegm(:,cni)); %Number of intermediate points for every controller
            if cni==1
                b(1,cni)=nnz(eegm(:,cni)); %Number of intermediate points for both controllers (first controller)
            else
                b(1,cni)=nnz(eegm(:,cni))+b(1,cni-1); %Number of intermediate points for both controllers (adding up)
            end
        else
            cc(1,1)=1;
            cc(1,2)=0;
        end
    end
    i=0;
    eg1=eg;
    eg2=eg;
    while i<sg
        i=i+1;
        if i<cc(1,1)+1
            eg1(1,1)=eegm(i,1);
            eg2(1,1)=eegm(i,1)+delta;
        else
            eg1=eg;
            eg2=eg;
            eg1(1,2)=eegm(i,2);
            eg2(1,2)=eegm(i,2)+delta;            
        end
        f0=func(h,xx,xpre1,xpre2,alpha1,alpha2,alpha0,betha1,betha2,tpre1,eg1,egpre,t,dae);
        fiter=fiter+1;
        f1=func(h,xx,xpre1,xpre2,alpha1,alpha2,alpha0,betha1,betha2,tpre1,eg2,egpre,t,dae);
        fiter=fiter+1;
        J2(:,i)=(f1-f0)/delta;
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

function [J,fiter]=Jacob3(h,xx,xpre1,xpre2,alpha1,alpha2,alpha0,betha1,betha2,tpre1,eg,egpre,t,dae,sg,eegm,eegv,cn,avrovpre,avrovmm,tg,tpre2,T,holdd,cq)
s=size(xx);
s=s(1);
delta=1e-6;
i=0;
J=zeros(s,s);
fiter=0;
while i<s
    i=i+1;
    f0=func(h,xx,xpre1,xpre2,alpha1,alpha2,alpha0,betha1,betha2,tpre1,eg,egpre,t,dae);
    fiter=fiter+1;
    yy=xx;
    yy(i)=xx(i)+delta;
    f1=func(h,yy,xpre1,xpre2,alpha1,alpha2,alpha0,betha1,betha2,tpre1,eg,egpre,t,dae);
    fiter=fiter+1;
    J(:,i)=(f1-f0)/delta;
end
if sg~=0
    %B
    cni=0;
    while cni<cn
        cni=cni+1;
        if size(eegm)~=1    
            cc(1,cni)=nnz(eegm(:,cni)); %Number of intermediate points for every controller
            if cni==1
                b(1,cni)=nnz(eegm(:,cni)); %Number of intermediate points for both controllers (first controller)
            else
                b(1,cni)=nnz(eegm(:,cni))+b(1,cni-1); %Number of intermediate points for both controllers (adding up)
            end
        else
            cc(1,1)=1;
            cc(1,2)=0;
        end
    end
    i=0;
    eg1=eg;
    eg2=eg;
    while i<sg
        i=i+1;
        if i<cc(1,1)+1
            eg1(1,1)=eegm(i,1);
            eg2(1,1)=eegm(i,1)+delta;
        else
            eg1=eg;
            eg2=eg;
            eg1(1,2)=eegm(i,2);
            eg2(1,2)=eegm(i,2)+delta;            
        end
        f0=func(h,xx,xpre1,xpre2,alpha1,alpha2,alpha0,betha1,betha2,tpre1,eg1,egpre,t,dae);
        fiter=fiter+1;
        f1=func(h,xx,xpre1,xpre2,alpha1,alpha2,alpha0,betha1,betha2,tpre1,eg2,egpre,t,dae);
        fiter=fiter+1;
        J2(:,i)=(f1-f0)/delta;
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
    i=0;
    %C
    while i<s
        i=i+1;
        yy=xx;
        yy(i)=xx(i)+delta;
        tgsave=tg;
        avrovmmsave=avrovmm;
        [feeg,tg,cnhist,avrovmm]=funcg(xx,xpre1,xpre2,eg,t,tpre1,tpre2,T,eegm,eegv,egpre,tg,h,holdd,cn,cq,avrovpre,avrovmm,dae); %Corrector of intermediate points (Implicit corrector)
        tg=tgsave;
        avrovmm=avrovmmsave;
        f0=feeg;
        fiter=fiter+2;
        [feeg,tg,cnhist,avrovmm]=funcg(yy,xpre1,xpre2,eg,t,tpre1,tpre2,T,eegm,eegv,egpre,tg,h,holdd,cn,cq,avrovpre,avrovmm,dae); %Corrector of intermediate points (Implicit corrector)
        tg=tgsave;
        avrovmm=avrovmmsave;
        f1=feeg;
        fiter=fiter+2;
        J3(:,i)=(f1-f0)/delta;
    end
    %D
%     while cni<cn
%         cni=cni+1;
%         if size(eegm)~=1    
%             cc(1,cni)=nnz(eegm(:,cni)); %Number of intermediate points for every controller
%             if cni==1
%                 b(1,cni)=nnz(eegm(:,cni)); %Number of intermediate points for both controllers (first controller)
%             else
%                 b(1,cni)=nnz(eegm(:,cni))+b(1,cni-1); %Number of intermediate points for both controllers (adding up)
%             end
%         else
%             cc(1,1)=1;
%             cc(1,2)=0;
%         end
%     end
%     i=0;
%     eg1=eg;
%     eg2=eg;
%     while i<sg
%         i=i+1;
%         if i<cc(1,1)+1
%             eg1(1,1)=eegm(i,1);
%             eg2(1,1)=eegm(i,1)+delta;
%         else
%             eg1=eg;
%             eg2=eg;
%             eg1(1,2)=eegm(i,2);
%             eg2(1,2)=eegm(i,2)+delta;            
%         end
%         tgsave=tg;
%         avrovmmsave=avrovmm;
%         [feeg,tg,cnhist,avrovmm]=funcg(xx,xpre1,xpre2,eg1,t,tpre1,tpre2,T,eegm,eegv,egpre,tg,h,holdd,cn,cq,avrovpre,avrovmm,dae); %Corrector of intermediate points (Implicit corrector)
%         f0=feeg;
%         fiter=fiter+2;
%         tg=tgsave;
%         avrovmm=avrovmmsave;
%         [feeg,tg,cnhist,avrovmm]=funcg(xx,xpre1,xpre2,eg2,t,tpre1,tpre2,T,eegm,eegv,egpre,tg,h,holdd,cn,cq,avrovpre,avrovmm,dae); %Corrector of intermediate points (Implicit corrector)
%         tg=tgsave;
%         avrovmm=avrovmmsave;
%         f1=feeg;
%         fiter=fiter+2;
%         J4(:,i)=(f1-f0)/delta;
%     end
%     J=[J J2;J3 J4];
    J44=eye(sg);
    J=[J J2;J3 J44];
end

% si=s+sg; %increasing the size of the Jacobian by the number of intermediate points
% delta=1e-3;
% i=0;
% J=zeros(si,si);
% iter=0;
% f0=func(h,xx,xpre1,xpre2,alpha1,alpha2,alpha0,betha1,betha2,tpre1,eg,egpre,t,dae);
% [feeg,tg,cnhist,avrovmm]=funcg(xx,xpre1,xpre2,eg,t,tpre1,tpre2,T,eegm,eegv,egpre,tg,h,holdd,cn,cq,avrovpre,avrovmm,dae); %Corrector of intermediate points (Implicit corrector)
%     if sg==0
%         f0J=f0;
%     else
%         f0J=[f0;feeg']; %extending the mismatch for DC signals
%     end
%     if sg==0
%         yy=xx;
%     else
%         yy=[xx;eegv]; %extending the state vector
%     end
%     yy0=yy; %Saving yy for every increament
% while i<si
%     i=i+1;
%     yy=yy0;
%     yy(i)=yy(i)+delta; %increasing state variables by delta one by one
%     xxJ=yy(1:s,1); %splitting state vector to feed solution and DC signals separetely to functions func and funcg
%     if sg~=0
%         eegvJ=yy(s+1:end,1); %splitting state vector to feed solution and DC signals separetely to functions func and funcg
%     else
%         eegvJ=0;
%     end
%     f1=func(h,xxJ,xpre1,xpre2,alpha1,alpha2,alpha0,betha1,betha2,tpre1,eg,egpre,t,dae);
%     iter=iter+1;
%     [e1g,tg]=funcg(xxJ,xpre1,eg,t,tpre1,T,eegvJ,egpre,tg,G,u,h);
%     [feeg,tg,cnhist,avrovmm]=funcg(xx,xpre1,xpre2,eg,t,tpre1,tpre2,T,eegm,eegv,egpre,tg,h,holdd,cn,cq,avrovpre,avrovmm,dae); %Corrector of intermediate points (Implicit corrector)
%     iter=iter+2;
%     if sg==0
%         f1J=f1;
%     else
%         f1J=[f1;e1g']; %extending the mismatches for increased values by delta
%     end
%     J(:,i)=(f1J-f0J)/delta; %compute the Jacobian for column i
% end
end

function [J,fiter]=Jacob4(h,xx,xpre1,xpre2,alpha1,alpha2,alpha0,betha1,betha2,tpre1,eg,egpre,t,dae,sg,eegm,eegv,cn,avrovpre,avrovmm,tg,tpre2,T,holdd,cq)
s=size(xx);
s=s(1);
delta=1e-6;
i=0;
J=zeros(s,s);
fiter=0;
while i<s
    i=i+1;
    f0=func(h,xx,xpre1,xpre2,alpha1,alpha2,alpha0,betha1,betha2,tpre1,eg,egpre,t,dae);
    fiter=fiter+1;
    yy=xx;
    yy(i)=xx(i)+delta;
    f1=func(h,yy,xpre1,xpre2,alpha1,alpha2,alpha0,betha1,betha2,tpre1,eg,egpre,t,dae);
    fiter=fiter+1;
    J(:,i)=(f1-f0)/delta;
end
if sg~=0
    %B

    %C
    i=0;
    while i<s
        i=i+1;
        yy=xx;
        yy(i)=xx(i)+delta;
        tgsave=tg;
        avrovmmsave=avrovmm;
        [feeg,tg,cnhist,avrovmm]=funcg(xx,xpre1,xpre2,eg,t,tpre1,tpre2,T,eegm,eegv,egpre,tg,h,holdd,cn,cq,avrovpre,avrovmm,dae); %Corrector of intermediate points (Implicit corrector)
        tg=tgsave;
        avrovmm=avrovmmsave;
        f0=feeg;
        fiter=fiter+2;
        [feeg,tg,cnhist,avrovmm]=funcg(yy,xpre1,xpre2,eg,t,tpre1,tpre2,T,eegm,eegv,egpre,tg,h,holdd,cn,cq,avrovpre,avrovmm,dae); %Corrector of intermediate points (Implicit corrector)
        tg=tgsave;
        avrovmm=avrovmmsave;
        f1=feeg;
        fiter=fiter+2;
        J3(:,i)=(f1-f0)/delta;
    end
    %D
%     J=[J J2;J3 J4];
    sj=size(J);
    sj=sj(1);
    J44=eye(sg);
    J22=zeros(sj(1),sg);
    J=[J J22;J3 J44];
end
end

function [J,fiter]=Jacob5(h,xx,xpre1,xpre2,alpha1,alpha2,alpha0,betha1,betha2,tpre1,eg,egpre,t,dae,sg,eegm,eegv,cn,avrovpre,avrovmm,tg,tpre2,T,holdd,cq)
s=size(xx);
s=s(1);
delta=1e-6;
i=0;
J=zeros(s,s);
fiter=0;
while i<s
    i=i+1;
    f0=func(h,xx,xpre1,xpre2,alpha1,alpha2,alpha0,betha1,betha2,tpre1,eg,egpre,t,dae);
    fiter=fiter+1;
    yy=xx;
    yy(i)=xx(i)+delta;
    f1=func(h,yy,xpre1,xpre2,alpha1,alpha2,alpha0,betha1,betha2,tpre1,eg,egpre,t,dae);
    fiter=fiter+1;
    J(:,i)=(f1-f0)/delta;
end
if sg~=0
    %B
    cni=0;
    while cni<cn
        cni=cni+1;
        if size(eegm)~=1    
            cc(1,cni)=nnz(eegm(:,cni)); %Number of intermediate points for every controller
            if cni==1
                b(1,cni)=nnz(eegm(:,cni)); %Number of intermediate points for both controllers (first controller)
            else
                b(1,cni)=nnz(eegm(:,cni))+b(1,cni-1); %Number of intermediate points for both controllers (adding up)
            end
        else
            cc(1,1)=1;
            cc(1,2)=0;
        end
    end
    i=0;
    eg1=eg;
    eg2=eg;
    while i<sg
        i=i+1;
        if i<cc(1,1)+1
            eg1(1,1)=eegm(i,1);
            eg2(1,1)=eegm(i,1)+delta;
        else
            eg1=eg;
            eg2=eg;
            eg1(1,2)=eegm(i,2);
            eg2(1,2)=eegm(i,2)+delta;            
        end
        f0=func(h,xx,xpre1,xpre2,alpha1,alpha2,alpha0,betha1,betha2,tpre1,eg1,egpre,t,dae);
        fiter=fiter+1;
        f1=func(h,xx,xpre1,xpre2,alpha1,alpha2,alpha0,betha1,betha2,tpre1,eg2,egpre,t,dae);
        fiter=fiter+1;
        J2(:,i)=(f1-f0)/delta;
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
    i=0;
    %C
    while i<s
        i=i+1;
        yy=xx;
        yy(i)=xx(i)+delta;
        tgsave=tg;
        avrovmmsave=avrovmm;
        [feeg,tg,cnhist,avrovmm]=funcg(xx,xpre1,xpre2,eg,t,tpre1,tpre2,T,eegm,eegv,egpre,tg,h,holdd,cn,cq,avrovpre,avrovmm,dae); %Corrector of intermediate points (Implicit corrector)
        tg=tgsave;
        avrovmm=avrovmmsave;
        f0=feeg;
        fiter=fiter+2;
        [feeg,tg,cnhist,avrovmm]=funcg(yy,xpre1,xpre2,eg,t,tpre1,tpre2,T,eegm,eegv,egpre,tg,h,holdd,cn,cq,avrovpre,avrovmm,dae); %Corrector of intermediate points (Implicit corrector)
        tg=tgsave;
        avrovmm=avrovmmsave;
        f1=feeg;
        fiter=fiter+2;
        J3(:,i)=(f1-f0)/delta;
    end
    %D
    while cni<cn
        cni=cni+1;
        if size(eegm)~=1    
            cc(1,cni)=nnz(eegm(:,cni)); %Number of intermediate points for every controller
            if cni==1
                b(1,cni)=nnz(eegm(:,cni)); %Number of intermediate points for both controllers (first controller)
            else
                b(1,cni)=nnz(eegm(:,cni))+b(1,cni-1); %Number of intermediate points for both controllers (adding up)
            end
        else
            cc(1,1)=1;
            cc(1,2)=0;
        end
    end
    i=0;
    eg1=eg;
    eg2=eg;
    while i<sg
        i=i+1;
        if i<cc(1,1)+1
            eg1(1,1)=eegm(i,1);
            eg2(1,1)=eegm(i,1)+delta;
        else
            eg1=eg;
            eg2=eg;
            eg1(1,2)=eegm(i,2);
            eg2(1,2)=eegm(i,2)+delta;            
        end
        tgsave=tg;
        avrovmmsave=avrovmm;
        [feeg,tg,cnhist,avrovmm]=funcg(xx,xpre1,xpre2,eg1,t,tpre1,tpre2,T,eegm,eegv,egpre,tg,h,holdd,cn,cq,avrovpre,avrovmm,dae); %Corrector of intermediate points (Implicit corrector)
        f0=feeg;
        fiter=fiter+2;
        tg=tgsave;
        avrovmm=avrovmmsave;
        [feeg,tg,cnhist,avrovmm]=funcg(xx,xpre1,xpre2,eg2,t,tpre1,tpre2,T,eegm,eegv,egpre,tg,h,holdd,cn,cq,avrovpre,avrovmm,dae); %Corrector of intermediate points (Implicit corrector)
        tg=tgsave;
        avrovmm=avrovmmsave;
        f1=feeg;
        fiter=fiter+2;
        J4(:,i)=(f1-f0)/delta;
    end
    J=[J J2;J3 J4];
end
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
function [eeeg,tg,cnhist,avrovmm,gvrovmm]=funcg(xx,xpre1,xpre2,eg,t,tpre1,tpre2,T,eegm,eegv,egpre,tg,h,holdd,cn,cq,avrovpre,avrovmm,gvrovpre,gvrovmm,dae)
tgpre=tg;
cni=0;
mm=0;
avrc=0;
gvrc=0;
awc=0;
cnhist=zeros(1,cn);
while cni<cn
    cni=cni+1;
    if tg(1,cni)>tpre1 && tg(1,cni)<=t+0.000001
        m=0;
        while (tg(1,cni)<=t+0.000001)       
            m=m+1;
            mm=mm+1;
            hg=tg(1,cni)-tpre1; %Calculating the time step for the intermediate points
%             xg=(hg*(evaleval(xx,t,eg))+((1+hg/holdd)*xpre1)+((-(hg^2)/((hg+holdd)*holdd))*xpre2))/((2*hg+holdd)/(hg+holdd)); %Second order BDF implicit interpolator
            xg=xpre1+hg*evaleval(xpre1,tpre1,egpre)+((hg^2/h)*(xx-xpre1-(h*evaleval(xpre1,tpre1,egpre))));
            xg=quz(xg,cq);
            [f2,gvrov]=contim(eegm,egpre,T,xg,m,mm,cni,avrovpre,avrovmm,gvrovpre,gvrovmm); %This function contains the controllers' equations
            eeeg(1,mm)=f2;
            temp(1,cni)=tg(1,cni);
            tg(1,cni)=tg(1,cni)+T(1,cni);
            cnhist(mm,cni)=1;
            if cni==1
                gvrovmm(:,m)=gvrov;
                gvrc=gvrc+1;
            elseif cni==2
               %avrovmm(:,m)=avrov;
               vrovmm(:,m)=avrovmm(:,m);
               avrc=avrc+1;
            end
        end
    end
end
if mm==0
        eeeg=0;
        maxl=0;
        cnhist=cnhist;
        awc=0;
end
tg=tgpre; %reset the counter of the intermediate points to be counted again in the next iterations
end

%This function computes the mismatch for intermediate points regarding the
%controllers' signals (cni is used to varry the controller equation)
function [f2,gvrov]=contim(eegm,egpre,T,xg,m,mm,cni,avrovpre,avrovmm,gvrovpre,gvrov)
if cni==1
            if m==1
                [aw,gvrov]=torque(egpre(1,cni),T(1,cni),xg,gvrovpre);
                aw=max(aw,0); %upper limnit anti wind up
%                 aw=min(aw,1); %lower limnit anti wind up
                aw=min(aw,1); %lower limnit anti wind up
                f2=eegm(mm,cni)-aw;
                avrov=0;
            else
                [aw,gvrov]=torque(eegm(mm-1,cni),T(1,cni),xg,gvrov(:,m-1));
                aw=max(aw,0); %upper limnit anti wind up
%                 aw=min(aw,1); %lower limnit anti wind up
                aw=min(aw,1); %lower limnit anti wind up
                f2=eegm(mm,cni)-aw;
                avrov=0;
            end
elseif cni==2
            global fuzzyAVR3Bus;
            Vter=sqrt(xg(11,1)^2+xg(12,1)^2);
            if m==1
                aw=evalfis(fuzzyAVR3Bus,Vter);
                %[aw,avrov]=avr(avrovpre,T(1,cni),xg); %Controller's equation
                aw=max(aw,0); %upper limnit anti wind up
%                 aw=min(aw,0.01); %lower limnit anti wind up
                aw=min(aw,0.01); %lower limnit anti wind up
                f2=eegm(mm,cni)-aw;
                gvrov=0;
            else
                aw=evalfis(fuzzyAVR3Bus,Vter);
                %[aw,avrov]=avr(avrovmm(:,m-1),T(1,cni),xg); %Controller's equation
                aw=max(aw,0); %upper limnit anti wind up
%                 aw=min(aw,0.01); %lower limnit anti wind up
                aw=min(aw,0.01); %lower limnit anti wind up
                f2=eegm(mm,cni)-aw; 
                gvrov=0;
            end
end
end

%Explicit interpolator
function [ee0gm,ee0gv,eg,tg,sg,avrovmm,gvrovmm]=predictorg(xx,xpre1,xpre2,t,tpre1,tpre2,holdd,T,eg,egpre,tg,cp,cn,cq,avrov,gvrov,dae)
tgpre=tg;  
cni=0;
mm=0;
avrc=0;
gvrc=0;
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
            
            [eg(1,cni),avrov,gvrov]=contex(eg,T,xg,cni,avrov,gvrov); %This function contains the controllers' equations
            ee0gm(mm,cni)=eg(1,cni);
            ee0gv(mm,1)=eg(1,cni);
            temp(1,cni)=tg(1,cni);
            tg(1,cni)=tg(1,cni)+T(1,cni);
            if cni==1
                gvrc=gvrc+1;
                gvrovmm(:,m)=gvrov;
            elseif cni==2
                avrc=avrc+1;
                avrovmm(:,m)=avrov;
            end
        end     
    end
end
if avrc==0
    avrovmm=avrov;
end
if gvrc==0
    gvrovmm=gvrov;
end
if mm~=0
    ee0gm=quz(ee0gm,cq);
    ee0gv=quz(ee0gv,cq);
    sg=mm;
else
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
function [egtemp,avrov,gvrov]=contex(eg,T,xg,cni,avrov,gvrov)
if cni==1
    [egtemp,gvrov]=torque(eg(1,cni),T(1,cni),xg,gvrov);
    egtemp=max(egtemp,0); %upper limnit anti wind up
%     egtemp=min(egtemp,1); %lower limnit anti wind up
    egtemp=min(egtemp,1); %lower limnit anti wind up
elseif cni==2
    global fuzzyAVR3Bus;
    Vter=sqrt(xg(11,1)^2+xg(12,1)^2);
    egtemp=evalfis(fuzzyAVR3Bus,Vter);
%     [egtemp,avrov]=avr(avrov,T(1,cni),xg); %Controller's equation
    egtemp=max(egtemp,0); %upper limnit anti wind up
%     egtemp=min(egtemp,0.01); %lower limnit anti wind up
    egtemp=min(egtemp,0.01); %lower limnit anti wind up

end
end

function [et2,gvrov]=torque(et2,T,xx,gvrov)
%Data
kp=12;
ki=0.2;
Droop=1;
    zo=0.606781005859375;
    wS=1+xx(3,1);
    Te=xx(17,1)*xx(14,1)-xx(18,1)*xx(13,1);
%initialization (old values)
T2=T;
x5=gvrov;
%equations
x3=zo-wS+1-(zo*(Droop));
% x3=zo-wS+1-(zo*(Droop));
% x3=zo-wS+1;
x4=x3*kp;
x5=x5+T2*x3*ki;
x6=x5+x4;
%assigning old values for the next sample
gvrov=x5;
et2=x6;
end

function [et1,avrov]=avr(avrov,T1,xx)
%Data
Tmm=0.2; %0.2
Te=0.12; %0.12
vo=1;
kp=0.003; %0.01
ki=0.001; %0.01
Ge=1;
kf=0.8; %0.1
Tf=0.9; %0.9
Vter=sqrt(xx(11,1)^2+xx(12,1)^2);
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
x5=x4+x3;
x7=(kf/T1)*(x62-x61);
x8=(Tf*x8-T1*x8+T1*x7)/Tf;
et1=x5;
%assigning old values for the next sample
avrov=[x1 x2 x3 x4 x5 x62 x7 x8]';
end

% function et2=torque(et2,T,xg)
% omega=xg(3,1);
% %Data
% Tsm=0.9;
% sigma=0.04;
% zo=0.606781005859375;
% wS=1+omega;
% %initialization (old values)
% T2=T;
% x9=et2;
% %equations
% % x9=(x9*(Tsm*sigma*wN-sigma*wN*T2)+T2*(zo*sigma*wN-wS+wN))/(Tsm*sigma*wN);
% x9=(x9*(Tsm*sigma-sigma*T2)+T2*(zo*sigma-wS+1))/(Tsm*sigma);
% %assigning old values for the next sample
% et2=x9;
% end
% 
% function [et1,avrov]=avr(avrov,T1,xg)
% %Data
% Tmm=0.2;
% Te=0.12;
% vo=1;
% kp=0.003; %0.01
% ki=0.001; %0.01
% Ge=1;
% kf=0.8; %0.1 %500
% Tf=0.9; %0.9
% Vter=sqrt(xg(9,1)^2+xg(10,1)^2);
% %initialization (old values)
% x1=avrov(1,1);
% x2=avrov(2,1);
% x3=avrov(3,1);
% x4=avrov(4,1);
% x5=avrov(5,1);
% x61=avrov(6,1);
% x7=avrov(7,1);
% x8=avrov(8,1);
% %equations
% x1=((Tmm-T1)*x1+T1*Vter)/Tmm;
% x2=vo-x1-x8;
% x3=kp*x2;
% x4=x4+T1*ki*x2;
% x5=x4+x3;
% x62=(x61*(Te-T1)+T1*Ge*x5)/Te;
% x7=(kf/T1)*(x62-x61);
% % x7=(kf/T1)*(x62);
% x8=(Tf*x8-T1*x8+T1*x7)/Tf;
% et1=x62;
% %assigning old values for the next sample
% avrov=[x1 x2 x3 x4 x5 x62 x7 x8]';
% end

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
j=sqrt(-1);
Ge=1;
Tee=0.12;
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
x51=et(1,2); %0.001196163166093
% Vf=0.0010873;
x61=et(1,1); %0.709060632610000
% Tm=0.709060632610000;
%Assiginng variables
psif=xx(1,1);
psiq1=xx(2,1);
omega=xx(3,1);
deltas=xx(4,1);
Vf=xx(5,1);
Tm=xx(6,1);
ifd=xx(7,1);
iq1=xx(8,1);
vd=xx(9,1);
vq=xx(10,1);
vx1=xx(11,1);
vy1=xx(12,1);
id=xx(13,1);
iq=xx(14,1);
ix1=xx(15,1);
iy1=xx(16,1);
psid=xx(17,1);
psiq=xx(18,1);
vx2=xx(19,1);
vy2=xx(20,1);
ix2=xx(21,1);
iy2=xx(22,1);
vx3=xx(23,1);
vy3=xx(24,1);
ix3=xx(25,1);
iy3=xx(26,1);
ix32=xx(27,1);
iy32=xx(28,1);
%Extra
if t>=3 && t<=3.1
    B=-0.75; 
    G=0;
end
% if t>3 && t<6
%     Ppv=0.6;
% end
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

%FieldVoltage and Torque of Generators
f5=(x51*Ge-Vf)/Tee;
f6=(x61-Tm)/0.3;

%%%Algebraic equations
%flux-current
f7=Ldd*id+Ldf*ifd-psid; 
f8=Lqq*iq+Lqq1*iq1-psiq; 
f9=Lff*ifd+Ldf*id-psif; 
f10=Lq1q1*iq1+Lqq1*iq-psiq1;
%stator park
f11=vd+Ra*id+psiq; 
f12=vq+Ra*iq-psid;
%network line1
f13=vx1-vx3-Re1*ix1+Xe1*iy1; 
f14=vy1-vy3-Re1*iy1-Xe1*ix1;
%dqxy
f15=vd-vx1*cos(theta0)-vy1*sin(theta0); 
f16=vq-vx1*sin(theta0)+vy1*cos(theta0); 
f17=id-ix1*cos(theta0)-iy1*sin(theta0); 
f18=iq-ix1*sin(theta0)+iy1*cos(theta0);
%networkPV PV
f19=-Ppv+vx2*ix2+vy2*iy2;
f20=-Qpv-vx2*iy2+vy2*ix2;
%networkpv line2
f21=vx2-vx3-Re2*ix2+Xe2*iy2;
f22=vy2-vy3-Re2*iy2-Xe2*ix2;
%network load
f23=-Pload+vx3*ix3+vy3*iy3;
f24=-Qload-vx3*iy3+vy3*ix3;
f25=-ix32+G*vx3-B*vy3;
f26=-iy32+B*vx3+G*vy3;
%network KCL
f27=ix1+ix2-ix3-ix32;
f28=iy1+iy2-iy3-iy32;
%sum
f=[f1 f2 f3 f4 f5 f6 f7 f8 f9 f10 f11 f12 f13 f14 f15 f16 f17 f18 f19 f20 f21 f22 f23 f24 f25 f26 f27 f28]';
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
ix3=-real(phasIload);
iy3=-imag(phasIload);
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
ssi=[0.994883515646896;0.670058836367880;6.50365295768179e-06;0.805918434474973;0.000991821289062500;0.606781005859375;0.905548718909198;1.36398958407194e-05;-0.728250367707204;0.685363244088517;1.00003568650393;-1.14814909617123e-05;-0.543883810960619;0.304558241246316;0.604793293296409;-0.150964878431310;0.686886035294749;0.730969786762007;0.996795749182842;0.000390119303810904;0.300983950783438;-0.0500429304514542;0.993685823814105;0.000288580706758570;0.905777244079846;-0.201007808882764;0;0];
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