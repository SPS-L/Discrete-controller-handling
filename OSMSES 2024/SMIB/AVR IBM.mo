within ;
model AVRIBM
  Modelica.Blocks.Interfaces.RealInput uter
    annotation (Placement(transformation(extent={{-140,20},{-100,60}})));
  Modelica.Blocks.Interfaces.RealOutput y
    annotation (Placement(transformation(extent={{100,-10},{120,10}})));
  parameter Real Kp=11;
  parameter Real Ki=3;
  parameter Real Tsm=0.3;
  parameter Real Tf=0.9;
  parameter Real Kf=0;
  parameter Real T=0.001;
  parameter Real t1=T,t2=2*T,t3=3*T,t4=4*T,t5=5*T;
  Real x1(start=1),x2(start=0),x3(start=0),x4(start=1),x5(start=1),x7(start=0),x8(start=0);
  //Real x11(start=1),x12(start=1),x13(start=1),x14(start=1),x15(start=1);
  //Real x21(start=0),x22(start=0),x23(start=0),x24(start=0),x25(start=0);
  //Real x31(start=0),x32(start=0),x33(start=0),x34(start=0),x35(start=0);
  //Real x41(start=1),x42(start=1),x43(start=1),x44(start=1),x45(start=1);
  //Real x51(start=1),x52(start=1),x53(start=1),x54(start=1),x55(start=1);
  //Real x71(start=0),x72(start=0),x73(start=0),x74(start=0),x75(start=0);
  //Real x81(start=0),x82(start=0),x83(start=0),x84(start=0),x85(start=0);
  //Real w1,w2,w3,w4,w5;
  //Real w10,w20,w30,w40,w50;
    //The second 5
  Real zf=uter;
  //Real zd=ud;
  Clock c1=Clock(1,80);
  Real z1,z2;
  Real w1;
  Real w2;
  Real tt=time;
  Real prett;
  Real ttdiff;
  Real ttfloor;
  Real prettfloor;
  Real ttfloordiff;
  Integer i;
  Real earray[1,17];
  //Real zd0=ud0;
  Modelica.Blocks.Interfaces.RealInput uf0
    annotation (Placement(transformation(extent={{-140,-60},{-100,-20}})));
    Real z0=uf0;
    //Real z0d=uf0d;
algorithm
   when Clock(c1,"ImplicitTrapezoid") then
       i:=0;
       prett:=pre(tt);
       ttdiff:=tt - prett;
       ttfloor:=floor(tt/T);
       prettfloor:=floor(prett/T);
       ttfloordiff:=ttfloor - prettfloor;
     if ttdiff==0 then
       z1:=(zf - pre(zf))/0.005;
       z2:=(z0 - pre(z0))/0.005;
     else
       z1:=(zf - pre(zf))/ttdiff;
       z2:=(z0 - pre(z0))/ttdiff;
     end if;
       while i < ttfloordiff loop
           i:=i + 1;
           w1:=(prettfloor*T+i*T-prett)*pre(z1)+pre(zf);
           w2:=(prettfloor*T+i*T-prett)*pre(z2)+pre(z0);
           x1:=(((Tsm - T)*x1) + (T*w1))/(Tsm);
           x2:=1 - x1 - x8;
           x3:=Kp*x2;
           x4:=x4 + T*Ki*x2;
           x5:=x3 + x4;
           x7:=(Kf/T)*(w2 - pre(w2));
           x8:=(Tf*x8 - T*x8 + T*x7)/Tf;
           earray[1,i]:=x5;
           if i==ttfloordiff then
             earray[1,7]:=earray[1,i];
           end if;
       end while;


     //w1:=t1*pre(zd)+pre(zf);
     //w10:=t1*pre(z0d)+pre(z0);
     //x11:=(((Tsm - T)*x15) + (T*w1))/(Tsm);
     //x21:=1 - x11 - x85;
     //x31:=Kp*x21;
     //x41:=x45 + T*Ki*x21;
     //x51:=x31 + x41;
     //x71:=(Kf/T)*(w10 - pre(w10));
     //x81:=(Tf*x85 - T*x85 + T*x71)/Tf;

     //w2:=t2*pre(zd)+pre(zf);
     //w20:=t2*pre(z0d)+pre(z0);
     //x12:=(((Tsm - T)*x11) + (T*w2))/(Tsm);
     //x22:=1 - x12 - x81;
     //x32:=Kp*x22;
     //x42:=x41 + T*Ki*x22;
     //x52:=x32 + x42;
     //x72:=(Kf/T)*(w20 - pre(w20));
     //x82:=(Tf*x81 - T*x81 + T*x72)/Tf;

     //w3:=t3*pre(zd)+pre(zf);
     //w30:=t3*pre(z0d)+pre(z0);
     //x13:=(((Tsm - T)*x12) + (T*w3))/(Tsm);
     //x23:=1 - x13 - x82;
     //x33:=Kp*x23;
     //x43:=x42 + T*Ki*x23;
     //x53:=x33 + x43;
     //x73:=(Kf/T)*(w30 - pre(w30));
     //x83:=(Tf*x82 - T*x82 + T*x73)/Tf;

     //w4:=t4*pre(zd)+pre(zf);
     //w40:=t4*pre(z0d)+pre(z0);
     //x14:=(((Tsm - T)*x13) + (T*w4))/(Tsm);
     //x24:=1 - x14 - x83;
     //x34:=Kp*x24;
     //x44:=x43 + T*Ki*x24;
     //x54:=x34 + x44;
     //x74:=(Kf/T)*(w40 - pre(w40));
     //x84:=(Tf*x83 - T*x83 + T*x74)/Tf;

     //w5:=t5*pre(zd)+pre(zf);
     //w50:=t5*pre(z0d)+pre(z0);
     //x15:=(((Tsm - T)*x14) + (T*w5))/(Tsm);
     //x25:=1 - x15 - x84;
     //x35:=Kp*x25;
     //x45:=x44 + T*Ki*x25;
     //x55:=x35 + x45;
     //x75:=(Kf/T)*(w50 - pre(w50));
     //x85:=(Tf*x84 - T*x84 + T*x75)/Tf;



   //w10:=t1*pre(zd0) + pre(z0);
   //x11:=1-((((Tsm - T)*x51) + (T*w1))/(Tsm));
   //x11:=1-w1-y10;
   //x12:=x52+(T*Ki*(1-w1));
   //e1:=x12;
   //e1:=x11*Kp + x12;
   //x17:=(Kf/T)*(w10 - w50);
   //x18:=(Tf*x58 - T*x58 + T*x17)/Tf;

   //w2:=t2*pre(zd)+pre(zf);
   //w20:=t2*pre(zd0) + pre(z0);
   //x21:=1-((((Tsm - T)*x11) + (T*w2))/(Tsm));
   //x21:=1-w2-y20;
   //x22:=x12+(T*Ki*(1-w2));
   //e2:=x22;
   //e2:=x21*Kp + x22;
   //x27:=(Kf/T)*(w20 - w10);
   //x28:=(Tf*x18 - T*x18 + T*x27)/Tf;

   //w3:=t3*pre(zd)+pre(zf);
   //w30:=t3*pre(zd0) + pre(z0);
   //x31:=1-((((Tsm - T)*x21) + (T*w3))/(Tsm));
   //x31:=1-w3-y30;
   //x32:=x22+(T*Ki*(1-w3));
   //e3:=x32;
   //e3:=x31*Kp + x32;
   //x37:=(Kf/T)*(w30 - w20);
   //x38:=(Tf*x28 - T*x28 + T*x37)/Tf;

   //w4:=t4*pre(zd)+pre(zf);
   //w40:=t4*pre(zd0) + pre(z0);
   //x41:=1-((((Tsm - T)*x31) + (T*w4))/(Tsm));
   //x41:=1-w4-y40;
   //x42:=x32+(T*Ki*(1-w4));
   //e4:=x42;
   //e4:=x41*Kp + x42;
   //x47:=(Kf/T)*(w40 - w30);
   //x48:=(Tf*x38 - T*x38 + T*x47)/Tf;

   //w5:=t5*pre(zd)+pre(zf);
   //w50:=t5*pre(zd0) + pre(z0);
   //x51:=1-((((Tsm - T)*x41) + (T*w5))/(Tsm));
   //x51:=1-w5-y50;
   //x52:=x42+(T*Ki*(1-w5));
   //e5:=x52;
   //e5:=x51*Kp + x52;
   //x57:=(Kf/T)*(w50 - w40);
   //x58:=(Tf*x48 - T*x48 + T*x57)/Tf;
  end when;
equation
  y=earray[1,7];
  annotation (uses(Modelica(version="4.0.0")), Icon(graphics={
       Rectangle(
          extent={{-100,100},{100,-100}},
          lineColor={28,108,200},
          lineThickness=1),
        Rectangle(
          extent={{-100,100},{100,-100}},
          lineColor={0,127,0},
          fillColor={0,127,127},
          fillPattern=FillPattern.Solid),
        Text(
          extent={{-58,32},{60,-34}},
          textColor={244,125,35},
          textString="IBM
AVR")}),
    experiment(
      StopTime=30,
      __Dymola_fixedstepsize=0.05,
      __Dymola_Algorithm="Rkfix2"));
end AVRIBM;
