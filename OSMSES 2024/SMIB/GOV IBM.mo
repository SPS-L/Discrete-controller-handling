within ;
model GOVIBM
  Modelica.Blocks.Interfaces.RealInput u0
    annotation (Placement(transformation(extent={{-140,20},{-100,60}})));
  Modelica.Blocks.Interfaces.RealInput uf
    annotation (Placement(transformation(extent={{-140,-60},{-100,-20}})));
  Modelica.Blocks.Interfaces.RealOutput y
    annotation (Placement(transformation(extent={{100,-10},{120,10}})));
  parameter Real sigma=0.02;
  parameter Real Tsm=0.1;
  parameter Real T=0.001;
  parameter Real t1=T,t2=2*T,t3=3*T,t4=4*T,t5=5*T;
  //Real e1,e2,e3,e4,e5(start=20);
  Real z0=u0;
  Real zf=uf;
  //Real zd(start=0)=ud;
  Real z1;
  Real hg;
  //Real y1,y2,y3,y4,y5;
  Real e(start=0.01);
  Real w(start=1);
  Real tt=time;
  Real prett;
  Real ttdiff;
  Real ttfloor;
  Real prettfloor;
  Real ttfloordiff;
  Integer i;
  Real earray[1,17];
  Clock c1=Clock(1,80);
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
     else
       z1:=(zf - pre(zf))/ttdiff;
     end if;
     while i < ttfloordiff loop
       i:=i + 1;
       hg:=prettfloor*T + i*T - prett;
       w:=(prettfloor*T+i*T-prett)*pre(z1)+pre(zf);
       e:=(e*((Tsm*sigma) - (sigma*T)) + T*(z0*sigma - w + 1))/(Tsm*sigma);
       earray[1,i]:=e;
       if i==ttfloordiff then
         earray[1,7]:=earray[1,i];
       end if;
     end while;



   //y1:=t1*pre(zd)+pre(zf);
   //e1:=(e5*((Tsm*sigma) - (sigma*T)) + T*(z0*sigma - zf + 1))/(Tsm*sigma);

  //y2:=t2*pre(zd)+pre(zf);
   //e2:=(e1*((Tsm*sigma) - (sigma*T)) + T*(z0*sigma - zf + 1))/(Tsm*sigma);

  //y3:=t3*pre(zd)+pre(zf);
   //e3:=(e2*((Tsm*sigma) - (sigma*T)) + T*(z0*sigma - zf + 1))/(Tsm*sigma);

  //y4:=t4*pre(zd)+pre(zf);
   //e4:=(e3*((Tsm*sigma) - (sigma*T)) + T*(z0*sigma - zf + 1))/(Tsm*sigma);

  //y5:=t5*pre(zd)+pre(zf);
   //e5:=(e4*((Tsm*sigma) - (sigma*T)) + T*(z0*sigma - zf + 1))/(Tsm*sigma);
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
          extent={{-74,40},{74,-44}},
          textColor={244,125,35},
          textString="IBM
GOV")}));
end GOVIBM;
