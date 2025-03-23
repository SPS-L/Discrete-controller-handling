within ;
model AVRCode

  Modelica.Blocks.Interfaces.RealInput Vter
    annotation (Placement(transformation(extent={{-180,20},{-140,60}}),
        iconTransformation(extent={{-180,20},{-140,60}})));
  Modelica.Blocks.Interfaces.RealInput Vf
    annotation (Placement(transformation(extent={{-180,-60},{-140,-20}}),
        iconTransformation(extent={{-180,-60},{-140,-20}})));
  Modelica.Blocks.Interfaces.RealOutput y
    annotation (Placement(transformation(extent={{140,-10},{160,10}}),
        iconTransformation(extent={{140,-10},{160,10}})));
  parameter Real Kp=11;
  parameter Real Ki=3;
  parameter Real Tsm=0.3;
  parameter Real Kf=0;
  parameter Real Tf=0.9;
  parameter Real T=0.001;
  Real x1(start=1),x2(start=0),x3(start=0),x4(start=1),x5(start=1),x7(start=0),x8(start=0);
  Real zf=Vter;
  Real z0=Vf;
  Clock c1=Clock(1,1000);
algorithm
  when Clock(c1, "ImplicitTrapezoid") then
  x1:=(((Tsm - T)*x1) + (T*zf))/(Tsm);
  x2:=1 - x1 - x8;
  x3:=Kp*x2;
  x4:=x4 + T*Ki*x2;
  x5:=x3 + x4;
  x7:=(Kf/T)*(z0 - pre(z0));
  x8:=(Tf*x8 - T*x8 + T*x7)/Tf;
  end when;
equation
  y=x5;
  annotation (uses(Modelica(version="4.0.0")), Icon(graphics={Text(
          extent={{-46,40},{48,-38}},
          textColor={28,108,200},
          textString="AVR"), Rectangle(
          extent={{-140,100},{140,-100}},
          lineColor={28,108,200},
          lineThickness=1)}));
end AVRCode;
