within ;
model GOVCODE

  Modelica.Blocks.Interfaces.RealInput u0
    annotation (Placement(transformation(extent={{-180,20},{-140,60}}),
        iconTransformation(extent={{-180,20},{-140,60}})));
  Modelica.Blocks.Interfaces.RealInput w
    annotation (Placement(transformation(extent={{-180,-60},{-140,-20}}),
        iconTransformation(extent={{-180,-60},{-140,-20}})));
  Modelica.Blocks.Interfaces.RealOutput Tm
    annotation (Placement(transformation(extent={{140,-8},{160,12}}),
        iconTransformation(extent={{140,-8},{160,12}})));

  parameter Real sigma=0.02;
  parameter Real Tsm=0.1;
  parameter Real T=0.001;
  //Real e(start=0.1);
  Real e;
  Real z0=u0;
  Real zf=w;
  Clock c2=Clock(1,1000);
algorithm
  when Clock(c2, "ImplicitTrapezoid") then
e:=(e*((Tsm*sigma) - (sigma*T)) + T*(z0*sigma - zf + 1))/(Tsm*sigma);
  end when;
equation
  Tm=e;
  annotation (uses(Modelica(version="4.0.0")), Icon(graphics={Text(
          extent={{-52,38},{60,-36}},
          textColor={28,108,200},
          textString="GOV"), Rectangle(
          extent={{-140,100},{140,-100}},
          lineColor={28,108,200},
          lineThickness=1)}));
end GOVCODE;
