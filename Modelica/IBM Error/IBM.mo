within ;
package IBM

  block IBMController
    Modelica.Blocks.Interfaces.RealOutput y annotation (Placement(transformation(
            extent={{60,-20},{98,18}}), iconTransformation(extent={{60,-20},{98,18}})));
    Modelica.Blocks.Interfaces.RealInput uf annotation (Placement(transformation(
            extent={{-100,-64},{-60,-24}}), iconTransformation(extent={{-100,-64},
              {-60,-24}})));
    Modelica.Blocks.Interfaces.RealInput uo annotation (Placement(transformation(
            extent={{-100,-22},{-60,18}}),iconTransformation(extent={{-100,-22},{-60,
              18}})));
    Modelica.Blocks.Interfaces.RealInput ud annotation (Placement(transformation(
            extent={{-100,-18},{-60,22}}),iconTransformation(extent={{-100,20},{-60,
              60}})));
    parameter Real Ki=0.07;
    parameter Real T=0.01;
    parameter Real t1=4*T,t2=3*T,t3=2*T,t4=T,t5=0;
    Real e1(start=0),e2,e3,e4,e5(start=0);
    Real y1,y2,y3,y4,y5;
    Real zf=uf;
    Real z2=uo;
    Real z1;
    Clock c1=Clock(1,20);
  equation
    z1=der(z2);
  algorithm
     when Clock(c1,"ImplicitTrapezoid") then
     //y1:=t1*pre(z1)+pre(z2);
     y1:=t1*pre(z1)+pre(z2)+(t1^2/0.05^2)*(z2-pre(z2)-0.05*pre(z1));
     e1:=e5 + Ki*T*(zf-y1);

    //y2:=t2*pre(z1)+pre(z2);
    y2:=t2*pre(z1)+pre(z2)+(t2^2/0.05^2)*(z2-pre(z2)-0.05*pre(z1));
    e2:=e1 + Ki*T*(zf-y2);

    //y3:=t3*pre(z1)+pre(z2);
    y3:=t3*pre(z1)+pre(z2)+(t3^2/0.05^2)*(z2-pre(z2)-0.05*pre(z1));
    e3:=e2 + Ki*T*(zf-y3);

    //y4:=t4*pre(z1)+pre(z2);
    y4:=t4*pre(z1)+pre(z2)+(t4^2/0.05^2)*(z2-pre(z2)-0.05*pre(z1));
    e4:=e3 + Ki*T*(zf-y4);

    //y5:=t5*pre(z1)+pre(z2);
    y5:=t5*pre(z1)+pre(z2)+(t5^2/0.05^2)*(z2-pre(z2)-0.05*pre(z1));
    e5:=e4 + Ki*T*(zf-y5);
     end when;
  equation
    y=e5;
    annotation (Icon(coordinateSystem(preserveAspectRatio=false), graphics={
          Rectangle(
            extent={{-60,60},{60,-60}},
            lineColor={28,108,200},
            fillColor={127,0,0},
            fillPattern=FillPattern.Solid),
          Text(
            extent={{-50,50},{54,-50}},
            textColor={28,108,200},
            textString="IBM
Digital Controller"),
          Text(
            extent={{-32,76},{34,66}},
            textColor={28,108,200},
            textString="%name")}), Diagram(coordinateSystem(preserveAspectRatio=false)),
      experiment(__Dymola_fixedstepsize=0.5, __Dymola_Algorithm="Rkfix2"));
  end IBMController;
  annotation (uses(Modelica(version="4.0.0")));
end IBM;
