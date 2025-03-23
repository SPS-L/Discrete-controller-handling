within ;
package IBM

  block IBMController
    Modelica.Blocks.Interfaces.RealOutput y annotation (Placement(transformation(
            extent={{142,-18},{180,20}}),
                                        iconTransformation(extent={{142,-18},{
              180,20}})));
    Modelica.Blocks.Interfaces.RealInput uf annotation (Placement(transformation(
            extent={{-180,-60},{-140,-20}}),iconTransformation(extent={{-180,-60},
              {-140,-20}})));
    Modelica.Blocks.Interfaces.RealInput uo annotation (Placement(transformation(
            extent={{-180,20},{-140,60}}),iconTransformation(extent={{-180,20},{-140,
              60}})));
    parameter Real Ki=0.07;
    parameter Real T=0.01;
    parameter Real t1=T,t2=2*T,t3=3*T,t4=4*T,t5=5*T;
    //Real e1(start=0),e2,e3,e4,e5(start=0);
    Real e(start=0);
    //Real y1,y2,y3,y4,y5,y6;
    Real w;
    Real zf=uf;
    Real z1;
    Real z2=uo;
    Real tt=time;
    Real prett;
    Real ttdiff;
    Real ttfloor;
    Real prettfloor;
    Real ttfloordiff;
    Integer i;
    Real earray[1,7];
  //Clock c1=Clock(tt);
    Clock c1=Clock(1,20);
  algorithm
     when Clock(c1,"ImplicitTrapezoid") then
     //when a then
     //when tt>prett then
     //when sample(time) then
     i:=0;
     prett:=pre(tt);
     ttdiff:=tt - prett;
     ttfloor:=floor(tt/T);
     prettfloor:=floor(prett/T);
     ttfloordiff:=ttfloor - prettfloor;
     if ttdiff==0 then
         z1:=zf - pre(zf)/0.05;
     else
       z1:=(zf - pre(zf))/(tt-prett);
     end if;
     while i < ttfloordiff loop
       i:=i + 1;
       //w:=(prettfloor*T+i*T-prett)*pre(z1)+pre(z2);
       //w:=(prettfloor*T+i*T-prett)*z3+pre(z2);
       w:=(prettfloor*T+i*T-prett)*pre(z1)+pre(z2)+((prettfloor*T+i*T-prett)^2/(tt-prett)^2)*(z2-pre(z2)-(tt-prett)*pre(z1));
       e:=e + Ki*T*(zf - w);
       earray[1,i]:=e;
       if i==ttfloordiff then
         earray[1,7]:=earray[1,i];
       end if;
     end while;
     //y1:=t1*pre(z1)+pre(z2);
     //y1:=t1*pre(z1)+pre(z2)+(t1^2/(tt-prett)^2)*(z2-pre(z2)-(tt-prett)*pre(z1));
     //e1:=e5 + Ki*T*(zf-y1);

    //y2:=t2*pre(z1)+pre(z2);
    //y2:=t2*pre(z1)+pre(z2)+(t2^2/(tt-prett)^2)*(z2-pre(z2)-(tt-prett)*pre(z1));
    //e2:=e1 + Ki*T*(zf-y2);

    //y3:=t3*pre(z1)+pre(z2);
    //y3:=t3*pre(z1)+pre(z2)+(t3^2/(tt-prett)^2)*(z2-pre(z2)-(tt-prett)*pre(z1));
    //e3:=e2 + Ki*T*(zf-y3);

    //y4:=t4*pre(z1)+pre(z2);
    //y4:=t4*pre(z1)+pre(z2)+(t4^2/(tt-prett)^2)*(z2-pre(z2)-(tt-prett)*pre(z1));
    //e4:=e3 + Ki*T*(zf-y4);

    //y5:=t5*pre(z1)+pre(z2);
    //y5:=t5*pre(z1)+pre(z2)+(t5^2/(tt-prett)^2)*(z2-pre(z2)-(tt-prett)*pre(z1));
    //e5:=e4 + Ki*T*(zf-y5);
     end when;
  equation
    //y=e;
    y=earray[1,7];
    annotation (Icon(coordinateSystem(preserveAspectRatio=false), graphics={
          Rectangle(
            extent={{-140,100},{142,-100}},
            lineColor={28,108,200},
            fillColor={0,127,127},
            fillPattern=FillPattern.Solid),
          Text(
            extent={{-54,50},{50,-50}},
            textColor={244,125,35},
            textString="IBM
Integral
Controller")}),                    Diagram(coordinateSystem(preserveAspectRatio=false)),
      experiment(__Dymola_fixedstepsize=0.5, __Dymola_Algorithm="Rkfix2"));
  end IBMController;
  annotation (uses(Modelica(version="4.0.0")));
end IBM;
