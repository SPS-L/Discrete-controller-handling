within ;
model SMIBAVRCodeGOVCode4

  OpenIPSL.Electrical.Buses.InfiniteBus infinite_bus(
    P_0=-1998000000,
    Q_0=87066000,
    angle_0=0,
    displayPF=false,
    v_0=0.90081)                                                                                                                       annotation (
    Placement(visible = true, transformation(origin={178,12},   extent = {{-10, 10}, {10, -10}}, rotation = 180)));
  OpenIPSL.Electrical.Branches.PwLine line_1(
    B=0,
    G=0,
    R=0,
    S_b=100000000,
    X=0.5*100/2220,
    displayPF=false)                                                                                                    annotation (
    Placement(visible = true, transformation(origin={94,12},   extent = {{20, 10}, {40, 30}}, rotation = 0)));
  OpenIPSL.Electrical.Branches.PwLine line_2(
    B=0,
    G=0,
    R=0,
    S_b=100000000,
    X=0.93*100/2220)                                                                                  annotation (
    Placement(visible = true, transformation(origin={94,12},   extent = {{20, -30}, {40, -10}}, rotation = 0)));
  OpenIPSL.Electrical.Events.PwFault fault(
    R=0,
    X=0.01*100/2220,
    t1=2,
    t2=2.1)                                                                               annotation (
    Placement(visible = true, transformation(origin={94,14},   extent = {{20, -60}, {40, -40}}, rotation = 0)));
  OpenIPSL.Electrical.Buses.Bus B3(displayPF=false, v_0=0.900810000000000)     annotation (
    Placement(visible = true, transformation(origin={86,12},   extent = {{60, -10}, {80, 10}}, rotation = 0)));
  OpenIPSL.Electrical.Buses.Bus B2(
    angle_0=0.351222533717268,
    displayPF=false,
    v_0=0.944299492912195)                                                                                  annotation (
    Placement(visible = true, transformation(origin={88,12},   extent = {{-20, -10}, {0, 10}}, rotation = 0)));
  OpenIPSL.Electrical.Buses.Bus B1(angle_0=0.494677176989154, displayPF=false)     annotation (
    Placement(visible = true, transformation(origin={96,12},   extent = {{-80, -10}, {-60, 10}}, rotation = 0)));
  OpenIPSL.Electrical.Branches.PSAT.TwoWindingTransformer transformer(
    Sn=2220000000,
    V_b=400000,
    Vn=400000,
    rT=0,
    xT=0.15)                                                                                                                         annotation (
    Placement(visible = true, transformation(origin={92,12},   extent = {{-50, -10}, {-30, 10}}, rotation = 0)));
  inner OpenIPSL.Electrical.SystemBase SysData(S_b=100000000, fn=60)     annotation (
    Placement(visible = true, transformation(origin={-24,142},  extent = {{60, -100}, {120, -60}}, rotation = 0)));
  Modelica.Blocks.Continuous.FirstOrder firstOrder(
    k=1,
    T=0.12,
    y_start=1)
    annotation (Placement(transformation(extent={{-58,6},{-38,26}})));
  GOVCODE gOVCODE
    annotation (Placement(transformation(extent={{-58,-38},{-38,-18}})));
  AVRCode aVRCode
    annotation (Placement(transformation(extent={{-60,50},{-40,70}})));
  OpenIPSL.Electrical.Machines.PSAT.Order4 order4_1(
    Sn=2220000000,
    Vn=400000,
    ra=0.003,
    x1d=0.3,
    M=7,
    D=0,
    w(start=1),
    xd=1.81,
    xq=1.76,
    x1q=0.65,
    T1d0=8,
    T1q0=1) annotation (Placement(transformation(extent={{-8,2},{12,22}})));
equation
  connect(line_1.n,B3. p) annotation (
    Line(points={{133,32},{140,32},{140,12},{156,12}},        color = {0, 0, 255}));
  connect(B3.p,infinite_bus. p) annotation (
    Line(points={{156,12},{168,12}},                        color = {0, 0, 255}));
  connect(line_2.n,B3. p) annotation (
    Line(points={{133,-8},{140,-8},{140,12},{156,12}},          color = {0, 0, 255}));
  connect(line_1.p,B2. p) annotation (
    Line(points={{115,32},{104,32},{104,12},{78,12}},     color = {0, 0, 255}));
  connect(line_2.p,B2. p) annotation (
    Line(points={{115,-8},{104,-8},{104,12},{78,12}},       color = {0, 0, 255}));
  connect(fault.p,B2. p) annotation (
    Line(points={{112.333,-36},{100,-36},{100,12},{78,12}},      color = {0, 0, 255}));
  connect(B1.p,transformer. p) annotation (
    Line(points={{26,12},{41,12}},           color = {0, 0, 255}));
  connect(transformer.n,B2. p) annotation (
    Line(points={{63,12},{78,12}},    color = {0, 0, 255}));
  connect(aVRCode.y, firstOrder.u)
    annotation (Line(points={{-35,60},{-32,60},{-32,74},{-72,74},{-72,16},{-60,
          16}},                                    color={0,0,127}));
  connect(firstOrder.y, aVRCode.Vf) annotation (Line(points={{-37,16},{-32,16},
          {-32,42},{-68,42},{-68,56},{-66,56}},    color={0,0,127}));
  connect(order4_1.p, B1.p)
    annotation (Line(points={{12,12},{26,12}}, color={0,0,255}));
  connect(order4_1.pm0, gOVCODE.u0) annotation (Line(points={{-6,1},{-6,-44},{
          -68,-44},{-68,-24},{-64,-24}},              color={0,0,127}));
  connect(gOVCODE.Tm, order4_1.pm) annotation (Line(points={{-33,-27.8},{-20,
          -27.8},{-20,7},{-10,7}},
                            color={0,0,127}));
  connect(firstOrder.y, order4_1.vf) annotation (Line(points={{-37,16},{-18,16},
          {-18,17},{-10,17}}, color={0,0,127}));
  connect(order4_1.v, aVRCode.Vter) annotation (Line(points={{13,15},{13,76},{
          -74,76},{-74,64},{-66,64}}, color={0,0,127}));
  connect(order4_1.w, gOVCODE.w) annotation (Line(points={{13,21},{14,21},{14,
          -48},{-74,-48},{-74,-32},{-64,-32}}, color={0,0,127}));
  annotation (uses(OpenIPSL(version="3.1.0-dev"), Modelica(version="4.0.0")),
      experiment(
      StopTime=30,
      __Dymola_fixedstepsize=0.001,
      __Dymola_Algorithm="Dassl"),
    Diagram(coordinateSystem(extent={{-120,-80},{200,100}})),
    Icon(coordinateSystem(extent={{-120,-80},{200,100}})));
end SMIBAVRCodeGOVCode4;
