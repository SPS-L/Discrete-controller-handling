within ;
model SMIBAVRIBMGOVIBM4

  OpenIPSL.Electrical.Buses.InfiniteBus infinite_bus(
    P_0=-1998000000,
    Q_0=87066000,
    angle_0=0,
    displayPF=false,
    v_0=0.90081)                                                                                                                       annotation (
    Placement(visible = true, transformation(origin={176,22},   extent = {{-10, 10}, {10, -10}}, rotation = 180)));
  OpenIPSL.Electrical.Branches.PwLine line_1(
    B=0,
    G=0,
    R=0,
    S_b=100000000,
    X=0.5*100/2220,
    displayPF=false)                                                                                                    annotation (
    Placement(visible = true, transformation(origin={90,22},   extent = {{20, 10}, {40, 30}}, rotation = 0)));
  OpenIPSL.Electrical.Branches.PwLine line_2(
    B=0,
    G=0,
    R=0,
    S_b=100000000,
    X=0.93*100/2220)                                                                                  annotation (
    Placement(visible = true, transformation(origin={92,22},   extent = {{20, -30}, {40, -10}}, rotation = 0)));
  OpenIPSL.Electrical.Events.PwFault fault(
    R=0,
    X=0.01*100/2220,
    t1=2,
    t2=2.1)                                                                               annotation (
    Placement(visible = true, transformation(origin={94,24},   extent = {{20, -60}, {40, -40}}, rotation = 0)));
  OpenIPSL.Electrical.Buses.Bus B3(displayPF=false, v_0=0.900810000000000)     annotation (
    Placement(visible = true, transformation(origin={88,22},   extent = {{60, -10}, {80, 10}}, rotation = 0)));
  OpenIPSL.Electrical.Buses.Bus B2(
    angle_0=0.351222533717268,
    displayPF=false,
    v_0=0.944299492912195)                                                                                  annotation (
    Placement(visible = true, transformation(origin={88,22},   extent = {{-20, -10}, {0, 10}}, rotation = 0)));
  OpenIPSL.Electrical.Buses.Bus B1(angle_0=0.494677176989154, displayPF=false)     annotation (
    Placement(visible = true, transformation(origin={108,22},  extent = {{-80, -10}, {-60, 10}}, rotation = 0)));
  OpenIPSL.Electrical.Branches.PSAT.TwoWindingTransformer transformer(
    Sn=2220000000,
    V_b=400000,
    Vn=400000,
    rT=0,
    xT=0.15)                                                                                                                         annotation (
    Placement(visible = true, transformation(origin={100,22},  extent = {{-50, -10}, {-30, 10}}, rotation = 0)));
  inner OpenIPSL.Electrical.SystemBase SysData(S_b=100000000, fn=60)     annotation (
    Placement(visible = true, transformation(origin={-26,64},   extent = {{60, -100}, {120, -60}}, rotation = 0)));
  Modelica.Blocks.Continuous.FirstOrder firstOrder(
    k=1,
    T=0.12,
    y_start=1)
    annotation (Placement(transformation(extent={{-32,18},{-12,38}})));
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
    T1q0=1) annotation (Placement(transformation(extent={{0,12},{20,32}})));
  GOVIBM gOVIBM
    annotation (Placement(transformation(extent={{-34,-14},{-14,6}})));
  AVRIBM aVRIBM
    annotation (Placement(transformation(extent={{-70,18},{-50,38}})));
equation
  connect(line_1.n,B3. p) annotation (
    Line(points={{129,42},{142,42},{142,22},{158,22}},        color = {0, 0, 255}));
  connect(B3.p,infinite_bus. p) annotation (
    Line(points={{158,22},{166,22}},                        color = {0, 0, 255}));
  connect(line_2.n,B3. p) annotation (
    Line(points={{131,2},{142,2},{142,22},{158,22}},            color = {0, 0, 255}));
  connect(line_1.p,B2. p) annotation (
    Line(points={{111,42},{104,42},{104,22},{78,22}},     color = {0, 0, 255}));
  connect(line_2.p,B2. p) annotation (
    Line(points={{113,2},{104,2},{104,22},{78,22}},         color = {0, 0, 255}));
  connect(fault.p,B2. p) annotation (
    Line(points={{112.333,-26},{98,-26},{98,22},{78,22}},        color = {0, 0, 255}));
  connect(B1.p,transformer. p) annotation (
    Line(points={{38,22},{49,22}},           color = {0, 0, 255}));
  connect(transformer.n,B2. p) annotation (
    Line(points={{71,22},{78,22}},    color = {0, 0, 255}));
  connect(order4_1.p, B1.p)
    annotation (Line(points={{20,22},{38,22}},  color={0,0,255}));
  connect(firstOrder.y, order4_1.vf) annotation (Line(points={{-11,28},{-11,27},
          {-2,27}},                color={0,0,127}));
  connect(gOVIBM.y, order4_1.pm) annotation (Line(points={{-13,-4},{-4,-4},{-4,
          18},{-2,18},{-2,17}},
                           color={0,0,127}));
  connect(aVRIBM.y, firstOrder.u)
    annotation (Line(points={{-49,28},{-34,28}},   color={0,0,127}));
  connect(order4_1.pm0, gOVIBM.u0) annotation (Line(points={{2,11},{2,-24},{-48,
          -24},{-48,0},{-36,0}},           color={0,0,127}));
  connect(firstOrder.y, aVRIBM.uf0) annotation (Line(points={{-11,28},{-8,28},{
          -8,10},{-78,10},{-78,24},{-72,24}},     color={0,0,127}));
  connect(order4_1.v, aVRIBM.uter) annotation (Line(points={{21,25},{32,25},{32,
          44},{-78,44},{-78,32},{-72,32}}, color={0,0,127}));
  connect(order4_1.w, gOVIBM.uf) annotation (Line(points={{21,31},{26,31},{26,
          -20},{-42,-20},{-42,-8},{-36,-8}}, color={0,0,127}));
  annotation (uses(OpenIPSL(version="3.1.0-dev"), Modelica(version="4.0.0")),
      experiment(
      StopTime=30,
      __Dymola_fixedstepsize=0.005,
      __Dymola_Algorithm="Dassl"),
    Diagram(coordinateSystem(extent={{-180,-100},{200,100}})),
    Icon(coordinateSystem(extent={{-180,-100},{200,100}})));
end SMIBAVRIBMGOVIBM4;
