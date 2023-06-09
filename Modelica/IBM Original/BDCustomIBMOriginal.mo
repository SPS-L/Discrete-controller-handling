within ;
model BDCustomIBM
  Modelica.Blocks.Sources.Step step
    annotation (Placement(transformation(extent={{-78,-14},{-58,8}})));
  Modelica.Blocks.Continuous.StateSpace stateSpace(
    A=[-0.04,sqrt(0.02^2 + 0.9^2); -sqrt(0.02^2 + 0.9^2),0],
    B=[-sqrt(0.02^2 + 0.9^2); 0],
    C=[1,0; 0,1],
    D=[0; 0],
    x_start={0,0})
           annotation (Placement(transformation(extent={{44,-8},{64,12}})));
  IBM.IBMController iBMController
    annotation (Placement(transformation(extent={{0,-8},{20,12}})));
  Modelica.Blocks.Continuous.Derivative derivative(k=1, T=0.0001)
    annotation (Placement(transformation(extent={{6,30},{26,50}})));
equation
  connect(iBMController.y, stateSpace.u[1]) annotation (Line(points={{17.9,1.9},
          {28.95,1.9},{28.95,2},{42,2}}, color={0,0,127}));
  connect(iBMController.uo, stateSpace.y[2]) annotation (Line(points={{2,1.8},{
          -6,1.8},{-6,16},{65,16},{65,2}}, color={0,0,127}));
  connect(stateSpace.y[2], derivative.u) annotation (Line(points={{65,2},{65,16},
          {-6,16},{-6,40},{4,40}}, color={0,0,127}));
  connect(derivative.y, iBMController.ud) annotation (Line(points={{27,40},{40,
          40},{40,66},{-20,66},{-20,6},{2,6}}, color={0,0,127}));
  connect(step.y, iBMController.uf) annotation (Line(points={{-57,-3},{-27.5,-3},
          {-27.5,-2.4},{2,-2.4}}, color={0,0,127}));
  annotation (
    Icon(coordinateSystem(preserveAspectRatio=false)),
    Diagram(coordinateSystem(preserveAspectRatio=false)),
    uses(Modelica(version="4.0.0")),
    experiment(
      StopTime=50,
      __Dymola_fixedstepsize=0.05,
      __Dymola_Algorithm="Rkfix2"));
end BDCustomIBM;
