within ;
model BDCustomIBM
  Modelica.Blocks.Sources.Step step
    annotation (Placement(transformation(extent={{-60,-16},{-40,6}})));
  Modelica.Blocks.Continuous.StateSpace stateSpace(
    A=[-0.04,sqrt(0.02^2 + 0.9^2); -sqrt(0.02^2 + 0.9^2),0],
    B=[-sqrt(0.02^2 + 0.9^2); 0],
    C=[1,0; 0,1],
    D=[0; 0],
    x_start={0,0})
           annotation (Placement(transformation(extent={{34,-10},{54,10}})));
  IBM.IBMController iBMController
    annotation (Placement(transformation(extent={{-10,-10},{10,10}})));
equation
  connect(iBMController.y, stateSpace.u[1]) annotation (Line(points={{16.1,0.1},
          {18.95,0.1},{18.95,0},{32,0}}, color={0,0,127}));
  connect(iBMController.uo, stateSpace.y[2]) annotation (Line(points={{-16,4},{
          -26,4},{-26,-20},{60,-20},{60,0},{55,0}},
                                           color={0,0,127}));
  connect(step.y, iBMController.uf) annotation (Line(points={{-39,-5},{-36,-4},
          {-16,-4}},              color={0,0,127}));
  annotation (
    Icon(coordinateSystem(preserveAspectRatio=false)),
    Diagram(coordinateSystem(preserveAspectRatio=false)),
    uses(Modelica(version="4.0.0")),
    experiment(
      StopTime=100,
      __Dymola_fixedstepsize=0.05,
      __Dymola_Algorithm="Dassl"));
end BDCustomIBM;
