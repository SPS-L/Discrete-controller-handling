within ;
model SRM1

  Modelica.Blocks.Sources.Step step
    annotation (Placement(transformation(extent={{-92,0},{-72,22}})));
  Modelica.Blocks.Math.Feedback feedback
    annotation (Placement(transformation(extent={{-56,2},{-36,22}})));
  Modelica.Blocks.Continuous.StateSpace stateSpace(
    A=[-0.04,sqrt(0.02^2 + 0.9^2); -sqrt(0.02^2 + 0.9^2),0],
    B=[-sqrt(0.02^2 + 0.9^2); 0],
    C=[1,0; 0,1],
    D=[0; 0],
    x_start={0,0})
           annotation (Placement(transformation(extent={{54,2},{74,22}})));
  Modelica.Blocks.Discrete.Sampler sampler(samplePeriod=0.01)
    annotation (Placement(transformation(extent={{-26,2},{-6,22}})));
  Modelica.Blocks.Continuous.Integrator integrator(k=0.07)
    annotation (Placement(transformation(extent={{14,2},{34,22}})));
equation
  connect(step.y,feedback. u1) annotation (Line(points={{-71,11},{-70,12},{-54,
          12}},          color={0,0,127}));
  connect(stateSpace.y[2],feedback. u2) annotation (Line(points={{75,12},{80,12},
          {80,-6},{-46,-6},{-46,4}},    color={0,0,127}));
  connect(feedback.y, sampler.u)
    annotation (Line(points={{-37,12},{-28,12}}, color={0,0,127}));
  connect(sampler.y, integrator.u)
    annotation (Line(points={{-5,12},{12,12}}, color={0,0,127}));
  connect(integrator.y, stateSpace.u[1])
    annotation (Line(points={{35,12},{52,12}}, color={0,0,127}));
  annotation (uses(Modelica(version="4.0.0")), experiment(
      StopTime=100,
      __Dymola_fixedstepsize=0.01,
      __Dymola_Algorithm="Rkfix2"));
end SRM1;
