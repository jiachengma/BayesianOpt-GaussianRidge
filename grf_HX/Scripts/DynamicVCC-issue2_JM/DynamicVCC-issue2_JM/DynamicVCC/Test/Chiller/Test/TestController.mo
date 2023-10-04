within DynamicVCC.Test.Chiller.Test;
model TestController
  extends Modelica.Icons.Example;

  DynamicVCC.Test.Chiller.Controllernew controller(gamma_init=0.9096);

  Modelica.Blocks.Sources.CombiTimeTable Mea(tableOnFile=true,smoothness=Modelica.Blocks.Types.Smoothness.ContinuousDerivative,tableName="Mea",fileName="C:/Jiacheng Ma/BoundaryCondition/Chiller/Mea.mat",columns=2:10);
  Modelica.Blocks.Sources.CombiTimeTable BC(tableOnFile=true,smoothness=Modelica.Blocks.Types.Smoothness.ContinuousDerivative,tableName="BC",fileName="C:/Jiacheng Ma/BoundaryCondition/Chiller/BC.mat",columns=2:10);

equation
  controller.u1=Mea.y[3];
  controller.u2=BC.y[9];

  annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(coordinateSystem(preserveAspectRatio=false)),
    experiment(
      StartTime=2000,
      StopTime=50000,
      __Dymola_Algorithm="Dassl"));
end TestController;
