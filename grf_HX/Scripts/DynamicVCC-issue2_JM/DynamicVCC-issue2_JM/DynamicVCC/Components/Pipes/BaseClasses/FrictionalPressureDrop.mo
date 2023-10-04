within DynamicVCC.Components.Pipes.BaseClasses;
package FrictionalPressureDrop "fin-and-tube coil air-side frictional pressure drop"
  extends Modelica.Icons.VariantsPackage;

  partial model PartialFrictionalPressureDrop
    replaceable package Medium=Modelica.Media.Interfaces.PartialMedium;

    parameter Integer n=1 "Number of control volumes";

    input Medium.ThermodynamicState states[n] "Fluid states";

    input SI.Area surfaceAreas[n] "Heat transfer area";

    input SI.Area crossAreas[n] "Cross-sectional area";

    input SI.Length dimensions[n] "Characteristic dimension, typically plate length or tube diameter";

    input SI.Length lengths[n] "Length along flow path";

    input SI.Velocity vs[n] "Mean velocity of each control volume";

    output SI.Pressure dps[n] "Frictional pressure drop";
    annotation (Icon(coordinateSystem(preserveAspectRatio=false), graphics={
                            Line(
            points={{-80,-60},{-80,60},{80,-60},{80,62}},
            color={0,0,255},
            thickness=1)}),                                        Diagram(coordinateSystem(preserveAspectRatio=false)));
  end PartialFrictionalPressureDrop;

  model PlainFinAndTube_Wang "Air flow friction of fin-and-tube heat exchangers (Wang, 2000)"

    import Modelica.Fluid.Utilities.regPow;
    import Modelica.Fluid.Pipes.BaseClasses.CharacteristicNumbers.ReynoldsNumber;
    extends DynamicVCC.Components.Pipes.BaseClasses.FrictionalPressureDrop.PartialFrictionalPressureDrop;

    // Coil geometry
    parameter Integer nRow=1 "Number of tube rows";

    parameter SI.Length pf=1 "Fin pitch";

    parameter SI.Length pl=1 "Tube longitudinal tube pitch";

    parameter SI.Length pt=1 "Tube transverse tube pitch";

    parameter SI.Thickness t_fin=1 "Fin thickness";

    parameter SI.Diameter Dh "hydraulic diamters";

    parameter SI.CoefficientOfFriction f0=0.1;

  protected
    Medium.Density rhos[n]=Medium.density(states);
    Medium.DynamicViscosity mus[n]=Medium.dynamicViscosity(states);
    SI.Diameter Dc[n];
    Real F1[n],F2[n],F3[n] "Correlation parameters";
    SI.CoefficientOfFriction f[n](each start=f0) "Friction factor";
    SI.ReynoldsNumber Re_Dc[n];
    Real G[n] "Mass flux";
  equation
    Dc=dimensions+2*fill(t_fin,n);
    for i in 1:n loop
      Re_Dc[i]=noEvent(max(1e-6,ReynoldsNumber(vs[i],rhos[i],mus[i],Dc[i])));
      G[i]=vs[i]*rhos[i];
      F1[i]=-0.764 + 0.739*pt/pl + 0.177*pf/Dc[i] - 0.00758/nRow;
      F2[i]=-15.689 + 64.021/log(Re_Dc[i]);
      F3[i]=1.696 - 15.695/log(Re_Dc[i]);
      f[i]=0.0267*regPow(Re_Dc[i],F1[i],1)*(pt/pl)^F2[i]*(pf/Dc[i])^F3[i];
      dps[i]=f[i]*surfaceAreas[i]/crossAreas[i]*G[i]^2/(2*rhos[i]);
    end for;

    annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(coordinateSystem(preserveAspectRatio=false)));
  end PlainFinAndTube_Wang;

  model ConstantFriction "Constant friction factor"
    extends DynamicVCC.Components.Pipes.BaseClasses.FrictionalPressureDrop.PartialFrictionalPressureDrop;

    parameter SI.CoefficientOfFriction f0=0.1;
  protected
    Medium.Density rhos[n]=Medium.density(states);
    Real G[n] "Mass flux";
  equation
    G=vs.*rhos;
    dps={f0*surfaceAreas[i]/crossAreas[i]*G[i]^2/(2*rhos[i]) for i in 1:n};
    annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(coordinateSystem(preserveAspectRatio=false)));
  end ConstantFriction;
end FrictionalPressureDrop;
