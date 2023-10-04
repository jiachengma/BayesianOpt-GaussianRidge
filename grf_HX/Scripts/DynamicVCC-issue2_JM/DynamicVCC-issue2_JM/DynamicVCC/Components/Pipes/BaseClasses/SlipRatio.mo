within DynamicVCC.Components.Pipes.BaseClasses;
package SlipRatio "Slip ratio model for two-phase slip flow"
  extends Modelica.Icons.VariantsPackage;

  partial model PartialSlipRatio "Partial model for slip ratio"
    replaceable package Medium=Modelica.Media.Interfaces.PartialTwoPhaseMedium;

    parameter Integer n=1 "Number of control volumes";

    input Medium.ThermodynamicState states[n];

    output Real S[n] "Slip ratio";

  protected
    Medium.Density rho_l[n]=Medium.bubbleDensity(Medium.setSat_p(Medium.pressure(states)));
    Medium.Density rho_v[n]=Medium.dewDensity(Medium.setSat_p(Medium.pressure(states)));

    annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(coordinateSystem(preserveAspectRatio=false)));
  end PartialSlipRatio;

  model Zivi "Zivi (1964)"
    import Modelica.Fluid.Utilities.regPow;
    extends DynamicVCC.Components.Pipes.BaseClasses.SlipRatio.PartialSlipRatio;
  equation
    S={regPow(rho_l[i]/rho_v[i],1/3) for i in 1:n};
    annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(coordinateSystem(preserveAspectRatio=false)));
  end Zivi;

  model Homogeneous "Homogeneous two-phase flow"
    extends DynamicVCC.Components.Pipes.BaseClasses.SlipRatio.PartialSlipRatio;
  equation
    S=ones(n);
    annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(coordinateSystem(preserveAspectRatio=false)));
  end Homogeneous;

  model Smith "Smith correlation of slip ratio (1969)"

    extends DynamicVCC.Components.Pipes.BaseClasses.SlipRatio.PartialSlipRatio;

    parameter Real K=0.4;
    parameter Real x_small=1e-3;
  protected
    Real S[n] "Slip ratio";
    Real x_2phase[n] "Set to 1 for single-phase";

  equation
    x_2phase= {smooth(0,noEvent(if x[i]<x_small then x_small else x[i])) for i in 1:n};
    S={K+(1-K)*sqrt((rho_f[i]/rho_g[i]+K*(1-x_2phase[i])/x_2phase[i])/(1+K*(1-x_2phase[i])/x_2phase[i])) for i in 1:n};
    gamma={x[i]/(x[i]+(1-x[i])*rho_g[i]/rho_f[i]*S[i]) for i in 1:n};

    annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(coordinateSystem(preserveAspectRatio=false)));
  end Smith;
end SlipRatio;
