within DynamicVCC.Components.Pipes;
model MoistAirCrossFlow "Moist air heat and mass transfer for cross flow HX"

  import Modelica.Constants.T_zero;
  import Modelica.Fluid.Utilities.regStep;

  parameter Integer Ncell=1 "Number of control volumes";

  // extending dry coil model
  extends DynamicVCC.Components.Pipes.BaseClasses.DryAirCoil(final n=Ncell);

  output SI.MassFlowRate m_flows_dehumid[n];

  parameter SI.LewisNumber Le=0.854 "Lewis number";
  parameter SI.SpecificEnthalpy delta_h_ig=2836.6e+3 "Latent heat of sublimation";

  SI.SpecificEnthalpy delta_h_fg[n]=Medium.enthalpyOfVaporization(T_a) "Latent heat of condensation";
  Medium.MassFraction w_s[n];

protected
    Real Ntu_mass[n];
    SI.HeatFlowRate Q_flows_lat[n] "Latent heat transfer rate";

equation

/******************** Mass transfer of moist air ************************/
    for i in 1:n loop
      w_s[i]=Medium.xsaturation_pT(ports_b[n].p, T_s[n]);
    end for;
    m_flows_dehumid={m_flows[i]*(w_a[i]-w_b[i]) for i in 1:n};
    Ntu_mass=Ntu/Le^(2/3);
    w_b={w_a[i]+min(0,(w_s[i]-w_a[i]))*(1-exp(-Ntu_mass[i])) for i in 1:n};

/************************ Heat transfer **************************/
    Q_flows_lat={m_flows[i]*(w_b[i]-w_a[i])*regStep(T_s[i]+T_zero,delta_h_fg[i],delta_h_ig,1e-2) for i in 1:n};
    Q_flows_forced=Q_flows_sen;//+Q_flows_lat;

  annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(coordinateSystem(preserveAspectRatio=false)),
    experiment(
      StartTime=1700,
      StopTime=5000,
      __Dymola_Algorithm="Radau"));
end MoistAirCrossFlow;
