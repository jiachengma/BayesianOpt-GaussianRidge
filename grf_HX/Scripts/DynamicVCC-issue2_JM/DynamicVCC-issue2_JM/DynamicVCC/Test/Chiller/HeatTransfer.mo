within DynamicVCC.Test.Chiller;
package HeatTransfer "Heat transfer correlations for chiller"
  extends Modelica.Icons.VariantsPackage;
  model SinglePhase
    extends DynamicVCC.Components.Pipes.BaseClasses.HeatTransfer_old.PartialHeatTransferCorrelation;

    parameter Real C_sf=1 "Surface enhancement factor";

  equation

    Nu=C_sf*sqrt(Re).*Pr;

    annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(coordinateSystem(preserveAspectRatio=false)));
  end SinglePhase;

  model WaterEvap "Evaporator water-side heat transfer correlation"
    extends DynamicVCC.Test.Chiller.HeatTransfer.IncompressibleCorrelation;

    parameter Real C_sf=1 "Surface enhancement factor";

  protected
    SI.CoefficientOfFriction f[n];

  equation
    f={-0.5087+0.2768*log10(Re[i])-0.0339*log10(Re[i])*log10(Re[i]) for i in 1:n};
    Nu={C_sf*(f[i]/8*(Re[i]-1000)*Pr[i])/(1.07+12.7*sqrt(f[i]/8)*(Pr[i]^0.67-1)) for i in 1:n};
    annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(coordinateSystem(preserveAspectRatio=false)));
  end WaterEvap;

  model WaterCond "Condenser water-side heat transfer correlation"
    extends DynamicVCC.Test.Chiller.HeatTransfer.IncompressibleCorrelation;

    parameter Real C_sf=1 "Surface enhancement factor";

  protected
    SI.CoefficientOfFriction f[n];

  equation
    f={1/(0.79*log(Re[i])-1.64)^2 for i in 1:n};
    Nu={C_sf*f[i]*Re[i]*Pr[i]/8/(1.07+12.7*sqrt(f[i]/8)*(Pr[i]^0.67-1)) for i in 1:n};
    annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(coordinateSystem(preserveAspectRatio=false)));
  end WaterCond;

  partial model IncompressibleCorrelation
    extends DynamicVCC.Components.Pipes.BaseClasses.HeatTransfer_old.PartialFlowHeatTransfer;

     // Fluid properties
    Medium.Density rho[n]=Medium.density(states);
    Medium.DynamicViscosity mu[n]=Medium.dynamicViscosity(states);
    Medium.PrandtlNumber Pr[n]=Medium.prandtlNumber(states);
    Medium.ThermalConductivity k[n]=Medium.thermalConductivity(states);


    // Variables
    SI.ReynoldsNumber Re[n];
    SI.NusseltNumber Nu[n];

  equation
    Re=Modelica.Fluid.Pipes.BaseClasses.CharacteristicNumbers.ReynoldsNumber(v,rho,mu,dimension);
    Nu=Modelica.Fluid.Pipes.BaseClasses.CharacteristicNumbers.NusseltNumber(alphas,dimension,k);
    annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(coordinateSystem(preserveAspectRatio=false)));
  end IncompressibleCorrelation;
end HeatTransfer;
