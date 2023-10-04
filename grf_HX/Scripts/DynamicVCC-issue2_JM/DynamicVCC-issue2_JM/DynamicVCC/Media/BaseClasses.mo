within DynamicVCC.Media;
package BaseClasses

  extends Modelica.Media.Interfaces.PartialTwoPhaseMedium;

  extends Modelica.Icons.BasesPackage;
  redeclare replaceable record extends ThermodynamicState

    AbsolutePressure p "Pressure";
    SpecificEnthalpy h "Specific enthalpy";
    Density d "Density";
    Temperature T "Temperature";

  end ThermodynamicState;

  redeclare replaceable record extends SaturationProperties
    annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(coordinateSystem(preserveAspectRatio=false)));
  end SaturationProperties;

  redeclare function extends density
  algorithm
    d:=state.d;
    annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(coordinateSystem(preserveAspectRatio=false)));
  end density;

  redeclare replaceable function extends temperature
  algorithm
    T:=state.T;
    annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(coordinateSystem(preserveAspectRatio=false)));
  end temperature;

  redeclare function extends specificEnthalpy

  algorithm
    h:=state.h;
    annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(coordinateSystem(preserveAspectRatio=false)));
  end specificEnthalpy;

  redeclare function extends specificInternalEnergy
  algorithm
    u:=specificEnthalpy(state)-pressure(state)/density(state);
    annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(coordinateSystem(preserveAspectRatio=false)));
  end specificInternalEnergy;

  redeclare function extends pressure
  algorithm
    p:=state.p;
    annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(coordinateSystem(preserveAspectRatio=false)));
  end pressure;

  redeclare partial function specificEnthalpy_ps
    input AbsolutePressure p;
    input SpecificEntropy s;
    input FixedPhase phase=1;
    output SpecificEnthalpy h;
    annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(coordinateSystem(preserveAspectRatio=false)));
  end specificEnthalpy_ps;

  replaceable partial function pressure_du
    input Density d;
    input SpecificInternalEnergy u;
    output AbsolutePressure p;

    annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(coordinateSystem(preserveAspectRatio=false)));
  end pressure_du;

  redeclare replaceable partial function extends density_derp_h
    extends Modelica.Icons.Function;
  protected
    AbsolutePressure p=pressure(state);
    SpecificEnthalpy h=specificEnthalpy(state);
    Density rho=density(state);
    SaturationProperties sat(psat=p,Tsat=T_default);
    SpecificEnthalpy hl=bubbleEnthalpy(sat);
    SpecificEnthalpy hv=dewEnthalpy(sat);
    Density rhof=bubbleDensity(sat);
    Density rhog=dewDensity(sat);
    Real dxdp;
    DerDensityByPressure drhodp_l "subcooled liquid";
    DerDensityByPressure drhodp_v "superheated vapor";
    DerDensityByPressure drhodp_tp "two-phase";
    DerDensityByPressure drhodp_f=dBubbleDensity_dPressure(sat) "bubble line";
    DerDensityByPressure drhodp_g=dDewDensity_dPressure(sat) "dew line";
    DerEnthalpyByPressure dhdp_f=dBubbleEnthalpy_dPressure(sat) "bubble line";
    DerEnthalpyByPressure dhdp_g=dDewEnthalpy_dPressure(sat) "dew line";
    Real w1;
    Real w2;
    Real w3;

    annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(coordinateSystem(preserveAspectRatio=false)));
  end density_derp_h;

  redeclare replaceable partial function extends density_derh_p
    extends Modelica.Icons.Function;
  protected
    AbsolutePressure p=pressure(state);
    SpecificEnthalpy h=specificEnthalpy(state);
    Density rho=density(state);
    SaturationProperties sat(psat=p,Tsat=T_default);
    SpecificEnthalpy hl=bubbleEnthalpy(sat);
    SpecificEnthalpy hv=dewEnthalpy(sat);
    Density rhof=bubbleDensity(sat);
    Density rhog=dewDensity(sat);
    DerDensityByEnthalpy drhodh_l "Subcooled liquid";
    DerDensityByEnthalpy drhodh_v "Superheated vapor";
    DerDensityByEnthalpy drhodh_tp;
    Real w1;
    Real w2;
    Real w3;
    annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(coordinateSystem(preserveAspectRatio=false)));
  end density_derh_p;

  redeclare function extends setState_phX

  protected
    SaturationProperties sat(psat=p,Tsat=T_default);
    SpecificEnthalpy hl=bubbleEnthalpy(sat);
    SpecificEnthalpy hv=dewEnthalpy(sat);
  algorithm
    state.p:=p;
    state.h:=h;
    state.d := density_ph(p, h);
    state.T := temperature_ph(p, h);
    state.phase:=if (h < hl) or (h > hv) then 1 else 2;
    annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(coordinateSystem(preserveAspectRatio=false)));
  end setState_phX;

  redeclare replaceable partial function extends thermalConductivity
  protected
    ThermalConductivity lambda_l "subcooled liquid";
    ThermalConductivity lambda_v "superheated vapor";
    ThermalConductivity lambda_tp "two-phase";
    ThermalConductivity lambda_f "saturated liquid";
    ThermalConductivity lambda_g "saturated vapor";
    AbsolutePressure p=pressure(state);
    SpecificEnthalpy h=specificEnthalpy(state);
    SaturationProperties sat(psat=p,Tsat=T_default);
    SpecificEnthalpy hl=bubbleEnthalpy(sat);
    SpecificEnthalpy hv=dewEnthalpy(sat);
    Real x=(h-hl)/(hv-hl) "Thermodynamic equilibrium quality";
    Real w1;
    Real w2;
    Real w3;

    annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(coordinateSystem(preserveAspectRatio=false)));
  end thermalConductivity;

  redeclare replaceable partial function extends dynamicViscosity
  protected
    DynamicViscosity eta_l "subcooled liquid";
    DynamicViscosity eta_tp "two-phase";
    DynamicViscosity eta_v "superheated vapor";
    DynamicViscosity eta_f "saturated liquid";
    DynamicViscosity eta_g "saturated vapor";
    AbsolutePressure p=pressure(state);
    SpecificEnthalpy h=specificEnthalpy(state);
    SaturationProperties sat(psat=p,Tsat=T_default);
    SpecificEnthalpy hl=bubbleEnthalpy(sat);
    SpecificEnthalpy hv=dewEnthalpy(sat);
    Real x=(h-hl)/(hv-hl) "Thermodynamic equilibrium quality";
    Real w1;
    Real w2;
    Real w3;
    annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(coordinateSystem(preserveAspectRatio=false)));
  end dynamicViscosity;

  redeclare replaceable partial function extends specificHeatCapacityCp
  protected
    SpecificHeatCapacity cp_f "Saturated liquid";
    SpecificHeatCapacity cp_g "Saturated vapor";
    SpecificHeatCapacity cp_l "subcooled liquid";
    SpecificHeatCapacity cp_tp "two-phase";
    SpecificHeatCapacity cp_v "superheated vapor";
    AbsolutePressure p=pressure(state);
    SpecificEnthalpy h=specificEnthalpy(state);
    SaturationProperties sat(psat=p,Tsat=T_default);
    SpecificEnthalpy hl=bubbleEnthalpy(sat);
    SpecificEnthalpy hv=dewEnthalpy(sat);
    Real x=(h-hl)/(hv-hl) "Thermodynamic equilibrium quality";
    Real w1;
    Real w2;
    Real w3;

    annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(coordinateSystem(preserveAspectRatio=false)));
  end specificHeatCapacityCp;

  replaceable partial function vaporSpecificEnthalpy_pd
    extends Modelica.Icons.Function;
    input AbsolutePressure p;
    input Density d;
    output SpecificEnthalpy h;
    annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(coordinateSystem(preserveAspectRatio=false)));
  end vaporSpecificEnthalpy_pd;

  replaceable partial function liquidSpecificEnthalpy_pd
    extends Modelica.Icons.Function;
    input AbsolutePressure p;
    input Density d;
    output SpecificEnthalpy h;
    annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(coordinateSystem(preserveAspectRatio=false)));
  end liquidSpecificEnthalpy_pd;
end BaseClasses;
