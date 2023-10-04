within DynamicVCC.Media;
package Utilities
  extends Modelica.Icons.UtilitiesPackage;
  function phaseTransition "Smooth transitions between three phases using Fuzzy Modeling"
    // Utilize thermodynamic equilibrium quality as a linguistic variable of three Fuzzy Numbers
    extends Modelica.Icons.Function;
    import Modelica.Fluid.Utilities.regStep;
    input Real x "Thermodynamic equilibrium quality";
    output Real w1 "Liquid phase weight";
    output Real w2 "Two-phase weight";
    output Real w3 "Vapor phase weight";
  protected
    Real mu_L "membership function of Fuzzy Number L";
    Real mu_T "membership function of Fuzzy Number T";
    Real mu_V "membership function of Fuzzy Number V";
  algorithm
    mu_L:=regStep(
      x,
      0,
      1,
      1e-7);
    mu_V:=regStep(
      x - 1,
      1,
      0,
      1e-7);
    mu_T:=1 - mu_L - mu_V;
    w1:=mu_L;
    w2:=mu_T;
    w3:=mu_V;

  end phaseTransition;
end Utilities;
