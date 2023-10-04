within DynamicVCC.Utilities;
function regPowGen
  "Generalized power approximation with given derivative in the origin (Qiao, 2022)"
  extends Modelica.Icons.Function;
  input Real x;
  input Real a;
  input Real delta=0.01 "Range of significant deviation from x^a*sgn(x)";
  input Real b=3;
  output Real y;
algorithm
  y := x^b*(x*x+delta*delta)^((a-b)/2);

end regPowGen;
