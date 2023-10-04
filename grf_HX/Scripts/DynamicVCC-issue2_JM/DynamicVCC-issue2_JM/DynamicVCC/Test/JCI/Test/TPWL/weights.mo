within DynamicVCC.Test.JCI.Test.TPWL;
function weights
  extends Modelica.Icons.Function;
  import Modelica.Math.Vectors.norm;
  input Real x[:];
  input Real xe1[:];
  input Real xe2[:];
  input Real beta=25;
  output Real w[2];
protected
  Real d[2];
  Real m;
  Real S;

algorithm
  d[1]:=norm(x - xe1, p=2);
  d[2]:=norm(x - xe2, p=2);
  m:=min(d)+1e-10;
  w[1]:=exp(-beta*d[1]/m);
  w[2]:=exp(-beta*d[2]/m);
  S:=sum(w);
  w[1]:=w[1]/(S);
  w[2]:=w[2]/(S);

end weights;
