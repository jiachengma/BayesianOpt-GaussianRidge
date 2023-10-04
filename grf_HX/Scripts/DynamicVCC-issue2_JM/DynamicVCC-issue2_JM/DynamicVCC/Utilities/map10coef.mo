within DynamicVCC.Utilities;
function map10coef
  input Real map[:,:];
  input Real Te "K";
  input Real Tc "K";
  output Real y[size(map,1)+1];
protected
  Real S;
  Real D;
  Real C[:];
algorithm
  y[1]:=0;
  for i in 1:size(map,1) loop
    C:=map[i,:];
    S := Utilities.K2F(Te);
    D := Utilities.K2F(Tc);
    y[i+1]:=C[1] + (C[2]*S) + (C[3]*D) + (C[4]*S^2) + (C[5]*S*D) + (C[6]*D^2) + (C[7]
      *S^3) + (C[8]*D*S^2) + (C[9]*S*D^2) + (C[10]*D^3);
  end for;
end map10coef;
