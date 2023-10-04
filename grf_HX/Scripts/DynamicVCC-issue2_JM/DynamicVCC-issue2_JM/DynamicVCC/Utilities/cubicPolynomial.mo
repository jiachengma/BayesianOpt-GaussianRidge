within DynamicVCC.Utilities;
function cubicPolynomial "Cubic polynomial interpolation given function values and derivatives at two end points"
  extends Modelica.Icons.Function;

  input Real x;
  input Real x1;
  input Real x2;
  input Real y1;
  input Real y2;
  input Real dy1;
  input Real dy2;
  output Real y;
protected
  Real diff_x=x2 - x1;
  Real m=(y2 - y1)/diff_x;
  Real c1=dy1;
  Real c2=(3*m - 2*dy1 - dy2)/diff_x;
  Real c3=(dy1 + dy2 - 2*m)/(diff_x*diff_x);
  Real dx;
algorithm
  dx:=x - x1;
  y:=y1 + c1*dx + c2*dx^2 + c3*dx^3;

  annotation (Documentation(info="<html>
<p><img src=\"modelica://DynamicVCC/Resources/Images/equations/equation-qz4yGt4b.png\" alt=\"y=y_1+c_1*(x-x_1)+c_2*(x-x_1)^2+c_3*(x-x_1)^3\"/></p>
</html>"));
end cubicPolynomial;
