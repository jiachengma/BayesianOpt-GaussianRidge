within DynamicVCC.Components.Units.HX;
model Piping "Connecting pipe"

    import Modelica.Constants.pi;

    input Medium_2.Temperature T_amb=298 "Ambient temperature";
    input Medium_2.AbsolutePressure p_amb=system.p_ambient "Ambient pressure";

    replaceable package Medium_2=Modelica.Media.Air.SimpleAir;

    extends DynamicVCC.Components.Units.HX.BaseClasses.PartialHX(
    final As_1=surfaceArea_i,
    final Ac_1=crossArea,
    final L_1=length,
    final diameter_1=d_i,
    final C_metalWall=C_tube);

    // Pipe geometry
    parameter SI.Diameter d_i "Tube outter diameter";
    parameter SI.Diameter d_o "Tube inner diameter";
    parameter SI.Length length "Tube length";
    parameter SI.Density rho_tube=8900 "Tube density";
    parameter SI.SpecificHeatCapacity cp_tube=385 "Tube material specific heat";

protected
    parameter SI.Area surfaceArea_i=pi*d_i*length "Tube inner surface area";
    parameter SI.Area surfaceArea_o=pi*d_o*length "Tube outter surface area";
    parameter SI.Area crossArea=d_i^2*pi/4 "cross-sectional area";
    parameter SI.Volume tubeVolume=(d_o^2-d_i^2)*pi/4*length "Tube wall volume";
    parameter SI.Mass m_tube=tubeVolume*rho_tube "Tube wall mass";
    parameter SI.HeatCapacity C_tube=m_tube*cp_tube "Tube wall heat capacity";
    Medium_2.ThermodynamicState states_amb[Ncell] "Ambient air states";

public
    replaceable model HeatTransfer_2 = DynamicVCC.Components.Pipes.BaseClasses.HeatTransfer.ConstantFlowHeatTransfer;

    // Heat loss to ambient
    HeatTransfer_2 heatLoss(
    redeclare final package Medium=Medium_2,
    final n=Ncell,
    final states=states_amb,
    final surfaceAreas=fill(surfaceArea_o/Ncell,Ncell),
    final dimensions=fill(d_o,Ncell),
    final lengths=fill(length/Ncell,Ncell),
    final vs=zeros(Ncell));

equation
    states_amb=Medium_2.setState_pT(fill(p_amb,Ncell),fill(T_amb,Ncell));

    connect(metalWall.heatPorts_b,heatLoss.heatPorts);
    Q_flow_2=sum(heatLoss.Q_flows);

    annotation (Icon(coordinateSystem(preserveAspectRatio=false), graphics={Rectangle(
            extent={{-100,-40},{100,40}},
            lineColor={0,0,0},
            lineThickness=1,
            fillColor={90,77,145},
            fillPattern=FillPattern.HorizontalCylinder),                    Text(
            extent={{-56,20},{66,142}},
            textColor={0,0,0},
            textString="%name")}),                         Diagram(coordinateSystem(preserveAspectRatio=false)),
              Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(coordinateSystem(preserveAspectRatio=false)));
end Piping;
