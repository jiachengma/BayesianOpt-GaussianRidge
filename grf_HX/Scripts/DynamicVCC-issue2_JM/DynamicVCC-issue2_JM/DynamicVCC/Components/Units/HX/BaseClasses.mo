within DynamicVCC.Components.Units.HX;
package BaseClasses
  extends Modelica.Icons.BasesPackage;

  partial model PartialHX "Partial HX model of 1-D working fluid flow and tube wall"

    outer DynamicVCC.Components.System system;

    import DynamicVCC.Components.Types.DifferentialState;

    // Number of control volumes
    parameter Integer Ncell=10;

    // Refrigerant medium
    replaceable package Medium_1 =
        Modelica.Media.Interfaces.PartialTwoPhaseMedium                          "Working fluid";

    // Geometry
    parameter SI.Area As_1 "Refrigerant side surface area";
    parameter SI.Area Ac_1 "Refrigerant (tube) cross area";
    parameter SI.Length L_1;
    parameter SI.Diameter diameter_1;
    parameter SI.HeatCapacity C_metalWall "Total heat capacity of tube wall";

    // Discretization
    parameter DynamicVCC.Components.Types.ModelStructure modelStructure=Types.ModelStructure.av_vb;
    parameter DifferentialState differentialState=DifferentialState.pdh;

    // Numerical
    parameter Boolean useLumpedPressure=false;
    parameter Boolean upStreamProperties=true;

    // Initialization
    parameter Medium_1.MassFlowRate m_flow_init=system.m_flow_init;
    parameter Medium_1.AbsolutePressure p_a_start=Medium_1.p_default;
    parameter Medium_1.AbsolutePressure p_b_start=p_a_start;
    parameter Medium_1.SpecificEnthalpy h_init[Ncell]=fill(Medium_1.h_default,Ncell);
    parameter SI.Temperature Tt_init[Ncell]=fill(Medium_1.T_default,Ncell);

    /********************* Connectors *****************************/
    DynamicVCC.Interfaces.FluidPort_a port_a1(
    redeclare final package Medium=Medium_1) annotation (Placement(transformation(extent={{-110,-10},{-90,10}})));
    DynamicVCC.Interfaces.FluidPort_b port_b1(
    redeclare final package Medium=Medium_1) annotation (Placement(transformation(extent={{90,-10},{110,10}}), iconTransformation(extent={{90,-10},{110,10}})));

    /********************* Components *****************************/

    // Heat transfer model
    replaceable model HeatTransfer_1 =DynamicVCC.Components.Pipes.BaseClasses.HeatTransfer.ConstantFlowHeatTransfer
      constrainedby DynamicVCC.Components.Pipes.BaseClasses.HeatTransfer.PartialFlowHeatTransfer
      "Refrigerant side heat transfer model";


    // Slip ratio
    replaceable model SlipRatio=DynamicVCC.Components.Pipes.BaseClasses.SlipRatio.Homogeneous
      constrainedby DynamicVCC.Components.Pipes.BaseClasses.SlipRatio.PartialSlipRatio;

    // Flow model
    replaceable model FlowModel_1 = DynamicVCC.Components.Pipes.BaseClasses.FlowModels.DetailedFlow (
    final upStreamProperties=upStreamProperties)
      constrainedby DynamicVCC.Components.Pipes.BaseClasses.FlowModels.PartialStaggeredGridFlowModel
    "Momentum balances at staggered grid";

    replaceable model RefFlow1D = DynamicVCC.Components.Pipes.RefFlow1D;
    // Refrigerant flow
    RefFlow1D refFlow(
      redeclare package Medium = Medium_1,
      final Ncell=Ncell,
      final lengths=fill(L_1/Ncell, Ncell),
      final crossAreas=fill(Ac_1, Ncell),
      final surfaceAreas=fill(As_1/Ncell, Ncell),
      final diameters=fill(diameter_1, Ncell),
      final modelStructure=modelStructure,
      final differentialState=differentialState,
      final m_flow_init=m_flow_init,
      redeclare model HeatTransfer = HeatTransfer_1,
      redeclare model SlipRatio = SlipRatio,
      redeclare model FlowModel = FlowModel_1,
      final p_a_start=p_a_start,
      final p_b_start=p_b_start,
      final h_init=h_init,
      final useLumpedPressure=useLumpedPressure);

    // Tube wall
    DynamicVCC.Components.Pipes.MetalWall metalWall(
      final Ncell=Ncell,
      final C=fill(C_metalWall/Ncell, Ncell),
      final T_init=Tt_init);

  protected
    SI.HeatFlowRate Q_flow_1 "Refrigerant capacity";
    SI.HeatFlowRate Q_flow_2 "Secondary fluid capacity";
    SI.Mass charge "Refrigerant charge";
    SI.Energy energy "Energy conservation";
    SI.Pressure dp "Pressure drop across heat exchanger";

  equation
    connect(refFlow.port_a,port_a1);
    connect(refFlow.port_b,port_b1);
    connect(refFlow.heatPorts,metalWall.heatPorts_a);
    Q_flow_1=sum(refFlow.heatTransfer.Q_flows);
    charge=sum(refFlow.rho.*refFlow.fluidVolumes);
    energy=sum(refFlow.rho.*refFlow.u.*refFlow.fluidVolumes);
    dp=abs(port_a1.p-port_b1.p);

    annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(coordinateSystem(preserveAspectRatio=false)));
  end PartialHX;

  record FinTubeCoil "Geometry of finned tube coil"

    extends Modelica.Icons.Record;

    import Modelica.Constants.pi;
    parameter SI.Diameter d_o=1 "tube outter diameter";
    parameter SI.Diameter d_i=1 "tube inner diameter";
    parameter SI.SpecificHeatCapacity cp_tube=385;
    parameter SI.Density rho_tube=8900;
    parameter SI.SpecificHeatCapacity cp_fin=900;
    parameter SI.Density rho_fin=2700;
    parameter SI.Length L_tube=1 "Length of single tube";
    parameter Integer N_row=1 "Number of tube rows";
    parameter Integer N_prow=1 "Number of tubes per row";
    parameter Integer N_circuits=1 "Number of circuits";
    parameter SI.Length pf=1 "distance between fins";
    parameter SI.Length pt=1 "Tube transverse spacing";
    parameter SI.Length pl=1 "Tube longitudinal spacing";
    parameter SI.Thickness t_fin=1 "Fin thickness";
    parameter Real Eta_fin_overall=1;
    parameter SI.Area Ac_e=1 "air-side cross sectional area";
    parameter SI.Area HTA_e=1 "total air side heat transfer area";
    parameter SI.Mass M_fin=1 "fin mass";

    parameter SI.Area HTA_r=N_row*N_prow*pi*d_i*L_tube "tube inner surface area";
    parameter SI.Area Ac_r=N_circuits*pi*d_i^2/4 "tube cross-sectional area";
    parameter SI.Mass M_tube=N_row*N_prow*pi/4*(d_o^2-d_i^2)*L_tube*rho_tube "tube mass";
    parameter SI.Length L_circuit=N_row*N_prow*L_tube/N_circuits "average circuit length";
    parameter SI.Length L_fin=pt*N_row "fin depth";
    parameter SI.HeatCapacity C_FinTube=M_tube*cp_tube+M_fin*cp_fin;



    annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(coordinateSystem(preserveAspectRatio=false)));
  end FinTubeCoil;

  package AirPressureDrop "Air-side frictional pressure drop"
    extends Modelica.Icons.VariantsPackage;

    model PartialPressureDrop "Base class for air-side pressure drop"
      replaceable package Medium=Modelica.Media.Air.MoistAir;

      parameter Integer n=1;

      input Medium.ThermodynamicState states[n] "Fluid states";


      annotation (Icon(coordinateSystem(preserveAspectRatio=false), graphics={
                              Line(
              points={{-80,-60},{-80,60},{80,-60},{80,62}},
              color={0,0,255},
              thickness=1)}), Diagram(coordinateSystem(preserveAspectRatio=false)));
    end PartialPressureDrop;
  end AirPressureDrop;
end BaseClasses;
