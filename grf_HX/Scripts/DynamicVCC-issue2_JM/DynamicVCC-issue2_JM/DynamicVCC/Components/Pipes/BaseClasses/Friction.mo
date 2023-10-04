within DynamicVCC.Components.Pipes.BaseClasses;
package Friction "Friction terms in momentum balances"

  extends Modelica.Icons.VariantsPackage;
  partial model PartialFrictionalPressureDrop "Partial model of frictional pressure drop"

    replaceable package Medium=Modelica.Media.Interfaces.PartialTwoPhaseMedium;

    parameter Integer m=2 "Number of momentum balances";

    parameter SI.Pressure dp_nominal=1e3;

    parameter SI.MassFlowRate m_flow_nominal=1;

    // Inputs
    input Medium.Density rho[m];

    input Medium.DynamicViscosity mu[m];

    input Medium.SaturationProperties sat[m];

    input SI.Velocity v[m];

    input SI.Length L[m] "Flow length";

    input SI.Diameter dimension[m] "Tube diameter";

    output SI.Pressure dp[m](each start=dp_nominal,displayUnit="Pa");

  protected
    SI.ReynoldsNumber Re[m];
  equation
    Re=Modelica.Fluid.Pipes.BaseClasses.CharacteristicNumbers.ReynoldsNumber(v,rho,mu,dimension);

    annotation (Icon(coordinateSystem(preserveAspectRatio=false), graphics={Line(
            points={{-80,36}},
            color={42,22,200},
            thickness=0.5), Line(
            points={{-80,-62},{-80,58},{80,-62},{80,60}},
            color={0,0,255},
            thickness=1)}),                                        Diagram(coordinateSystem(preserveAspectRatio=false)));
  end PartialFrictionalPressureDrop;

  partial model PartialAirCoilDP "Partial firctional pressure drop of air coil"

    extends DynamicVCC.Components.Pipes.BaseClasses.Friction.PartialFrictionalPressureDrop;

    input SI.Area Ac[m];

    input SI.Area As[m];

  protected
    Real G[m] "Mass flux";
  equation
    G={v[i]*rho[i] for i in 1:m};

    annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(coordinateSystem(preserveAspectRatio=false)));
  end PartialAirCoilDP;

  model AirCoilDP_ConstFactor "Constant frictional factor"

    extends DynamicVCC.Components.Pipes.BaseClasses.Friction.PartialAirCoilDP;

    parameter Real f0=1 "Frictional factor";

  equation
    dp={f0*G[i]^2*As[i]/Ac[i]/rho[i]/2 for i in 1:m};

    annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(coordinateSystem(preserveAspectRatio=false)));
  end AirCoilDP_ConstFactor;

  model AirCoilDP_Wang "Air coil frictional pressure drop correlation (Wang, 2000)"

    import Modelica.Fluid.Utilities.regStep;

    extends DynamicVCC.Components.Pipes.BaseClasses.Friction.PartialAirCoilDP;

    // Detailed fin geometry
    parameter Integer Nt=1 "Number of tube rows";

    parameter SI.Diameter D=1 "Tube outter diameter";

    parameter SI.Length pf=1 "Fin pitch";

    parameter SI.Length pl=1 "Tube longitudinal tube pitch";

    parameter SI.Length pt=1 "Tube transverse tube pitch";

    parameter SI.Thickness t_fin=1 "Fin thickness";

  protected
    SI.Diameter Dh[m];
    SI.Diameter Dc;
    Real F1[m],F2[m],F3[m] "Correlation parameters";
    Real f[m] "Friction factor";
    SI.ReynoldsNumber Re_Dc[m];

  equation
    Dc=D+2*t_fin;
    Dh={4*Ac[i]*(pl*(Nt-1)+D)/As[i] for i in 1:m};
    Re_Dc=Modelica.Fluid.Pipes.BaseClasses.CharacteristicNumbers.ReynoldsNumber(v,rho,mu,Dc);

    for i in 1:m loop
      F1[i]=-0.764 + 0.739*pt/pl + 0.177*pf/Dc - 0.00758/Nt;
      F2[i]=-15.689 + 64.021/log(Re_Dc[i]);
      F3[i]=1.696 - 15.695/log(Re_Dc[i]);
      f[i]=regStep(Re_Dc[i]-50,0.0267*Re_Dc[i]^F1[i]*(pt/pl)^F2[i]*(pf/Dc)^F3[i],0.1,1);
      dp[i]=f[i]*G[i]^2*As[i]/Ac[i]/rho[i]/2;
    end for;

    annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(coordinateSystem(preserveAspectRatio=false)));
  end AirCoilDP_Wang;

  model AirCoilDP_Wang1998 "Wavy fin-and-tube heat exchanger"

    import Modelica.Fluid.Utilities.regStep;

    extends DynamicVCC.Components.Pipes.BaseClasses.Friction.PartialAirCoilDP;

    // Detailed fin geometry

    parameter SI.Diameter D=1 "Tube outter diameter";

    parameter SI.Length pf=1 "Fin pitch";

    parameter SI.Area At[m] "External tube surface area";

    parameter SI.Thickness t_fin=1 "Fin thickness";

  protected
    SI.Diameter Dc;
    Real f[m] "Friction factor";
    SI.ReynoldsNumber Re_Dc[m];

  equation
    Dc=D+2*t_fin;
    Re_Dc=Modelica.Fluid.Pipes.BaseClasses.CharacteristicNumbers.ReynoldsNumber(v,rho,mu,Dc);

    for i in 1:m loop
      f[i]=regStep(Re_Dc[i]-50,0.264*(0.105+0.708*exp(-Re_Dc[i]/225))*Re_Dc[i]^(-0.637)*(As[i]/At[i])^0.263*(pf/Dc)^(-0.317),0.1,1);
      dp[i]=f[i]*G[i]^2*As[i]/Ac[i]/rho[i]/2;
    end for;

    annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(coordinateSystem(preserveAspectRatio=false)),
                Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(coordinateSystem(preserveAspectRatio=false)));
  end AirCoilDP_Wang1998;

  package Correlations "Correlations for refrigerant frictional coefficient"
    extends Modelica.Icons.VariantsPackage;
    partial model PartialTwoPhase "Partial model of two-phase pressure drop multiplier"

      outer DynamicsVCC.Components.System system;

      extends PartialFrictionalPressureDrop;

      parameter SI.Acceleration g=system.g;

      parameter Boolean liquid=true;

      replaceable model LiquidPhase=SaturatedPressureDrop (
      final liquid=liquid);

      LiquidPhase liquidPhase(
      redeclare final package Medium=Medium,
      final m=m,
      final rho=rho,
      final mu=mu_l,
      final sat=sat,
      final v=v,
      final L=L,
      final dimension=dimension);

      Real x[m] "Vapor quality";
      Medium.ThermodynamicState bubbleStates[m];
      Medium.ThermodynamicState dewStates[m];
      Medium.Density rho_l[m]=Medium.density(bubbleStates);
      Medium.Density rho_v[m]=Medium.density(dewStates);
      Medium.DynamicViscosity mu_l[m]=Medium.dynamicViscosity(bubbleStates);
      Medium.DynamicViscosity mu_v[m]=Medium.dynamicViscosity(dewStates);

    protected
      Real Phi[m] "Two-phase multiplier";

    equation
      for i in 1:m loop
        bubbleStates[i]=Medium.setBubbleState(sat[i]);
        dewStates[i]=Medium.setDewState(sat[i]);
        x[i]=max(0.0,min(1.0,(1.0/rho[i]-1.0/rho_l[i])/(1.0/rho_v[i]-1.0/rho_l[i])));
      end for;

      annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(coordinateSystem(preserveAspectRatio=false)));
    end PartialTwoPhase;

    model SinglePhase "Frictional factor for single phase flow in pipes"

      extends PartialFrictionalPressureDrop;

      parameter SI.ReynoldsNumber Re1=2000 "Re leaving laminar region";
      parameter SI.ReynoldsNumber Re2=6000 "Re entering turbulent curve";

    protected
      Real f_modify[m] "Modified friction factor Re^2*f";

    equation

      // Cubic interpolation in transition region 2000<Re<6000
      for i in 1:m loop
        //f_modify[i]= if Re[i]<Re1 then 64*Re[i] elseif Re[i]>Re2 then 0.3164/Re[i]^0.25 else -8.85155347892110e-06*Re[i]^3+0.145402991488032*Re[i]^2-411.393324205076*Re[i]+439987.110289392;
          f_modify[i]=smooth(1, noEvent(if Re[i]<Re1 then 64*Re[i] elseif Re[i]>Re2 then 0.3164*Re[i]^1.75 else -8.85155347892110e-06*Re[i]^3+0.145402991488032*Re[i]^2-411.393324205076*Re[i]+439987.110289392));
      end for;
      dp={f_modify[i]*L[i]*mu[i]^2/dimension[i]^3/rho[i]/2 for i in 1:m};

    end SinglePhase;

    model Constant "Calculate frictional pressure drop from constant Darcy friction factor"
      extends PartialFrictionalPressureDrop;

      parameter Real f0=1;

    equation
      dp={f0*L[i]*v[i]^2*rho[i]/dimension[i]/2 for i in 1:m};

    end Constant;

    model SaturatedPressureDrop "Liquid-phase or vapor-phase frictional pressure drop"

      extends PartialFrictionalPressureDrop;

      parameter Boolean liquid=true;

    protected
      Real f_modify[m] "Modified friction factor Re^2*f";
      Medium.Density rho_l[m]=Medium.bubbleDensity(sat) "Bubble line density";
      Medium.Density rho_v[m]=Medium.dewDensity(sat) "Dew line density";
      Medium.Density rho_sat[m]=if liquid then rho_l else rho_v;
    equation

      for i in 1:m loop
        f_modify[i]=0.079*Re[i]^1.75;
        dp[i]=f_modify[i]*L[i]*mu[i]^2/dimension[i]^3/rho_sat[i]/2;
      end for;

      annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(coordinateSystem(preserveAspectRatio=false)));
    end SaturatedPressureDrop;

    model Boiling_Gronnerud "Gronnerud correlation (1979)"

      extends PartialTwoPhase(
      final liquid=true);

    protected
      Real Fr[m];
      Real f[m];
      Real dpdz[m];

    equation

      for i in 1:m loop
        Fr[i]=(rho[i]*v[i])^2/(g*dimension[i]*rho_l[i]^2);
        f[i]=smooth(0,noEvent(if Fr[i]<1.0 then Fr[i]^0.3+0.0055*(log(1/Fr[i]))^2 else 1.0));
        dpdz[i]=f[i]*(x[i]+4*(x[i]^1.8-x[i]^10*f[i]^0.5));
        Phi[i]=1+dpdz[i]*(rho_l[i]/rho_v[i]/(mu_l[i]/mu_v[i])^0.25-1);
        dp[i]=Phi[i]*liquidPhase.dp[i];
      end for;

      annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(coordinateSystem(preserveAspectRatio=false)));
    end Boiling_Gronnerud;

    model Condensing_LockhartMartinelli "Lockhart and Martinelli correlation (1949)"

      extends PartialTwoPhase(final liquid=true);

      replaceable model VaporPhase=SaturatedPressureDrop (
      final liquid=false);

      VaporPhase vaporPhase(
      redeclare final package Medium=Medium,
      final m=m,
      final rho=rho,
      final mu=mu_v,
      final sat=sat,
      final v=v,
      final L=L,
      final dimension=dimension);

    equation

      annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(coordinateSystem(preserveAspectRatio=false)));
    end Condensing_LockhartMartinelli;

    model TwoPhase_SinglePhase "Transition of single-phase and two-phase correlations"

      import Modelica.Media.Air.MoistAir.Utilities.spliceFunction;

      extends PartialFrictionalPressureDrop;

      replaceable model TwoPhase=Boiling_Gronnerud;

      TwoPhase twoPhase(
      redeclare final package Medium=Medium,
      final m=m,
      final rho=rho,
      final mu=mu,
      final sat=sat,
      final v=v,
      final L=L,
      final dimension=dimension);

      replaceable model SinglePhase=Correlations.SinglePhase;

      SinglePhase singlePhase(
      redeclare final package Medium=Medium,
      final m=m,
      final rho=rho,
      final mu=mu,
      final sat=sat,
      final v=v,
      final L=L,
      final dimension=dimension);

    protected
      SI.Pressure dp_l_tp[m];
      SI.Pressure dp_tp_v[m];
    equation

      for i in 1:m loop
        dp_l_tp[i]=spliceFunction(twoPhase.dp[i],singlePhase.dp[i],twoPhase.rho_l[i]-rho[i],1);
        dp_tp_v[i]=spliceFunction(twoPhase.dp[i],singlePhase.dp[i],rho[i]-twoPhase.rho_v[i],1);
        //dp[i]=smooth(0,noEvent(if rho[i]<(twoPhase.rho_l[i]+twoPhase.rho_v[i])/2 then dp_tp_v[i] else dp_l_tp[i]));
        dp[i]=if twoPhase.x[i]<0.0 then singlePhase.dp[i] elseif twoPhase.x[i]>1.0 then singlePhase.dp[i] else twoPhase.dp[i];
      end for;

      annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(coordinateSystem(preserveAspectRatio=false)));
    end TwoPhase_SinglePhase;
  end Correlations;
end Friction;
