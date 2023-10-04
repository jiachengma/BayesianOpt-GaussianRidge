within DynamicVCC.Components.Pipes.BaseClasses;
package WallFriction
  extends Modelica.Icons.VariantsPackage;
  import Modelica.Constants.pi;

  package PartialWallFriction
    replaceable partial function massFlowRate_dp

      extends Modelica.Icons.Function;
      input SI.Pressure dp "Pressure drop dp = port_a.p - port_b.p";
      input SI.Density rho;
      input SI.DynamicViscosity mu;
      input SI.Length length "Pipe length";
      input SI.Diameter diameter "Pipe diameter";
      input SI.Area crossArea=pi*diameter^2/4 "Cross-sectional area";
      input SI.AbsolutePressure dp_small
      "Regularization of zero flow if |dp| < dp_small";
      output SI.MassFlowRate m_flow;
      annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(coordinateSystem(preserveAspectRatio=false)));
    end massFlowRate_dp;

    replaceable partial function pressureLoss_m_flow
      extends Modelica.Icons.Function;
      input SI.MassFlowRate m_flow;
      input SI.Density rho;
      input SI.DynamicViscosity mu;
      input SI.Length length "Pipe length";
      input SI.Diameter diameter "Pipe diameter";
      input SI.Area crossArea=pi*diameter^2/4 "Cross-sectional area";
      input SI.MassFlowRate m_flow_small
      "Regularization neighborhood for zeros flow if |m_flow| < m_flow_small";
      output SI.Pressure dp "Frictional pressure drop (dp=port_a.p-port_b.p)";

    end pressureLoss_m_flow;
  end PartialWallFriction;

  package Correlation_SinglePhase "Friction factor correlations for single-phase flow"
    // Functions should be monotonically increasing so that a unique inverse exists
    extends DynamicVCC.Components.Pipes.BaseClasses.WallFriction.PartialWallFriction;

    redeclare function extends massFlowRate_dp
      import Modelica.Math;
      import DynamicVCC.Utilities.cubicPolynomial;
    protected
      parameter SI.ReynoldsNumber Re1=2000 "Re leaving laminar region";
      parameter SI.ReynoldsNumber Re2=6000 "Re entering turbulent region";
      parameter Real lambda21=64*Re1;
      parameter Real lambda22=0.3164*Re2^1.75;
      SI.ReynoldsNumber Re;
      Real lambda2 "Modified friction coefficient (= lambda*Re^2)";
      function interpolateTransition "Interpolate in the transition region between laminar and turbulent"
        input Real lambda2;
        input Real lambda21;
        input Real lambda22;
        output SI.ReynoldsNumber Re;
      protected
        Real x=Math.log10(lambda2);
        // Starting point in log() scale
        Real x1=Math.log10(lambda21);
        Real y1=Math.log10(lambda21/64);
        Real dy1=1;
        // End point in log() scale;
        Real x2=Math.log10(lambda22);
        Real y2=Math.log10((lambda22/0.3164)^(1/1.75));
        Real dy2=1/1.75;
      algorithm
        Re:=10^cubicPolynomial(
            x,
            x1,
            x2,
            y1,
            y2,
            dy1,
            dy2);
      end interpolateTransition;
    algorithm
      lambda2:=abs(dp)*2*rho*diameter^3/(length*mu^2);
      Re:=smooth(1,noEvent(if lambda2<lambda21 then
      lambda2/64 elseif lambda2>lambda22 then
      (lambda2/0.3164)^(1/1.75) else
      interpolateTransition(lambda2,lambda21,lambda22)));
      m_flow:=crossArea*mu/diameter*(if dp>=0 then Re else -Re);



      annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(coordinateSystem(preserveAspectRatio=false)));
    end massFlowRate_dp;

    redeclare function extends pressureLoss_m_flow
      import Modelica.Math;
      import DynamicVCC.Utilities.cubicPolynomial;
    protected
      parameter SI.ReynoldsNumber Re1=2000 "Re leaving laminar region";
      parameter SI.ReynoldsNumber Re2=6000 "Re entering turbulent region";
      SI.ReynoldsNumber Re;
      Real lambda2 "Modified friction coefficient (= lambda*Re^2)";

    function interpolateTransition "Interpolate in the transition region between laminar and turbulent"
      input SI.ReynoldsNumber Re;
      input SI.ReynoldsNumber Re1;
      input SI.ReynoldsNumber Re2;
      output Real lambda2;
      protected
      Real x=Math.log10(Re);
      // Starting point in log() scale
      Real x1=Math.log10(Re1);
      Real y1=Math.log10(64*Re1);
      Real dy1=1;
      // End point in log() scale
      Real x2=Math.log10(Re2);
      Real y2=Math.log10(0.3164*Re2^1.75);
      Real dy2=1.75;
    algorithm
      lambda2:=10^cubicPolynomial(
          x,
          x1,
          x2,
          y1,
          y2,
          dy1,
          dy2);
    end interpolateTransition;

    algorithm

      // Determine pressure drop
      Re:=Modelica.Fluid.Pipes.BaseClasses.CharacteristicNumbers.ReynoldsNumber_m_flow(
        m_flow,
        mu,
        diameter);
      // Cubic interpolation in transition region 2000<Re<6000
      lambda2:=smooth(1, noEvent(if Re < Re1 then 64*Re elseif
       Re > Re2 then 0.3164*Re^1.75 else
       interpolateTransition(Re,Re1,Re2)));
      dp:=length*mu^2/(2*rho*diameter^3)*(if m_flow>=0 then lambda2 else -lambda2);

      annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(coordinateSystem(preserveAspectRatio=false)));
    end pressureLoss_m_flow;
  end Correlation_SinglePhase;

  package Constant "Constant friction factor"
    extends DynamicVCC.Components.Pipes.BaseClasses.WallFriction.PartialWallFriction;


    redeclare function extends massFlowRate_dp
      import Modelica.Fluid.Utilities.regRoot;
      input Real lambda0;
    algorithm
      m_flow:=crossArea*sqrt(2*diameter*rho/length)*regRoot(dp)/sqrt(lambda0);
      annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(coordinateSystem(preserveAspectRatio=false)));
    end massFlowRate_dp;

    redeclare function extends pressureLoss_m_flow
      input Real lambda0;
    algorithm
      dp:=lambda0*length/diameter*m_flow^2/(2*rho*crossArea^2);
      annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(coordinateSystem(preserveAspectRatio=false)));
    end pressureLoss_m_flow;
  end Constant;

  package Nominal
    extends DynamicVCC.Components.Pipes.BaseClasses.WallFriction.PartialWallFriction;

    replaceable function extends massFlowRate_dp
      import DynamicVCC.Utilities.regPowGen;
      input SI.MassFlowRate m_flow_nominal;
      input SI.Pressure dp_nominal;
      input Real k "multiplier";
      input Real n "power constant";
    algorithm
      m_flow:=m_flow_nominal*regPowGen(dp/k/dp_nominal, n, 1);

      annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(coordinateSystem(preserveAspectRatio=false)));
    end massFlowRate_dp;

    replaceable function extends pressureLoss_m_flow
      import DynamicVCC.Utilities.regPowGen;
      input SI.MassFlowRate m_flow_nominal;
      input SI.Pressure dp_nominal;
      input Real k "multiplier";
      input Real n "power constant";
    algorithm
      dp:=k*dp_nominal*regPowGen(m_flow/m_flow_nominal,n,1e-4);
      annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(coordinateSystem(preserveAspectRatio=false)));
    end pressureLoss_m_flow;
  end Nominal;

  package Correlation_Condensation
    extends DynamicVCC.Components.Pipes.BaseClasses.WallFriction.PartialWallFriction;


    redeclare function extends massFlowRate_dp
      annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(coordinateSystem(preserveAspectRatio=false)));
    end massFlowRate_dp;

    redeclare function extends pressureLoss_m_flow
      annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(coordinateSystem(preserveAspectRatio=false)));
    end pressureLoss_m_flow;
  end Correlation_Condensation;
end WallFriction;
