#include <moutil.c>
PreNonAliasDef(6)
PreNonAliasDef(7)
PreNonAliasDef(8)
PreNonAliasDef(9)
PreNonAliasDef(10)
StartNonAlias(5)
DeclareAlias2("sink_OD[8].medium.T", "Temperature of medium [K|degC]", \
"sink_OD[8].T", 1, 7, 391, 0)
DeclareAlias2("sink_OD[8].medium.X[1]", "Mass fractions (= (component mass)/total mass  m_i/m) [kg/kg]",\
 "sink_OD[8].X[1]", 1, 7, 392, 0)
DeclareVariable("sink_OD[8].medium.X[2]", "Mass fractions (= (component mass)/total mass  m_i/m) [kg/kg]",\
 0.99, 0.0,1.0,0.1,0,513)
DeclareVariable("sink_OD[8].medium.u", "Specific internal energy of medium [J/kg]",\
 0.0, -100000000.0,100000000.0,1000000.0,0,513)
DeclareVariable("sink_OD[8].medium.R_s", "Gas constant (of mixture if applicable) [J/(kg.K)]",\
 1000.0, 0.0,10000000.0,1000.0,0,513)
DeclareVariable("sink_OD[8].medium.MM", "Molar mass (of mixture or single fluid) [kg/mol]",\
 0.032, 0.001,0.25,0.032,0,513)
DeclareAlias2("sink_OD[8].medium.state.p", "Absolute pressure of medium [Pa|bar]",\
 "sink_OD[8].p", 1, 7, 390, 0)
DeclareAlias2("sink_OD[8].medium.state.T", "Temperature of medium [K|degC]", \
"sink_OD[8].T", 1, 7, 391, 0)
DeclareAlias2("sink_OD[8].medium.state.X[1]", "Mass fractions (= (component mass)/total mass  m_i/m) [kg/kg]",\
 "sink_OD[8].X[1]", 1, 7, 392, 0)
DeclareAlias2("sink_OD[8].medium.state.X[2]", "Mass fractions (= (component mass)/total mass  m_i/m) [kg/kg]",\
 "sink_OD[8].medium.X[2]", 1, 5, 5654, 0)
DeclareVariable("sink_OD[8].medium.preferredMediumStates", "= true if StateSelect.prefer shall be used for the independent property variables of the medium [:#(type=Boolean)]",\
 false, 0.0,0.0,0.0,0,515)
DeclareVariable("sink_OD[8].medium.standardOrderComponents", "If true, and reducedX = true, the last element of X will be computed from the other ones [:#(type=Boolean)]",\
 true, 0.0,0.0,0.0,0,515)
DeclareVariable("sink_OD[8].medium.T_degC", "Temperature of medium in [degC] [degC;]",\
 0.0, 0.0,0.0,0.0,0,513)
DeclareVariable("sink_OD[8].medium.p_bar", "Absolute pressure of medium in [bar] [bar]",\
 0.0, 0.0,0.0,0.0,0,513)
DeclareVariable("sink_OD[8].medium.x_water", "Mass of total water/mass of dry air [kg/kg]",\
 0.0, 0.0,1.0,0.1,0,513)
DeclareVariable("sink_OD[8].medium.phi", "Relative humidity", 0.0, 0.0,0.0,0.0,0,513)
DeclareVariable("sink_OD[8].medium.X_liquid", "Mass fraction of liquid or solid water [kg/kg]",\
 0.0, 0.0,1.0,0.1,0,2561)
DeclareVariable("sink_OD[8].medium.X_steam", "Mass fraction of steam water [kg/kg]",\
 0.0, 0.0,1.0,0.1,0,2561)
DeclareVariable("sink_OD[8].medium.X_air", "Mass fraction of air [kg/kg]", 0.0, \
0.0,1.0,0.1,0,2561)
DeclareVariable("sink_OD[8].medium.X_sat", "Steam water mass fraction of saturation boundary in kg_water/kg_moistair [kg/kg]",\
 0.0, 0.0,1.0,0.1,0,2561)
DeclareVariable("sink_OD[8].medium.x_sat", "Steam water mass content of saturation boundary in kg_water/kg_dryair [kg/kg]",\
 0.0, 0.0,1.0,0.1,0,2561)
DeclareVariable("sink_OD[8].medium.p_steam_sat", "Partial saturation pressure of steam [Pa|bar]",\
 100000.0, 0.0,100000000.0,100000.0,0,2561)
DeclareAlias2("sink_OD[8].ports[1].m_flow", "Mass flow rate from the connection point into the component [kg/s]",\
 "fan_OD[8].ports[1].m_flow", -1, 5, 5204, 132)
DeclareAlias2("sink_OD[8].ports[1].p", "Thermodynamic pressure in the connection point [Pa|bar]",\
 "sink_OD[8].p", 1, 7, 390, 4)
DeclareVariable("sink_OD[8].ports[1].h_outflow", "Specific thermodynamic enthalpy close to the connection point if m_flow < 0 [J/kg]",\
 0.0, -10000000000.0,10000000000.0,1000000.0,0,521)
DeclareAlias2("sink_OD[8].ports[1].Xi_outflow[1]", "Independent mixture mass fractions m_i/m close to the connection point if m_flow < 0 [kg/kg]",\
 "sink_OD[8].X[1]", 1, 7, 392, 4)
DeclareVariable("sink_OD[8].flowDirection", "Allowed flow direction [:#(type=Modelica.Fluid.Types.PortFlowDirection)]",\
 3, 1.0,3.0,0.0,0,2565)
DeclareVariable("sink_OD[8].use_p_in", "Get the pressure from the input connector [:#(type=Boolean)]",\
 false, 0.0,0.0,0.0,0,1539)
DeclareVariable("sink_OD[8].use_T_in", "Get the temperature from the input connector [:#(type=Boolean)]",\
 false, 0.0,0.0,0.0,0,1539)
DeclareVariable("sink_OD[8].use_X_in", "Get the composition from the input connector [:#(type=Boolean)]",\
 false, 0.0,0.0,0.0,0,1539)
DeclareVariable("sink_OD[8].use_C_in", "Get the trace substances from the input connector [:#(type=Boolean)]",\
 false, 0.0,0.0,0.0,0,1539)
DeclareParameter("sink_OD[8].p", "Fixed value of pressure [Pa|bar]", 390, 101325,\
 0.0,100000000.0,100000.0,0,560)
DeclareParameter("sink_OD[8].T", "Fixed value of temperature [K|degC]", 391, \
293.15, 190.0,647.0,300.0,0,560)
DeclareParameter("sink_OD[8].X[1]", "Fixed value of composition [kg/kg]", 392, \
0.01, 0.0,1.0,0.1,0,560)
DeclareParameter("sink_OD[8].X[2]", "Fixed value of composition [kg/kg]", 393, \
0.99, 0.0,1.0,0.1,0,560)
DeclareAlias2("sink_OD[8].p_in_internal", "Needed to connect to conditional connector [Pa]",\
 "sink_OD[8].p", 1, 7, 390, 1024)
DeclareAlias2("sink_OD[8].T_in_internal", "Needed to connect to conditional connector [K]",\
 "sink_OD[8].T", 1, 7, 391, 1024)
DeclareVariable("sink_OD[8].X_in_internal[1]", "Needed to connect to conditional connector [1]",\
 0.0, 0.0,0.0,0.0,0,2561)
DeclareVariable("sink_OD[8].X_in_internal[2]", "Needed to connect to conditional connector [1]",\
 0.0, 0.0,0.0,0.0,0,2561)
DeclareVariable("sink_OD[9].nPorts", "Number of ports [:#(type=Integer)]", 1, \
0.0,0.0,0.0,0,517)
DeclareAlias2("sink_OD[9].medium.p", "Absolute pressure of medium [Pa|bar]", \
"sink_OD[9].p", 1, 7, 394, 0)
DeclareVariable("sink_OD[9].medium.Xi[1]", "Structurally independent mass fractions [1]",\
 0.01, 0.0,1.0,0.0,0,513)
DeclareAlias2("sink_OD[9].medium.h", "Specific enthalpy of medium [J/kg]", \
"sink_OD[9].ports[1].h_outflow", 1, 5, 5697, 0)
DeclareVariable("sink_OD[9].medium.d", "Density of medium [kg/m3|g/cm3]", 1, 0.0,\
100000.0,1.0,0,513)
DeclareAlias2("sink_OD[9].medium.T", "Temperature of medium [K|degC]", \
"sink_OD[9].T", 1, 7, 395, 0)
DeclareAlias2("sink_OD[9].medium.X[1]", "Mass fractions (= (component mass)/total mass  m_i/m) [kg/kg]",\
 "sink_OD[9].X[1]", 1, 7, 396, 0)
DeclareVariable("sink_OD[9].medium.X[2]", "Mass fractions (= (component mass)/total mass  m_i/m) [kg/kg]",\
 0.99, 0.0,1.0,0.1,0,513)
DeclareVariable("sink_OD[9].medium.u", "Specific internal energy of medium [J/kg]",\
 0.0, -100000000.0,100000000.0,1000000.0,0,513)
DeclareVariable("sink_OD[9].medium.R_s", "Gas constant (of mixture if applicable) [J/(kg.K)]",\
 1000.0, 0.0,10000000.0,1000.0,0,513)
DeclareVariable("sink_OD[9].medium.MM", "Molar mass (of mixture or single fluid) [kg/mol]",\
 0.032, 0.001,0.25,0.032,0,513)
DeclareAlias2("sink_OD[9].medium.state.p", "Absolute pressure of medium [Pa|bar]",\
 "sink_OD[9].p", 1, 7, 394, 0)
DeclareAlias2("sink_OD[9].medium.state.T", "Temperature of medium [K|degC]", \
"sink_OD[9].T", 1, 7, 395, 0)
DeclareAlias2("sink_OD[9].medium.state.X[1]", "Mass fractions (= (component mass)/total mass  m_i/m) [kg/kg]",\
 "sink_OD[9].X[1]", 1, 7, 396, 0)
DeclareAlias2("sink_OD[9].medium.state.X[2]", "Mass fractions (= (component mass)/total mass  m_i/m) [kg/kg]",\
 "sink_OD[9].medium.X[2]", 1, 5, 5681, 0)
DeclareVariable("sink_OD[9].medium.preferredMediumStates", "= true if StateSelect.prefer shall be used for the independent property variables of the medium [:#(type=Boolean)]",\
 false, 0.0,0.0,0.0,0,515)
DeclareVariable("sink_OD[9].medium.standardOrderComponents", "If true, and reducedX = true, the last element of X will be computed from the other ones [:#(type=Boolean)]",\
 true, 0.0,0.0,0.0,0,515)
DeclareVariable("sink_OD[9].medium.T_degC", "Temperature of medium in [degC] [degC;]",\
 0.0, 0.0,0.0,0.0,0,513)
DeclareVariable("sink_OD[9].medium.p_bar", "Absolute pressure of medium in [bar] [bar]",\
 0.0, 0.0,0.0,0.0,0,513)
DeclareVariable("sink_OD[9].medium.x_water", "Mass of total water/mass of dry air [kg/kg]",\
 0.0, 0.0,1.0,0.1,0,513)
DeclareVariable("sink_OD[9].medium.phi", "Relative humidity", 0.0, 0.0,0.0,0.0,0,513)
DeclareVariable("sink_OD[9].medium.X_liquid", "Mass fraction of liquid or solid water [kg/kg]",\
 0.0, 0.0,1.0,0.1,0,2561)
DeclareVariable("sink_OD[9].medium.X_steam", "Mass fraction of steam water [kg/kg]",\
 0.0, 0.0,1.0,0.1,0,2561)
DeclareVariable("sink_OD[9].medium.X_air", "Mass fraction of air [kg/kg]", 0.0, \
0.0,1.0,0.1,0,2561)
DeclareVariable("sink_OD[9].medium.X_sat", "Steam water mass fraction of saturation boundary in kg_water/kg_moistair [kg/kg]",\
 0.0, 0.0,1.0,0.1,0,2561)
DeclareVariable("sink_OD[9].medium.x_sat", "Steam water mass content of saturation boundary in kg_water/kg_dryair [kg/kg]",\
 0.0, 0.0,1.0,0.1,0,2561)
DeclareVariable("sink_OD[9].medium.p_steam_sat", "Partial saturation pressure of steam [Pa|bar]",\
 100000.0, 0.0,100000000.0,100000.0,0,2561)
DeclareAlias2("sink_OD[9].ports[1].m_flow", "Mass flow rate from the connection point into the component [kg/s]",\
 "fan_OD[9].ports[1].m_flow", -1, 5, 5239, 132)
DeclareAlias2("sink_OD[9].ports[1].p", "Thermodynamic pressure in the connection point [Pa|bar]",\
 "sink_OD[9].p", 1, 7, 394, 4)
DeclareVariable("sink_OD[9].ports[1].h_outflow", "Specific thermodynamic enthalpy close to the connection point if m_flow < 0 [J/kg]",\
 0.0, -10000000000.0,10000000000.0,1000000.0,0,521)
DeclareAlias2("sink_OD[9].ports[1].Xi_outflow[1]", "Independent mixture mass fractions m_i/m close to the connection point if m_flow < 0 [kg/kg]",\
 "sink_OD[9].X[1]", 1, 7, 396, 4)
DeclareVariable("sink_OD[9].flowDirection", "Allowed flow direction [:#(type=Modelica.Fluid.Types.PortFlowDirection)]",\
 3, 1.0,3.0,0.0,0,2565)
DeclareVariable("sink_OD[9].use_p_in", "Get the pressure from the input connector [:#(type=Boolean)]",\
 false, 0.0,0.0,0.0,0,1539)
DeclareVariable("sink_OD[9].use_T_in", "Get the temperature from the input connector [:#(type=Boolean)]",\
 false, 0.0,0.0,0.0,0,1539)
DeclareVariable("sink_OD[9].use_X_in", "Get the composition from the input connector [:#(type=Boolean)]",\
 false, 0.0,0.0,0.0,0,1539)
DeclareVariable("sink_OD[9].use_C_in", "Get the trace substances from the input connector [:#(type=Boolean)]",\
 false, 0.0,0.0,0.0,0,1539)
DeclareParameter("sink_OD[9].p", "Fixed value of pressure [Pa|bar]", 394, 101325,\
 0.0,100000000.0,100000.0,0,560)
DeclareParameter("sink_OD[9].T", "Fixed value of temperature [K|degC]", 395, \
293.15, 190.0,647.0,300.0,0,560)
DeclareParameter("sink_OD[9].X[1]", "Fixed value of composition [kg/kg]", 396, \
0.01, 0.0,1.0,0.1,0,560)
DeclareParameter("sink_OD[9].X[2]", "Fixed value of composition [kg/kg]", 397, \
0.99, 0.0,1.0,0.1,0,560)
DeclareAlias2("sink_OD[9].p_in_internal", "Needed to connect to conditional connector [Pa]",\
 "sink_OD[9].p", 1, 7, 394, 1024)
DeclareAlias2("sink_OD[9].T_in_internal", "Needed to connect to conditional connector [K]",\
 "sink_OD[9].T", 1, 7, 395, 1024)
DeclareVariable("sink_OD[9].X_in_internal[1]", "Needed to connect to conditional connector [1]",\
 0.0, 0.0,0.0,0.0,0,2561)
DeclareVariable("sink_OD[9].X_in_internal[2]", "Needed to connect to conditional connector [1]",\
 0.0, 0.0,0.0,0.0,0,2561)
DeclareVariable("sink_OD[10].nPorts", "Number of ports [:#(type=Integer)]", 1, \
0.0,0.0,0.0,0,517)
DeclareAlias2("sink_OD[10].medium.p", "Absolute pressure of medium [Pa|bar]", \
"sink_OD[10].p", 1, 7, 398, 0)
DeclareVariable("sink_OD[10].medium.Xi[1]", "Structurally independent mass fractions [1]",\
 0.01, 0.0,1.0,0.0,0,513)
DeclareAlias2("sink_OD[10].medium.h", "Specific enthalpy of medium [J/kg]", \
"sink_OD[10].ports[1].h_outflow", 1, 5, 5724, 0)
DeclareVariable("sink_OD[10].medium.d", "Density of medium [kg/m3|g/cm3]", 1, \
0.0,100000.0,1.0,0,513)
DeclareAlias2("sink_OD[10].medium.T", "Temperature of medium [K|degC]", \
"sink_OD[10].T", 1, 7, 399, 0)
DeclareAlias2("sink_OD[10].medium.X[1]", "Mass fractions (= (component mass)/total mass  m_i/m) [kg/kg]",\
 "sink_OD[10].X[1]", 1, 7, 400, 0)
DeclareVariable("sink_OD[10].medium.X[2]", "Mass fractions (= (component mass)/total mass  m_i/m) [kg/kg]",\
 0.99, 0.0,1.0,0.1,0,513)
DeclareVariable("sink_OD[10].medium.u", "Specific internal energy of medium [J/kg]",\
 0.0, -100000000.0,100000000.0,1000000.0,0,513)
DeclareVariable("sink_OD[10].medium.R_s", "Gas constant (of mixture if applicable) [J/(kg.K)]",\
 1000.0, 0.0,10000000.0,1000.0,0,513)
DeclareVariable("sink_OD[10].medium.MM", "Molar mass (of mixture or single fluid) [kg/mol]",\
 0.032, 0.001,0.25,0.032,0,513)
DeclareAlias2("sink_OD[10].medium.state.p", "Absolute pressure of medium [Pa|bar]",\
 "sink_OD[10].p", 1, 7, 398, 0)
DeclareAlias2("sink_OD[10].medium.state.T", "Temperature of medium [K|degC]", \
"sink_OD[10].T", 1, 7, 399, 0)
DeclareAlias2("sink_OD[10].medium.state.X[1]", "Mass fractions (= (component mass)/total mass  m_i/m) [kg/kg]",\
 "sink_OD[10].X[1]", 1, 7, 400, 0)
DeclareAlias2("sink_OD[10].medium.state.X[2]", "Mass fractions (= (component mass)/total mass  m_i/m) [kg/kg]",\
 "sink_OD[10].medium.X[2]", 1, 5, 5708, 0)
DeclareVariable("sink_OD[10].medium.preferredMediumStates", "= true if StateSelect.prefer shall be used for the independent property variables of the medium [:#(type=Boolean)]",\
 false, 0.0,0.0,0.0,0,515)
DeclareVariable("sink_OD[10].medium.standardOrderComponents", "If true, and reducedX = true, the last element of X will be computed from the other ones [:#(type=Boolean)]",\
 true, 0.0,0.0,0.0,0,515)
DeclareVariable("sink_OD[10].medium.T_degC", "Temperature of medium in [degC] [degC;]",\
 0.0, 0.0,0.0,0.0,0,513)
DeclareVariable("sink_OD[10].medium.p_bar", "Absolute pressure of medium in [bar] [bar]",\
 0.0, 0.0,0.0,0.0,0,513)
DeclareVariable("sink_OD[10].medium.x_water", "Mass of total water/mass of dry air [kg/kg]",\
 0.0, 0.0,1.0,0.1,0,513)
DeclareVariable("sink_OD[10].medium.phi", "Relative humidity", 0.0, 0.0,0.0,0.0,\
0,513)
DeclareVariable("sink_OD[10].medium.X_liquid", "Mass fraction of liquid or solid water [kg/kg]",\
 0.0, 0.0,1.0,0.1,0,2561)
DeclareVariable("sink_OD[10].medium.X_steam", "Mass fraction of steam water [kg/kg]",\
 0.0, 0.0,1.0,0.1,0,2561)
DeclareVariable("sink_OD[10].medium.X_air", "Mass fraction of air [kg/kg]", 0.0,\
 0.0,1.0,0.1,0,2561)
DeclareVariable("sink_OD[10].medium.X_sat", "Steam water mass fraction of saturation boundary in kg_water/kg_moistair [kg/kg]",\
 0.0, 0.0,1.0,0.1,0,2561)
DeclareVariable("sink_OD[10].medium.x_sat", "Steam water mass content of saturation boundary in kg_water/kg_dryair [kg/kg]",\
 0.0, 0.0,1.0,0.1,0,2561)
DeclareVariable("sink_OD[10].medium.p_steam_sat", "Partial saturation pressure of steam [Pa|bar]",\
 100000.0, 0.0,100000000.0,100000.0,0,2561)
DeclareAlias2("sink_OD[10].ports[1].m_flow", "Mass flow rate from the connection point into the component [kg/s]",\
 "fan_OD[10].ports[1].m_flow", -1, 5, 5274, 132)
DeclareAlias2("sink_OD[10].ports[1].p", "Thermodynamic pressure in the connection point [Pa|bar]",\
 "sink_OD[10].p", 1, 7, 398, 4)
DeclareVariable("sink_OD[10].ports[1].h_outflow", "Specific thermodynamic enthalpy close to the connection point if m_flow < 0 [J/kg]",\
 0.0, -10000000000.0,10000000000.0,1000000.0,0,521)
DeclareAlias2("sink_OD[10].ports[1].Xi_outflow[1]", "Independent mixture mass fractions m_i/m close to the connection point if m_flow < 0 [kg/kg]",\
 "sink_OD[10].X[1]", 1, 7, 400, 4)
DeclareVariable("sink_OD[10].flowDirection", "Allowed flow direction [:#(type=Modelica.Fluid.Types.PortFlowDirection)]",\
 3, 1.0,3.0,0.0,0,2565)
DeclareVariable("sink_OD[10].use_p_in", "Get the pressure from the input connector [:#(type=Boolean)]",\
 false, 0.0,0.0,0.0,0,1539)
DeclareVariable("sink_OD[10].use_T_in", "Get the temperature from the input connector [:#(type=Boolean)]",\
 false, 0.0,0.0,0.0,0,1539)
DeclareVariable("sink_OD[10].use_X_in", "Get the composition from the input connector [:#(type=Boolean)]",\
 false, 0.0,0.0,0.0,0,1539)
DeclareVariable("sink_OD[10].use_C_in", "Get the trace substances from the input connector [:#(type=Boolean)]",\
 false, 0.0,0.0,0.0,0,1539)
DeclareParameter("sink_OD[10].p", "Fixed value of pressure [Pa|bar]", 398, 101325,\
 0.0,100000000.0,100000.0,0,560)
DeclareParameter("sink_OD[10].T", "Fixed value of temperature [K|degC]", 399, \
293.15, 190.0,647.0,300.0,0,560)
DeclareParameter("sink_OD[10].X[1]", "Fixed value of composition [kg/kg]", 400, \
0.01, 0.0,1.0,0.1,0,560)
DeclareParameter("sink_OD[10].X[2]", "Fixed value of composition [kg/kg]", 401, \
0.99, 0.0,1.0,0.1,0,560)
DeclareAlias2("sink_OD[10].p_in_internal", "Needed to connect to conditional connector [Pa]",\
 "sink_OD[10].p", 1, 7, 398, 1024)
DeclareAlias2("sink_OD[10].T_in_internal", "Needed to connect to conditional connector [K]",\
 "sink_OD[10].T", 1, 7, 399, 1024)
DeclareVariable("sink_OD[10].X_in_internal[1]", "Needed to connect to conditional connector [1]",\
 0.0, 0.0,0.0,0.0,0,2561)
DeclareVariable("sink_OD[10].X_in_internal[2]", "Needed to connect to conditional connector [1]",\
 0.0, 0.0,0.0,0.0,0,2561)
DeclareVariable("sink_OD[11].nPorts", "Number of ports [:#(type=Integer)]", 1, \
0.0,0.0,0.0,0,517)
DeclareAlias2("sink_OD[11].medium.p", "Absolute pressure of medium [Pa|bar]", \
"sink_OD[11].p", 1, 7, 402, 0)
DeclareVariable("sink_OD[11].medium.Xi[1]", "Structurally independent mass fractions [1]",\
 0.01, 0.0,1.0,0.0,0,513)
DeclareAlias2("sink_OD[11].medium.h", "Specific enthalpy of medium [J/kg]", \
"sink_OD[11].ports[1].h_outflow", 1, 5, 5751, 0)
DeclareVariable("sink_OD[11].medium.d", "Density of medium [kg/m3|g/cm3]", 1, \
0.0,100000.0,1.0,0,513)
DeclareAlias2("sink_OD[11].medium.T", "Temperature of medium [K|degC]", \
"sink_OD[11].T", 1, 7, 403, 0)
DeclareAlias2("sink_OD[11].medium.X[1]", "Mass fractions (= (component mass)/total mass  m_i/m) [kg/kg]",\
 "sink_OD[11].X[1]", 1, 7, 404, 0)
DeclareVariable("sink_OD[11].medium.X[2]", "Mass fractions (= (component mass)/total mass  m_i/m) [kg/kg]",\
 0.99, 0.0,1.0,0.1,0,513)
DeclareVariable("sink_OD[11].medium.u", "Specific internal energy of medium [J/kg]",\
 0.0, -100000000.0,100000000.0,1000000.0,0,513)
DeclareVariable("sink_OD[11].medium.R_s", "Gas constant (of mixture if applicable) [J/(kg.K)]",\
 1000.0, 0.0,10000000.0,1000.0,0,513)
DeclareVariable("sink_OD[11].medium.MM", "Molar mass (of mixture or single fluid) [kg/mol]",\
 0.032, 0.001,0.25,0.032,0,513)
DeclareAlias2("sink_OD[11].medium.state.p", "Absolute pressure of medium [Pa|bar]",\
 "sink_OD[11].p", 1, 7, 402, 0)
DeclareAlias2("sink_OD[11].medium.state.T", "Temperature of medium [K|degC]", \
"sink_OD[11].T", 1, 7, 403, 0)
DeclareAlias2("sink_OD[11].medium.state.X[1]", "Mass fractions (= (component mass)/total mass  m_i/m) [kg/kg]",\
 "sink_OD[11].X[1]", 1, 7, 404, 0)
DeclareAlias2("sink_OD[11].medium.state.X[2]", "Mass fractions (= (component mass)/total mass  m_i/m) [kg/kg]",\
 "sink_OD[11].medium.X[2]", 1, 5, 5735, 0)
DeclareVariable("sink_OD[11].medium.preferredMediumStates", "= true if StateSelect.prefer shall be used for the independent property variables of the medium [:#(type=Boolean)]",\
 false, 0.0,0.0,0.0,0,515)
DeclareVariable("sink_OD[11].medium.standardOrderComponents", "If true, and reducedX = true, the last element of X will be computed from the other ones [:#(type=Boolean)]",\
 true, 0.0,0.0,0.0,0,515)
DeclareVariable("sink_OD[11].medium.T_degC", "Temperature of medium in [degC] [degC;]",\
 0.0, 0.0,0.0,0.0,0,513)
DeclareVariable("sink_OD[11].medium.p_bar", "Absolute pressure of medium in [bar] [bar]",\
 0.0, 0.0,0.0,0.0,0,513)
DeclareVariable("sink_OD[11].medium.x_water", "Mass of total water/mass of dry air [kg/kg]",\
 0.0, 0.0,1.0,0.1,0,513)
DeclareVariable("sink_OD[11].medium.phi", "Relative humidity", 0.0, 0.0,0.0,0.0,\
0,513)
DeclareVariable("sink_OD[11].medium.X_liquid", "Mass fraction of liquid or solid water [kg/kg]",\
 0.0, 0.0,1.0,0.1,0,2561)
DeclareVariable("sink_OD[11].medium.X_steam", "Mass fraction of steam water [kg/kg]",\
 0.0, 0.0,1.0,0.1,0,2561)
DeclareVariable("sink_OD[11].medium.X_air", "Mass fraction of air [kg/kg]", 0.0,\
 0.0,1.0,0.1,0,2561)
DeclareVariable("sink_OD[11].medium.X_sat", "Steam water mass fraction of saturation boundary in kg_water/kg_moistair [kg/kg]",\
 0.0, 0.0,1.0,0.1,0,2561)
DeclareVariable("sink_OD[11].medium.x_sat", "Steam water mass content of saturation boundary in kg_water/kg_dryair [kg/kg]",\
 0.0, 0.0,1.0,0.1,0,2561)
DeclareVariable("sink_OD[11].medium.p_steam_sat", "Partial saturation pressure of steam [Pa|bar]",\
 100000.0, 0.0,100000000.0,100000.0,0,2561)
DeclareAlias2("sink_OD[11].ports[1].m_flow", "Mass flow rate from the connection point into the component [kg/s]",\
 "fan_OD[11].ports[1].m_flow", -1, 5, 5309, 132)
DeclareAlias2("sink_OD[11].ports[1].p", "Thermodynamic pressure in the connection point [Pa|bar]",\
 "sink_OD[11].p", 1, 7, 402, 4)
DeclareVariable("sink_OD[11].ports[1].h_outflow", "Specific thermodynamic enthalpy close to the connection point if m_flow < 0 [J/kg]",\
 0.0, -10000000000.0,10000000000.0,1000000.0,0,521)
DeclareAlias2("sink_OD[11].ports[1].Xi_outflow[1]", "Independent mixture mass fractions m_i/m close to the connection point if m_flow < 0 [kg/kg]",\
 "sink_OD[11].X[1]", 1, 7, 404, 4)
DeclareVariable("sink_OD[11].flowDirection", "Allowed flow direction [:#(type=Modelica.Fluid.Types.PortFlowDirection)]",\
 3, 1.0,3.0,0.0,0,2565)
DeclareVariable("sink_OD[11].use_p_in", "Get the pressure from the input connector [:#(type=Boolean)]",\
 false, 0.0,0.0,0.0,0,1539)
DeclareVariable("sink_OD[11].use_T_in", "Get the temperature from the input connector [:#(type=Boolean)]",\
 false, 0.0,0.0,0.0,0,1539)
DeclareVariable("sink_OD[11].use_X_in", "Get the composition from the input connector [:#(type=Boolean)]",\
 false, 0.0,0.0,0.0,0,1539)
DeclareVariable("sink_OD[11].use_C_in", "Get the trace substances from the input connector [:#(type=Boolean)]",\
 false, 0.0,0.0,0.0,0,1539)
DeclareParameter("sink_OD[11].p", "Fixed value of pressure [Pa|bar]", 402, 101325,\
 0.0,100000000.0,100000.0,0,560)
DeclareParameter("sink_OD[11].T", "Fixed value of temperature [K|degC]", 403, \
293.15, 190.0,647.0,300.0,0,560)
DeclareParameter("sink_OD[11].X[1]", "Fixed value of composition [kg/kg]", 404, \
0.01, 0.0,1.0,0.1,0,560)
DeclareParameter("sink_OD[11].X[2]", "Fixed value of composition [kg/kg]", 405, \
0.99, 0.0,1.0,0.1,0,560)
DeclareAlias2("sink_OD[11].p_in_internal", "Needed to connect to conditional connector [Pa]",\
 "sink_OD[11].p", 1, 7, 402, 1024)
DeclareAlias2("sink_OD[11].T_in_internal", "Needed to connect to conditional connector [K]",\
 "sink_OD[11].T", 1, 7, 403, 1024)
DeclareVariable("sink_OD[11].X_in_internal[1]", "Needed to connect to conditional connector [1]",\
 0.0, 0.0,0.0,0.0,0,2561)
DeclareVariable("sink_OD[11].X_in_internal[2]", "Needed to connect to conditional connector [1]",\
 0.0, 0.0,0.0,0.0,0,2561)
DeclareVariable("sink_OD[12].nPorts", "Number of ports [:#(type=Integer)]", 1, \
0.0,0.0,0.0,0,517)
DeclareAlias2("sink_OD[12].medium.p", "Absolute pressure of medium [Pa|bar]", \
"sink_OD[12].p", 1, 7, 406, 0)
DeclareVariable("sink_OD[12].medium.Xi[1]", "Structurally independent mass fractions [1]",\
 0.01, 0.0,1.0,0.0,0,513)
DeclareAlias2("sink_OD[12].medium.h", "Specific enthalpy of medium [J/kg]", \
"sink_OD[12].ports[1].h_outflow", 1, 5, 5778, 0)
DeclareVariable("sink_OD[12].medium.d", "Density of medium [kg/m3|g/cm3]", 1, \
0.0,100000.0,1.0,0,513)
DeclareAlias2("sink_OD[12].medium.T", "Temperature of medium [K|degC]", \
"sink_OD[12].T", 1, 7, 407, 0)
DeclareAlias2("sink_OD[12].medium.X[1]", "Mass fractions (= (component mass)/total mass  m_i/m) [kg/kg]",\
 "sink_OD[12].X[1]", 1, 7, 408, 0)
DeclareVariable("sink_OD[12].medium.X[2]", "Mass fractions (= (component mass)/total mass  m_i/m) [kg/kg]",\
 0.99, 0.0,1.0,0.1,0,513)
DeclareVariable("sink_OD[12].medium.u", "Specific internal energy of medium [J/kg]",\
 0.0, -100000000.0,100000000.0,1000000.0,0,513)
DeclareVariable("sink_OD[12].medium.R_s", "Gas constant (of mixture if applicable) [J/(kg.K)]",\
 1000.0, 0.0,10000000.0,1000.0,0,513)
DeclareVariable("sink_OD[12].medium.MM", "Molar mass (of mixture or single fluid) [kg/mol]",\
 0.032, 0.001,0.25,0.032,0,513)
DeclareAlias2("sink_OD[12].medium.state.p", "Absolute pressure of medium [Pa|bar]",\
 "sink_OD[12].p", 1, 7, 406, 0)
DeclareAlias2("sink_OD[12].medium.state.T", "Temperature of medium [K|degC]", \
"sink_OD[12].T", 1, 7, 407, 0)
DeclareAlias2("sink_OD[12].medium.state.X[1]", "Mass fractions (= (component mass)/total mass  m_i/m) [kg/kg]",\
 "sink_OD[12].X[1]", 1, 7, 408, 0)
DeclareAlias2("sink_OD[12].medium.state.X[2]", "Mass fractions (= (component mass)/total mass  m_i/m) [kg/kg]",\
 "sink_OD[12].medium.X[2]", 1, 5, 5762, 0)
DeclareVariable("sink_OD[12].medium.preferredMediumStates", "= true if StateSelect.prefer shall be used for the independent property variables of the medium [:#(type=Boolean)]",\
 false, 0.0,0.0,0.0,0,515)
DeclareVariable("sink_OD[12].medium.standardOrderComponents", "If true, and reducedX = true, the last element of X will be computed from the other ones [:#(type=Boolean)]",\
 true, 0.0,0.0,0.0,0,515)
DeclareVariable("sink_OD[12].medium.T_degC", "Temperature of medium in [degC] [degC;]",\
 0.0, 0.0,0.0,0.0,0,513)
DeclareVariable("sink_OD[12].medium.p_bar", "Absolute pressure of medium in [bar] [bar]",\
 0.0, 0.0,0.0,0.0,0,513)
DeclareVariable("sink_OD[12].medium.x_water", "Mass of total water/mass of dry air [kg/kg]",\
 0.0, 0.0,1.0,0.1,0,513)
DeclareVariable("sink_OD[12].medium.phi", "Relative humidity", 0.0, 0.0,0.0,0.0,\
0,513)
DeclareVariable("sink_OD[12].medium.X_liquid", "Mass fraction of liquid or solid water [kg/kg]",\
 0.0, 0.0,1.0,0.1,0,2561)
DeclareVariable("sink_OD[12].medium.X_steam", "Mass fraction of steam water [kg/kg]",\
 0.0, 0.0,1.0,0.1,0,2561)
DeclareVariable("sink_OD[12].medium.X_air", "Mass fraction of air [kg/kg]", 0.0,\
 0.0,1.0,0.1,0,2561)
DeclareVariable("sink_OD[12].medium.X_sat", "Steam water mass fraction of saturation boundary in kg_water/kg_moistair [kg/kg]",\
 0.0, 0.0,1.0,0.1,0,2561)
DeclareVariable("sink_OD[12].medium.x_sat", "Steam water mass content of saturation boundary in kg_water/kg_dryair [kg/kg]",\
 0.0, 0.0,1.0,0.1,0,2561)
DeclareVariable("sink_OD[12].medium.p_steam_sat", "Partial saturation pressure of steam [Pa|bar]",\
 100000.0, 0.0,100000000.0,100000.0,0,2561)
DeclareAlias2("sink_OD[12].ports[1].m_flow", "Mass flow rate from the connection point into the component [kg/s]",\
 "fan_OD[12].ports[1].m_flow", -1, 5, 5344, 132)
DeclareAlias2("sink_OD[12].ports[1].p", "Thermodynamic pressure in the connection point [Pa|bar]",\
 "sink_OD[12].p", 1, 7, 406, 4)
DeclareVariable("sink_OD[12].ports[1].h_outflow", "Specific thermodynamic enthalpy close to the connection point if m_flow < 0 [J/kg]",\
 0.0, -10000000000.0,10000000000.0,1000000.0,0,521)
DeclareAlias2("sink_OD[12].ports[1].Xi_outflow[1]", "Independent mixture mass fractions m_i/m close to the connection point if m_flow < 0 [kg/kg]",\
 "sink_OD[12].X[1]", 1, 7, 408, 4)
DeclareVariable("sink_OD[12].flowDirection", "Allowed flow direction [:#(type=Modelica.Fluid.Types.PortFlowDirection)]",\
 3, 1.0,3.0,0.0,0,2565)
DeclareVariable("sink_OD[12].use_p_in", "Get the pressure from the input connector [:#(type=Boolean)]",\
 false, 0.0,0.0,0.0,0,1539)
DeclareVariable("sink_OD[12].use_T_in", "Get the temperature from the input connector [:#(type=Boolean)]",\
 false, 0.0,0.0,0.0,0,1539)
DeclareVariable("sink_OD[12].use_X_in", "Get the composition from the input connector [:#(type=Boolean)]",\
 false, 0.0,0.0,0.0,0,1539)
DeclareVariable("sink_OD[12].use_C_in", "Get the trace substances from the input connector [:#(type=Boolean)]",\
 false, 0.0,0.0,0.0,0,1539)
DeclareParameter("sink_OD[12].p", "Fixed value of pressure [Pa|bar]", 406, 101325,\
 0.0,100000000.0,100000.0,0,560)
DeclareParameter("sink_OD[12].T", "Fixed value of temperature [K|degC]", 407, \
293.15, 190.0,647.0,300.0,0,560)
DeclareParameter("sink_OD[12].X[1]", "Fixed value of composition [kg/kg]", 408, \
0.01, 0.0,1.0,0.1,0,560)
DeclareParameter("sink_OD[12].X[2]", "Fixed value of composition [kg/kg]", 409, \
0.99, 0.0,1.0,0.1,0,560)
DeclareAlias2("sink_OD[12].p_in_internal", "Needed to connect to conditional connector [Pa]",\
 "sink_OD[12].p", 1, 7, 406, 1024)
DeclareAlias2("sink_OD[12].T_in_internal", "Needed to connect to conditional connector [K]",\
 "sink_OD[12].T", 1, 7, 407, 1024)
DeclareVariable("sink_OD[12].X_in_internal[1]", "Needed to connect to conditional connector [1]",\
 0.0, 0.0,0.0,0.0,0,2561)
DeclareVariable("sink_OD[12].X_in_internal[2]", "Needed to connect to conditional connector [1]",\
 0.0, 0.0,0.0,0.0,0,2561)
DeclareVariable("sink_OD[13].nPorts", "Number of ports [:#(type=Integer)]", 1, \
0.0,0.0,0.0,0,517)
DeclareAlias2("sink_OD[13].medium.p", "Absolute pressure of medium [Pa|bar]", \
"sink_OD[13].p", 1, 7, 410, 0)
DeclareVariable("sink_OD[13].medium.Xi[1]", "Structurally independent mass fractions [1]",\
 0.01, 0.0,1.0,0.0,0,513)
DeclareAlias2("sink_OD[13].medium.h", "Specific enthalpy of medium [J/kg]", \
"sink_OD[13].ports[1].h_outflow", 1, 5, 5805, 0)
DeclareVariable("sink_OD[13].medium.d", "Density of medium [kg/m3|g/cm3]", 1, \
0.0,100000.0,1.0,0,513)
DeclareAlias2("sink_OD[13].medium.T", "Temperature of medium [K|degC]", \
"sink_OD[13].T", 1, 7, 411, 0)
DeclareAlias2("sink_OD[13].medium.X[1]", "Mass fractions (= (component mass)/total mass  m_i/m) [kg/kg]",\
 "sink_OD[13].X[1]", 1, 7, 412, 0)
DeclareVariable("sink_OD[13].medium.X[2]", "Mass fractions (= (component mass)/total mass  m_i/m) [kg/kg]",\
 0.99, 0.0,1.0,0.1,0,513)
DeclareVariable("sink_OD[13].medium.u", "Specific internal energy of medium [J/kg]",\
 0.0, -100000000.0,100000000.0,1000000.0,0,513)
DeclareVariable("sink_OD[13].medium.R_s", "Gas constant (of mixture if applicable) [J/(kg.K)]",\
 1000.0, 0.0,10000000.0,1000.0,0,513)
DeclareVariable("sink_OD[13].medium.MM", "Molar mass (of mixture or single fluid) [kg/mol]",\
 0.032, 0.001,0.25,0.032,0,513)
DeclareAlias2("sink_OD[13].medium.state.p", "Absolute pressure of medium [Pa|bar]",\
 "sink_OD[13].p", 1, 7, 410, 0)
DeclareAlias2("sink_OD[13].medium.state.T", "Temperature of medium [K|degC]", \
"sink_OD[13].T", 1, 7, 411, 0)
DeclareAlias2("sink_OD[13].medium.state.X[1]", "Mass fractions (= (component mass)/total mass  m_i/m) [kg/kg]",\
 "sink_OD[13].X[1]", 1, 7, 412, 0)
DeclareAlias2("sink_OD[13].medium.state.X[2]", "Mass fractions (= (component mass)/total mass  m_i/m) [kg/kg]",\
 "sink_OD[13].medium.X[2]", 1, 5, 5789, 0)
DeclareVariable("sink_OD[13].medium.preferredMediumStates", "= true if StateSelect.prefer shall be used for the independent property variables of the medium [:#(type=Boolean)]",\
 false, 0.0,0.0,0.0,0,515)
DeclareVariable("sink_OD[13].medium.standardOrderComponents", "If true, and reducedX = true, the last element of X will be computed from the other ones [:#(type=Boolean)]",\
 true, 0.0,0.0,0.0,0,515)
DeclareVariable("sink_OD[13].medium.T_degC", "Temperature of medium in [degC] [degC;]",\
 0.0, 0.0,0.0,0.0,0,513)
DeclareVariable("sink_OD[13].medium.p_bar", "Absolute pressure of medium in [bar] [bar]",\
 0.0, 0.0,0.0,0.0,0,513)
DeclareVariable("sink_OD[13].medium.x_water", "Mass of total water/mass of dry air [kg/kg]",\
 0.0, 0.0,1.0,0.1,0,513)
DeclareVariable("sink_OD[13].medium.phi", "Relative humidity", 0.0, 0.0,0.0,0.0,\
0,513)
DeclareVariable("sink_OD[13].medium.X_liquid", "Mass fraction of liquid or solid water [kg/kg]",\
 0.0, 0.0,1.0,0.1,0,2561)
DeclareVariable("sink_OD[13].medium.X_steam", "Mass fraction of steam water [kg/kg]",\
 0.0, 0.0,1.0,0.1,0,2561)
DeclareVariable("sink_OD[13].medium.X_air", "Mass fraction of air [kg/kg]", 0.0,\
 0.0,1.0,0.1,0,2561)
DeclareVariable("sink_OD[13].medium.X_sat", "Steam water mass fraction of saturation boundary in kg_water/kg_moistair [kg/kg]",\
 0.0, 0.0,1.0,0.1,0,2561)
DeclareVariable("sink_OD[13].medium.x_sat", "Steam water mass content of saturation boundary in kg_water/kg_dryair [kg/kg]",\
 0.0, 0.0,1.0,0.1,0,2561)
DeclareVariable("sink_OD[13].medium.p_steam_sat", "Partial saturation pressure of steam [Pa|bar]",\
 100000.0, 0.0,100000000.0,100000.0,0,2561)
DeclareAlias2("sink_OD[13].ports[1].m_flow", "Mass flow rate from the connection point into the component [kg/s]",\
 "fan_OD[13].ports[1].m_flow", -1, 5, 5379, 132)
DeclareAlias2("sink_OD[13].ports[1].p", "Thermodynamic pressure in the connection point [Pa|bar]",\
 "sink_OD[13].p", 1, 7, 410, 4)
DeclareVariable("sink_OD[13].ports[1].h_outflow", "Specific thermodynamic enthalpy close to the connection point if m_flow < 0 [J/kg]",\
 0.0, -10000000000.0,10000000000.0,1000000.0,0,521)
DeclareAlias2("sink_OD[13].ports[1].Xi_outflow[1]", "Independent mixture mass fractions m_i/m close to the connection point if m_flow < 0 [kg/kg]",\
 "sink_OD[13].X[1]", 1, 7, 412, 4)
DeclareVariable("sink_OD[13].flowDirection", "Allowed flow direction [:#(type=Modelica.Fluid.Types.PortFlowDirection)]",\
 3, 1.0,3.0,0.0,0,2565)
DeclareVariable("sink_OD[13].use_p_in", "Get the pressure from the input connector [:#(type=Boolean)]",\
 false, 0.0,0.0,0.0,0,1539)
DeclareVariable("sink_OD[13].use_T_in", "Get the temperature from the input connector [:#(type=Boolean)]",\
 false, 0.0,0.0,0.0,0,1539)
DeclareVariable("sink_OD[13].use_X_in", "Get the composition from the input connector [:#(type=Boolean)]",\
 false, 0.0,0.0,0.0,0,1539)
DeclareVariable("sink_OD[13].use_C_in", "Get the trace substances from the input connector [:#(type=Boolean)]",\
 false, 0.0,0.0,0.0,0,1539)
DeclareParameter("sink_OD[13].p", "Fixed value of pressure [Pa|bar]", 410, 101325,\
 0.0,100000000.0,100000.0,0,560)
DeclareParameter("sink_OD[13].T", "Fixed value of temperature [K|degC]", 411, \
293.15, 190.0,647.0,300.0,0,560)
DeclareParameter("sink_OD[13].X[1]", "Fixed value of composition [kg/kg]", 412, \
0.01, 0.0,1.0,0.1,0,560)
DeclareParameter("sink_OD[13].X[2]", "Fixed value of composition [kg/kg]", 413, \
0.99, 0.0,1.0,0.1,0,560)
DeclareAlias2("sink_OD[13].p_in_internal", "Needed to connect to conditional connector [Pa]",\
 "sink_OD[13].p", 1, 7, 410, 1024)
DeclareAlias2("sink_OD[13].T_in_internal", "Needed to connect to conditional connector [K]",\
 "sink_OD[13].T", 1, 7, 411, 1024)
DeclareVariable("sink_OD[13].X_in_internal[1]", "Needed to connect to conditional connector [1]",\
 0.0, 0.0,0.0,0.0,0,2561)
DeclareVariable("sink_OD[13].X_in_internal[2]", "Needed to connect to conditional connector [1]",\
 0.0, 0.0,0.0,0.0,0,2561)
DeclareVariable("sink_OD[14].nPorts", "Number of ports [:#(type=Integer)]", 1, \
0.0,0.0,0.0,0,517)
DeclareAlias2("sink_OD[14].medium.p", "Absolute pressure of medium [Pa|bar]", \
"sink_OD[14].p", 1, 7, 414, 0)
DeclareVariable("sink_OD[14].medium.Xi[1]", "Structurally independent mass fractions [1]",\
 0.01, 0.0,1.0,0.0,0,513)
DeclareAlias2("sink_OD[14].medium.h", "Specific enthalpy of medium [J/kg]", \
"sink_OD[14].ports[1].h_outflow", 1, 5, 5832, 0)
DeclareVariable("sink_OD[14].medium.d", "Density of medium [kg/m3|g/cm3]", 1, \
0.0,100000.0,1.0,0,513)
DeclareAlias2("sink_OD[14].medium.T", "Temperature of medium [K|degC]", \
"sink_OD[14].T", 1, 7, 415, 0)
DeclareAlias2("sink_OD[14].medium.X[1]", "Mass fractions (= (component mass)/total mass  m_i/m) [kg/kg]",\
 "sink_OD[14].X[1]", 1, 7, 416, 0)
DeclareVariable("sink_OD[14].medium.X[2]", "Mass fractions (= (component mass)/total mass  m_i/m) [kg/kg]",\
 0.99, 0.0,1.0,0.1,0,513)
DeclareVariable("sink_OD[14].medium.u", "Specific internal energy of medium [J/kg]",\
 0.0, -100000000.0,100000000.0,1000000.0,0,513)
DeclareVariable("sink_OD[14].medium.R_s", "Gas constant (of mixture if applicable) [J/(kg.K)]",\
 1000.0, 0.0,10000000.0,1000.0,0,513)
DeclareVariable("sink_OD[14].medium.MM", "Molar mass (of mixture or single fluid) [kg/mol]",\
 0.032, 0.001,0.25,0.032,0,513)
DeclareAlias2("sink_OD[14].medium.state.p", "Absolute pressure of medium [Pa|bar]",\
 "sink_OD[14].p", 1, 7, 414, 0)
DeclareAlias2("sink_OD[14].medium.state.T", "Temperature of medium [K|degC]", \
"sink_OD[14].T", 1, 7, 415, 0)
DeclareAlias2("sink_OD[14].medium.state.X[1]", "Mass fractions (= (component mass)/total mass  m_i/m) [kg/kg]",\
 "sink_OD[14].X[1]", 1, 7, 416, 0)
DeclareAlias2("sink_OD[14].medium.state.X[2]", "Mass fractions (= (component mass)/total mass  m_i/m) [kg/kg]",\
 "sink_OD[14].medium.X[2]", 1, 5, 5816, 0)
DeclareVariable("sink_OD[14].medium.preferredMediumStates", "= true if StateSelect.prefer shall be used for the independent property variables of the medium [:#(type=Boolean)]",\
 false, 0.0,0.0,0.0,0,515)
DeclareVariable("sink_OD[14].medium.standardOrderComponents", "If true, and reducedX = true, the last element of X will be computed from the other ones [:#(type=Boolean)]",\
 true, 0.0,0.0,0.0,0,515)
DeclareVariable("sink_OD[14].medium.T_degC", "Temperature of medium in [degC] [degC;]",\
 0.0, 0.0,0.0,0.0,0,513)
DeclareVariable("sink_OD[14].medium.p_bar", "Absolute pressure of medium in [bar] [bar]",\
 0.0, 0.0,0.0,0.0,0,513)
DeclareVariable("sink_OD[14].medium.x_water", "Mass of total water/mass of dry air [kg/kg]",\
 0.0, 0.0,1.0,0.1,0,513)
DeclareVariable("sink_OD[14].medium.phi", "Relative humidity", 0.0, 0.0,0.0,0.0,\
0,513)
DeclareVariable("sink_OD[14].medium.X_liquid", "Mass fraction of liquid or solid water [kg/kg]",\
 0.0, 0.0,1.0,0.1,0,2561)
DeclareVariable("sink_OD[14].medium.X_steam", "Mass fraction of steam water [kg/kg]",\
 0.0, 0.0,1.0,0.1,0,2561)
DeclareVariable("sink_OD[14].medium.X_air", "Mass fraction of air [kg/kg]", 0.0,\
 0.0,1.0,0.1,0,2561)
DeclareVariable("sink_OD[14].medium.X_sat", "Steam water mass fraction of saturation boundary in kg_water/kg_moistair [kg/kg]",\
 0.0, 0.0,1.0,0.1,0,2561)
DeclareVariable("sink_OD[14].medium.x_sat", "Steam water mass content of saturation boundary in kg_water/kg_dryair [kg/kg]",\
 0.0, 0.0,1.0,0.1,0,2561)
DeclareVariable("sink_OD[14].medium.p_steam_sat", "Partial saturation pressure of steam [Pa|bar]",\
 100000.0, 0.0,100000000.0,100000.0,0,2561)
DeclareAlias2("sink_OD[14].ports[1].m_flow", "Mass flow rate from the connection point into the component [kg/s]",\
 "fan_OD[14].ports[1].m_flow", -1, 5, 5414, 132)
DeclareAlias2("sink_OD[14].ports[1].p", "Thermodynamic pressure in the connection point [Pa|bar]",\
 "sink_OD[14].p", 1, 7, 414, 4)
DeclareVariable("sink_OD[14].ports[1].h_outflow", "Specific thermodynamic enthalpy close to the connection point if m_flow < 0 [J/kg]",\
 0.0, -10000000000.0,10000000000.0,1000000.0,0,521)
DeclareAlias2("sink_OD[14].ports[1].Xi_outflow[1]", "Independent mixture mass fractions m_i/m close to the connection point if m_flow < 0 [kg/kg]",\
 "sink_OD[14].X[1]", 1, 7, 416, 4)
DeclareVariable("sink_OD[14].flowDirection", "Allowed flow direction [:#(type=Modelica.Fluid.Types.PortFlowDirection)]",\
 3, 1.0,3.0,0.0,0,2565)
DeclareVariable("sink_OD[14].use_p_in", "Get the pressure from the input connector [:#(type=Boolean)]",\
 false, 0.0,0.0,0.0,0,1539)
DeclareVariable("sink_OD[14].use_T_in", "Get the temperature from the input connector [:#(type=Boolean)]",\
 false, 0.0,0.0,0.0,0,1539)
DeclareVariable("sink_OD[14].use_X_in", "Get the composition from the input connector [:#(type=Boolean)]",\
 false, 0.0,0.0,0.0,0,1539)
DeclareVariable("sink_OD[14].use_C_in", "Get the trace substances from the input connector [:#(type=Boolean)]",\
 false, 0.0,0.0,0.0,0,1539)
DeclareParameter("sink_OD[14].p", "Fixed value of pressure [Pa|bar]", 414, 101325,\
 0.0,100000000.0,100000.0,0,560)
DeclareParameter("sink_OD[14].T", "Fixed value of temperature [K|degC]", 415, \
293.15, 190.0,647.0,300.0,0,560)
DeclareParameter("sink_OD[14].X[1]", "Fixed value of composition [kg/kg]", 416, \
0.01, 0.0,1.0,0.1,0,560)
DeclareParameter("sink_OD[14].X[2]", "Fixed value of composition [kg/kg]", 417, \
0.99, 0.0,1.0,0.1,0,560)
DeclareAlias2("sink_OD[14].p_in_internal", "Needed to connect to conditional connector [Pa]",\
 "sink_OD[14].p", 1, 7, 414, 1024)
DeclareAlias2("sink_OD[14].T_in_internal", "Needed to connect to conditional connector [K]",\
 "sink_OD[14].T", 1, 7, 415, 1024)
DeclareVariable("sink_OD[14].X_in_internal[1]", "Needed to connect to conditional connector [1]",\
 0.0, 0.0,0.0,0.0,0,2561)
DeclareVariable("sink_OD[14].X_in_internal[2]", "Needed to connect to conditional connector [1]",\
 0.0, 0.0,0.0,0.0,0,2561)
DeclareVariable("sink_OD[15].nPorts", "Number of ports [:#(type=Integer)]", 1, \
0.0,0.0,0.0,0,517)
DeclareAlias2("sink_OD[15].medium.p", "Absolute pressure of medium [Pa|bar]", \
"sink_OD[15].p", 1, 7, 418, 0)
DeclareVariable("sink_OD[15].medium.Xi[1]", "Structurally independent mass fractions [1]",\
 0.01, 0.0,1.0,0.0,0,513)
DeclareAlias2("sink_OD[15].medium.h", "Specific enthalpy of medium [J/kg]", \
"sink_OD[15].ports[1].h_outflow", 1, 5, 5859, 0)
DeclareVariable("sink_OD[15].medium.d", "Density of medium [kg/m3|g/cm3]", 1, \
0.0,100000.0,1.0,0,513)
DeclareAlias2("sink_OD[15].medium.T", "Temperature of medium [K|degC]", \
"sink_OD[15].T", 1, 7, 419, 0)
DeclareAlias2("sink_OD[15].medium.X[1]", "Mass fractions (= (component mass)/total mass  m_i/m) [kg/kg]",\
 "sink_OD[15].X[1]", 1, 7, 420, 0)
DeclareVariable("sink_OD[15].medium.X[2]", "Mass fractions (= (component mass)/total mass  m_i/m) [kg/kg]",\
 0.99, 0.0,1.0,0.1,0,513)
DeclareVariable("sink_OD[15].medium.u", "Specific internal energy of medium [J/kg]",\
 0.0, -100000000.0,100000000.0,1000000.0,0,513)
DeclareVariable("sink_OD[15].medium.R_s", "Gas constant (of mixture if applicable) [J/(kg.K)]",\
 1000.0, 0.0,10000000.0,1000.0,0,513)
DeclareVariable("sink_OD[15].medium.MM", "Molar mass (of mixture or single fluid) [kg/mol]",\
 0.032, 0.001,0.25,0.032,0,513)
DeclareAlias2("sink_OD[15].medium.state.p", "Absolute pressure of medium [Pa|bar]",\
 "sink_OD[15].p", 1, 7, 418, 0)
DeclareAlias2("sink_OD[15].medium.state.T", "Temperature of medium [K|degC]", \
"sink_OD[15].T", 1, 7, 419, 0)
DeclareAlias2("sink_OD[15].medium.state.X[1]", "Mass fractions (= (component mass)/total mass  m_i/m) [kg/kg]",\
 "sink_OD[15].X[1]", 1, 7, 420, 0)
DeclareAlias2("sink_OD[15].medium.state.X[2]", "Mass fractions (= (component mass)/total mass  m_i/m) [kg/kg]",\
 "sink_OD[15].medium.X[2]", 1, 5, 5843, 0)
DeclareVariable("sink_OD[15].medium.preferredMediumStates", "= true if StateSelect.prefer shall be used for the independent property variables of the medium [:#(type=Boolean)]",\
 false, 0.0,0.0,0.0,0,515)
DeclareVariable("sink_OD[15].medium.standardOrderComponents", "If true, and reducedX = true, the last element of X will be computed from the other ones [:#(type=Boolean)]",\
 true, 0.0,0.0,0.0,0,515)
DeclareVariable("sink_OD[15].medium.T_degC", "Temperature of medium in [degC] [degC;]",\
 0.0, 0.0,0.0,0.0,0,513)
DeclareVariable("sink_OD[15].medium.p_bar", "Absolute pressure of medium in [bar] [bar]",\
 0.0, 0.0,0.0,0.0,0,513)
DeclareVariable("sink_OD[15].medium.x_water", "Mass of total water/mass of dry air [kg/kg]",\
 0.0, 0.0,1.0,0.1,0,513)
DeclareVariable("sink_OD[15].medium.phi", "Relative humidity", 0.0, 0.0,0.0,0.0,\
0,513)
DeclareVariable("sink_OD[15].medium.X_liquid", "Mass fraction of liquid or solid water [kg/kg]",\
 0.0, 0.0,1.0,0.1,0,2561)
DeclareVariable("sink_OD[15].medium.X_steam", "Mass fraction of steam water [kg/kg]",\
 0.0, 0.0,1.0,0.1,0,2561)
DeclareVariable("sink_OD[15].medium.X_air", "Mass fraction of air [kg/kg]", 0.0,\
 0.0,1.0,0.1,0,2561)
DeclareVariable("sink_OD[15].medium.X_sat", "Steam water mass fraction of saturation boundary in kg_water/kg_moistair [kg/kg]",\
 0.0, 0.0,1.0,0.1,0,2561)
DeclareVariable("sink_OD[15].medium.x_sat", "Steam water mass content of saturation boundary in kg_water/kg_dryair [kg/kg]",\
 0.0, 0.0,1.0,0.1,0,2561)
DeclareVariable("sink_OD[15].medium.p_steam_sat", "Partial saturation pressure of steam [Pa|bar]",\
 100000.0, 0.0,100000000.0,100000.0,0,2561)
DeclareAlias2("sink_OD[15].ports[1].m_flow", "Mass flow rate from the connection point into the component [kg/s]",\
 "fan_OD[15].ports[1].m_flow", -1, 5, 5449, 132)
DeclareAlias2("sink_OD[15].ports[1].p", "Thermodynamic pressure in the connection point [Pa|bar]",\
 "sink_OD[15].p", 1, 7, 418, 4)
DeclareVariable("sink_OD[15].ports[1].h_outflow", "Specific thermodynamic enthalpy close to the connection point if m_flow < 0 [J/kg]",\
 0.0, -10000000000.0,10000000000.0,1000000.0,0,521)
DeclareAlias2("sink_OD[15].ports[1].Xi_outflow[1]", "Independent mixture mass fractions m_i/m close to the connection point if m_flow < 0 [kg/kg]",\
 "sink_OD[15].X[1]", 1, 7, 420, 4)
DeclareVariable("sink_OD[15].flowDirection", "Allowed flow direction [:#(type=Modelica.Fluid.Types.PortFlowDirection)]",\
 3, 1.0,3.0,0.0,0,2565)
DeclareVariable("sink_OD[15].use_p_in", "Get the pressure from the input connector [:#(type=Boolean)]",\
 false, 0.0,0.0,0.0,0,1539)
DeclareVariable("sink_OD[15].use_T_in", "Get the temperature from the input connector [:#(type=Boolean)]",\
 false, 0.0,0.0,0.0,0,1539)
DeclareVariable("sink_OD[15].use_X_in", "Get the composition from the input connector [:#(type=Boolean)]",\
 false, 0.0,0.0,0.0,0,1539)
DeclareVariable("sink_OD[15].use_C_in", "Get the trace substances from the input connector [:#(type=Boolean)]",\
 false, 0.0,0.0,0.0,0,1539)
DeclareParameter("sink_OD[15].p", "Fixed value of pressure [Pa|bar]", 418, 101325,\
 0.0,100000000.0,100000.0,0,560)
DeclareParameter("sink_OD[15].T", "Fixed value of temperature [K|degC]", 419, \
293.15, 190.0,647.0,300.0,0,560)
DeclareParameter("sink_OD[15].X[1]", "Fixed value of composition [kg/kg]", 420, \
0.01, 0.0,1.0,0.1,0,560)
DeclareParameter("sink_OD[15].X[2]", "Fixed value of composition [kg/kg]", 421, \
0.99, 0.0,1.0,0.1,0,560)
DeclareAlias2("sink_OD[15].p_in_internal", "Needed to connect to conditional connector [Pa]",\
 "sink_OD[15].p", 1, 7, 418, 1024)
DeclareAlias2("sink_OD[15].T_in_internal", "Needed to connect to conditional connector [K]",\
 "sink_OD[15].T", 1, 7, 419, 1024)
DeclareVariable("sink_OD[15].X_in_internal[1]", "Needed to connect to conditional connector [1]",\
 0.0, 0.0,0.0,0.0,0,2561)
DeclareVariable("sink_OD[15].X_in_internal[2]", "Needed to connect to conditional connector [1]",\
 0.0, 0.0,0.0,0.0,0,2561)
DeclareVariable("charge", "[kg]", 0.0, 0.0,1E+100,0.0,0,512)
DeclareVariable("BC.nout", "Number of outputs [:#(type=Integer)]", 6, 1.0,1E+100,\
0.0,0,517)
DeclareAlias2("BC.y[1]", "Connector of Real output signals", "compressor.speed", 1,\
 5, 3175, 0)
DeclareAlias2("BC.y[2]", "Connector of Real output signals", "eXV.opening", 1, 5,\
 3228, 0)
DeclareVariable("BC.y[3]", "Connector of Real output signals", 0.0, 0.0,0.0,0.0,\
0,512)
DeclareVariable("BC.y[4]", "Connector of Real output signals", 0.0, 0.0,0.0,0.0,\
0,512)
DeclareVariable("BC.y[5]", "Connector of Real output signals", 288.15, 190.0,\
647.0,300.0,0,512)
DeclareAlias2("BC.y[6]", "Connector of Real output signals", "compressor.T_amb", 1,\
 5, 3176, 0)
DeclareVariable("BC.tableOnFile", "= true, if table is defined on file or in function usertab [:#(type=Boolean)]",\
 true, 0.0,0.0,0.0,0,515)
DeclareParameter("BC.verboseRead", "= true, if info message that file is loading is to be printed [:#(type=Boolean)]",\
 422, true, 0.0,0.0,0.0,0,562)
DeclareParameter("BC.columns[1]", "Columns of table to be interpolated [:#(type=Integer)]",\
 423, 2, 0.0,0.0,0.0,0,564)
DeclareParameter("BC.columns[2]", "Columns of table to be interpolated [:#(type=Integer)]",\
 424, 3, 0.0,0.0,0.0,0,564)
DeclareParameter("BC.columns[3]", "Columns of table to be interpolated [:#(type=Integer)]",\
 425, 4, 0.0,0.0,0.0,0,564)
DeclareParameter("BC.columns[4]", "Columns of table to be interpolated [:#(type=Integer)]",\
 426, 5, 0.0,0.0,0.0,0,564)
DeclareParameter("BC.columns[5]", "Columns of table to be interpolated [:#(type=Integer)]",\
 427, 6, 0.0,0.0,0.0,0,564)
DeclareParameter("BC.columns[6]", "Columns of table to be interpolated [:#(type=Integer)]",\
 428, 7, 0.0,0.0,0.0,0,564)
DeclareVariable("BC.smoothness", "Smoothness of table interpolation [:#(type=Modelica.Blocks.Types.Smoothness)]",\
 6, 1.0,6.0,0.0,0,517)
DeclareVariable("BC.extrapolation", "Extrapolation of data outside the definition range [:#(type=Modelica.Blocks.Types.Extrapolation)]",\
 2, 1.0,4.0,0.0,0,517)
DeclareVariable("BC.timeScale", "Time scale of first table column [s]", 1, 1E-15,\
1E+100,0.0,0,513)
DeclareParameter("BC.offset[1]", "Offsets of output signals", 429, 0, 0.0,0.0,\
0.0,0,560)
DeclareParameter("BC.startTime", "Output = offset for time < startTime [s]", 430,\
 0, 0.0,0.0,0.0,0,560)
DeclareVariable("BC.shiftTime", "Shift time of first table column [s]", 0.0, \
0.0,0.0,0.0,0,513)
DeclareParameter("BC.timeEvents", "Time event handling of table interpolation [:#(type=Modelica.Blocks.Types.TimeEvents)]",\
 431, 1, 1.0,3.0,0.0,0,564)
DeclareVariable("BC.verboseExtrapolation", "= true, if warning messages are to be printed if time is outside the table definition range [:#(type=Boolean)]",\
 false, 0.0,0.0,0.0,0,515)
DeclareVariable("BC.t_min", "Minimum abscissa value defined in table [s]", 0.0, \
0.0,0.0,0.0,0,513)
DeclareVariable("BC.t_max", "Maximum abscissa value defined in table [s]", 0.0, \
0.0,0.0,0.0,0,513)
DeclareVariable("BC.t_minScaled", "Minimum (scaled) abscissa value defined in table [1]",\
 0.0, 0.0,0.0,0.0,0,513)
DeclareVariable("BC.t_maxScaled", "Maximum (scaled) abscissa value defined in table [1]",\
 0.0, 0.0,0.0,0.0,0,513)
DeclareVariable("BC.p_offset[1]", "Offsets of output signals", 0.0, 0.0,0.0,0.0,\
0,2561)
DeclareVariable("BC.p_offset[2]", "Offsets of output signals", 0.0, 0.0,0.0,0.0,\
0,2561)
DeclareVariable("BC.p_offset[3]", "Offsets of output signals", 0.0, 0.0,0.0,0.0,\
0,2561)
DeclareVariable("BC.p_offset[4]", "Offsets of output signals", 0.0, 0.0,0.0,0.0,\
0,2561)
DeclareVariable("BC.p_offset[5]", "Offsets of output signals", 0.0, 0.0,0.0,0.0,\
0,2561)
DeclareVariable("BC.p_offset[6]", "Offsets of output signals", 0.0, 0.0,0.0,0.0,\
0,2561)
DeclareVariable("BC.tableID.id", "[:#(type=Integer)]", 0, 0.0,0.0,0.0,0,2565)
DeclareVariable("BC.nextTimeEvent", "Next time event instant [s]", 0, 0.0,0.0,\
0.0,0,2704)
DeclareVariable("BC.nextTimeEventScaled", "Next scaled time event instant [1]", 0,\
 0.0,0.0,0.0,0,2704)
DeclareVariable("BC.timeScaled", "Scaled time [1]", 0.0, 0.0,0.0,0.0,0,2560)
DeclareVariable("Mea.nout", "Number of outputs [:#(type=Integer)]", 8, 1.0,\
1E+100,0.0,0,517)
DeclareAlias2("Mea.y[1]", "Connector of Real output signals", "y_mea[1]", 1, 3, 6,\
 0)
DeclareAlias2("Mea.y[2]", "Connector of Real output signals", "y_mea[2]", 1, 3, 7,\
 0)
DeclareAlias2("Mea.y[3]", "Connector of Real output signals", "y_mea[3]", 1, 3, 8,\
 0)
DeclareAlias2("Mea.y[4]", "Connector of Real output signals", "y_mea[4]", 1, 3, 9,\
 0)
DeclareVariable("Mea.y[5]", "Connector of Real output signals", 0.0, 0.0,0.0,0.0,\
0,512)
DeclareAlias2("Mea.y[6]", "Connector of Real output signals", "y_mea[5]", 1, 3, 10,\
 0)
DeclareAlias2("Mea.y[7]", "Connector of Real output signals", "y_mea[6]", 1, 3, 11,\
 0)
DeclareVariable("Mea.y[8]", "Connector of Real output signals", 0.0, 0.0,0.0,0.0,\
0,512)
DeclareVariable("Mea.tableOnFile", "= true, if table is defined on file or in function usertab [:#(type=Boolean)]",\
 true, 0.0,0.0,0.0,0,515)
DeclareParameter("Mea.verboseRead", "= true, if info message that file is loading is to be printed [:#(type=Boolean)]",\
 432, true, 0.0,0.0,0.0,0,562)
DeclareParameter("Mea.columns[1]", "Columns of table to be interpolated [:#(type=Integer)]",\
 433, 2, 0.0,0.0,0.0,0,564)
DeclareParameter("Mea.columns[2]", "Columns of table to be interpolated [:#(type=Integer)]",\
 434, 3, 0.0,0.0,0.0,0,564)
DeclareParameter("Mea.columns[3]", "Columns of table to be interpolated [:#(type=Integer)]",\
 435, 4, 0.0,0.0,0.0,0,564)
DeclareParameter("Mea.columns[4]", "Columns of table to be interpolated [:#(type=Integer)]",\
 436, 5, 0.0,0.0,0.0,0,564)
DeclareParameter("Mea.columns[5]", "Columns of table to be interpolated [:#(type=Integer)]",\
 437, 6, 0.0,0.0,0.0,0,564)
DeclareParameter("Mea.columns[6]", "Columns of table to be interpolated [:#(type=Integer)]",\
 438, 7, 0.0,0.0,0.0,0,564)
DeclareParameter("Mea.columns[7]", "Columns of table to be interpolated [:#(type=Integer)]",\
 439, 8, 0.0,0.0,0.0,0,564)
DeclareParameter("Mea.columns[8]", "Columns of table to be interpolated [:#(type=Integer)]",\
 440, 9, 0.0,0.0,0.0,0,564)
DeclareVariable("Mea.smoothness", "Smoothness of table interpolation [:#(type=Modelica.Blocks.Types.Smoothness)]",\
 6, 1.0,6.0,0.0,0,517)
DeclareVariable("Mea.extrapolation", "Extrapolation of data outside the definition range [:#(type=Modelica.Blocks.Types.Extrapolation)]",\
 2, 1.0,4.0,0.0,0,517)
DeclareVariable("Mea.timeScale", "Time scale of first table column [s]", 1, \
1E-15,1E+100,0.0,0,513)
DeclareParameter("Mea.offset[1]", "Offsets of output signals", 441, 0, 0.0,0.0,\
0.0,0,560)
DeclareParameter("Mea.startTime", "Output = offset for time < startTime [s]", 442,\
 0, 0.0,0.0,0.0,0,560)
DeclareVariable("Mea.shiftTime", "Shift time of first table column [s]", 0.0, \
0.0,0.0,0.0,0,513)
DeclareParameter("Mea.timeEvents", "Time event handling of table interpolation [:#(type=Modelica.Blocks.Types.TimeEvents)]",\
 443, 1, 1.0,3.0,0.0,0,564)
DeclareVariable("Mea.verboseExtrapolation", "= true, if warning messages are to be printed if time is outside the table definition range [:#(type=Boolean)]",\
 false, 0.0,0.0,0.0,0,515)
DeclareVariable("Mea.t_min", "Minimum abscissa value defined in table [s]", 0.0,\
 0.0,0.0,0.0,0,513)
DeclareVariable("Mea.t_max", "Maximum abscissa value defined in table [s]", 0.0,\
 0.0,0.0,0.0,0,513)
DeclareVariable("Mea.t_minScaled", "Minimum (scaled) abscissa value defined in table [1]",\
 0.0, 0.0,0.0,0.0,0,513)
DeclareVariable("Mea.t_maxScaled", "Maximum (scaled) abscissa value defined in table [1]",\
 0.0, 0.0,0.0,0.0,0,513)
DeclareVariable("Mea.p_offset[1]", "Offsets of output signals", 0.0, 0.0,0.0,0.0,\
0,2561)
DeclareVariable("Mea.p_offset[2]", "Offsets of output signals", 0.0, 0.0,0.0,0.0,\
0,2561)
DeclareVariable("Mea.p_offset[3]", "Offsets of output signals", 0.0, 0.0,0.0,0.0,\
0,2561)
DeclareVariable("Mea.p_offset[4]", "Offsets of output signals", 0.0, 0.0,0.0,0.0,\
0,2561)
DeclareVariable("Mea.p_offset[5]", "Offsets of output signals", 0.0, 0.0,0.0,0.0,\
0,2561)
DeclareVariable("Mea.p_offset[6]", "Offsets of output signals", 0.0, 0.0,0.0,0.0,\
0,2561)
DeclareVariable("Mea.p_offset[7]", "Offsets of output signals", 0.0, 0.0,0.0,0.0,\
0,2561)
DeclareVariable("Mea.p_offset[8]", "Offsets of output signals", 0.0, 0.0,0.0,0.0,\
0,2561)
DeclareVariable("Mea.tableID.id", "[:#(type=Integer)]", 0, 0.0,0.0,0.0,0,2565)
DeclareVariable("Mea.nextTimeEvent", "Next time event instant [s]", 0, 0.0,0.0,\
0.0,0,2704)
DeclareVariable("Mea.nextTimeEventScaled", "Next scaled time event instant [1]",\
 0, 0.0,0.0,0.0,0,2704)
DeclareVariable("Mea.timeScaled", "Scaled time [1]", 0.0, 0.0,0.0,0.0,0,2560)
DeclareVariable("subcooling.port.m_flow", "Mass flow rate from the connection point into the component [kg/s]",\
 0, 0.0,100000.0,0.0,0,777)
DeclareAlias2("subcooling.port.p", "Thermodynamic pressure in the connection point [Pa|bar]",\
 "indoorCoil.refFlow.p[15]", 1, 1, 73, 4)
DeclareVariable("subcooling.port.h_outflow", "Specific thermodynamic enthalpy close to the connection point if m_flow < 0 [J/kg]",\
 400000.0, 100000.0,500000.0,500000.0,0,521)
DeclareAlias2("subcooling.T", "superheating degree [K]", "y[5]", 1, 3, 4, 0)
DeclareVariable("subcooling.T_sat", "[K|degC]", 298, 1.0,10000.0,350.0,0,2560)
DeclareVariable("outdoorCoil.refFlow.states[1].phase", "Phase of the fluid: 1 for 1-phase, 2 for two-phase, 0 for not known, e.g., interactive use [:#(type=Integer)]",\
 0, 0.0,2.0,0.0,0,644)
DeclareVariable("outdoorCoil.refFlow.states[1].p", "Pressure [Pa|bar]", \
1000000.0, 100000.0,4800000.0,4000000.0,0,512)
DeclareVariable("outdoorCoil.refFlow.states[1].h", "Specific enthalpy [J/kg]", \
400000.0, 100000.0,500000.0,500000.0,0,512)
DeclareVariable("outdoorCoil.refFlow.states[1].d", "Density [kg/m3|g/cm3]", 500,\
 0.0,100000.0,500.0,0,512)
DeclareVariable("outdoorCoil.refFlow.states[1].T", "Temperature [K|degC]", 298, \
1.0,10000.0,350.0,0,512)
DeclareVariable("outdoorCoil.refFlow.states[2].phase", "Phase of the fluid: 1 for 1-phase, 2 for two-phase, 0 for not known, e.g., interactive use [:#(type=Integer)]",\
 0, 0.0,2.0,0.0,0,644)
DeclareVariable("outdoorCoil.refFlow.states[2].p", "Pressure [Pa|bar]", \
1000000.0, 100000.0,4800000.0,4000000.0,0,512)
DeclareVariable("outdoorCoil.refFlow.states[2].h", "Specific enthalpy [J/kg]", \
400000.0, 100000.0,500000.0,500000.0,0,512)
DeclareVariable("outdoorCoil.refFlow.states[2].d", "Density [kg/m3|g/cm3]", 500,\
 0.0,100000.0,500.0,0,512)
DeclareVariable("outdoorCoil.refFlow.states[2].T", "Temperature [K|degC]", 298, \
1.0,10000.0,350.0,0,512)
DeclareVariable("outdoorCoil.refFlow.states[3].phase", "Phase of the fluid: 1 for 1-phase, 2 for two-phase, 0 for not known, e.g., interactive use [:#(type=Integer)]",\
 0, 0.0,2.0,0.0,0,644)
DeclareVariable("outdoorCoil.refFlow.states[3].p", "Pressure [Pa|bar]", \
1000000.0, 100000.0,4800000.0,4000000.0,0,512)
DeclareVariable("outdoorCoil.refFlow.states[3].h", "Specific enthalpy [J/kg]", \
400000.0, 100000.0,500000.0,500000.0,0,512)
DeclareVariable("outdoorCoil.refFlow.states[3].d", "Density [kg/m3|g/cm3]", 500,\
 0.0,100000.0,500.0,0,512)
DeclareVariable("outdoorCoil.refFlow.states[3].T", "Temperature [K|degC]", 298, \
1.0,10000.0,350.0,0,512)
DeclareVariable("outdoorCoil.refFlow.states[4].phase", "Phase of the fluid: 1 for 1-phase, 2 for two-phase, 0 for not known, e.g., interactive use [:#(type=Integer)]",\
 0, 0.0,2.0,0.0,0,644)
DeclareVariable("outdoorCoil.refFlow.states[4].p", "Pressure [Pa|bar]", \
1000000.0, 100000.0,4800000.0,4000000.0,0,512)
DeclareVariable("outdoorCoil.refFlow.states[4].h", "Specific enthalpy [J/kg]", \
400000.0, 100000.0,500000.0,500000.0,0,512)
DeclareVariable("outdoorCoil.refFlow.states[4].d", "Density [kg/m3|g/cm3]", 500,\
 0.0,100000.0,500.0,0,512)
DeclareVariable("outdoorCoil.refFlow.states[4].T", "Temperature [K|degC]", 298, \
1.0,10000.0,350.0,0,512)
DeclareVariable("outdoorCoil.refFlow.states[5].phase", "Phase of the fluid: 1 for 1-phase, 2 for two-phase, 0 for not known, e.g., interactive use [:#(type=Integer)]",\
 0, 0.0,2.0,0.0,0,644)
DeclareVariable("outdoorCoil.refFlow.states[5].p", "Pressure [Pa|bar]", \
1000000.0, 100000.0,4800000.0,4000000.0,0,512)
DeclareVariable("outdoorCoil.refFlow.states[5].h", "Specific enthalpy [J/kg]", \
400000.0, 100000.0,500000.0,500000.0,0,512)
DeclareVariable("outdoorCoil.refFlow.states[5].d", "Density [kg/m3|g/cm3]", 500,\
 0.0,100000.0,500.0,0,512)
DeclareVariable("outdoorCoil.refFlow.states[5].T", "Temperature [K|degC]", 298, \
1.0,10000.0,350.0,0,512)
DeclareVariable("outdoorCoil.refFlow.states[6].phase", "Phase of the fluid: 1 for 1-phase, 2 for two-phase, 0 for not known, e.g., interactive use [:#(type=Integer)]",\
 0, 0.0,2.0,0.0,0,644)
DeclareVariable("outdoorCoil.refFlow.states[6].p", "Pressure [Pa|bar]", \
1000000.0, 100000.0,4800000.0,4000000.0,0,512)
DeclareVariable("outdoorCoil.refFlow.states[6].h", "Specific enthalpy [J/kg]", \
400000.0, 100000.0,500000.0,500000.0,0,512)
DeclareVariable("outdoorCoil.refFlow.states[6].d", "Density [kg/m3|g/cm3]", 500,\
 0.0,100000.0,500.0,0,512)
DeclareVariable("outdoorCoil.refFlow.states[6].T", "Temperature [K|degC]", 298, \
1.0,10000.0,350.0,0,512)
DeclareVariable("outdoorCoil.refFlow.states[7].phase", "Phase of the fluid: 1 for 1-phase, 2 for two-phase, 0 for not known, e.g., interactive use [:#(type=Integer)]",\
 0, 0.0,2.0,0.0,0,644)
DeclareVariable("outdoorCoil.refFlow.states[7].p", "Pressure [Pa|bar]", \
1000000.0, 100000.0,4800000.0,4000000.0,0,512)
DeclareVariable("outdoorCoil.refFlow.states[7].h", "Specific enthalpy [J/kg]", \
400000.0, 100000.0,500000.0,500000.0,0,512)
DeclareVariable("outdoorCoil.refFlow.states[7].d", "Density [kg/m3|g/cm3]", 500,\
 0.0,100000.0,500.0,0,512)
DeclareVariable("outdoorCoil.refFlow.states[7].T", "Temperature [K|degC]", 298, \
1.0,10000.0,350.0,0,512)
DeclareVariable("outdoorCoil.refFlow.states[8].phase", "Phase of the fluid: 1 for 1-phase, 2 for two-phase, 0 for not known, e.g., interactive use [:#(type=Integer)]",\
 0, 0.0,2.0,0.0,0,644)
DeclareVariable("outdoorCoil.refFlow.states[8].p", "Pressure [Pa|bar]", \
1000000.0, 100000.0,4800000.0,4000000.0,0,512)
DeclareVariable("outdoorCoil.refFlow.states[8].h", "Specific enthalpy [J/kg]", \
400000.0, 100000.0,500000.0,500000.0,0,512)
DeclareVariable("outdoorCoil.refFlow.states[8].d", "Density [kg/m3|g/cm3]", 500,\
 0.0,100000.0,500.0,0,512)
DeclareVariable("outdoorCoil.refFlow.states[8].T", "Temperature [K|degC]", 298, \
1.0,10000.0,350.0,0,512)
DeclareVariable("outdoorCoil.refFlow.states[9].phase", "Phase of the fluid: 1 for 1-phase, 2 for two-phase, 0 for not known, e.g., interactive use [:#(type=Integer)]",\
 0, 0.0,2.0,0.0,0,644)
DeclareVariable("outdoorCoil.refFlow.states[9].p", "Pressure [Pa|bar]", \
1000000.0, 100000.0,4800000.0,4000000.0,0,512)
DeclareVariable("outdoorCoil.refFlow.states[9].h", "Specific enthalpy [J/kg]", \
400000.0, 100000.0,500000.0,500000.0,0,512)
DeclareVariable("outdoorCoil.refFlow.states[9].d", "Density [kg/m3|g/cm3]", 500,\
 0.0,100000.0,500.0,0,512)
DeclareVariable("outdoorCoil.refFlow.states[9].T", "Temperature [K|degC]", 298, \
1.0,10000.0,350.0,0,512)
DeclareVariable("outdoorCoil.refFlow.states[10].phase", "Phase of the fluid: 1 for 1-phase, 2 for two-phase, 0 for not known, e.g., interactive use [:#(type=Integer)]",\
 0, 0.0,2.0,0.0,0,644)
DeclareVariable("outdoorCoil.refFlow.states[10].p", "Pressure [Pa|bar]", \
1000000.0, 100000.0,4800000.0,4000000.0,0,512)
DeclareVariable("outdoorCoil.refFlow.states[10].h", "Specific enthalpy [J/kg]", \
400000.0, 100000.0,500000.0,500000.0,0,512)
DeclareVariable("outdoorCoil.refFlow.states[10].d", "Density [kg/m3|g/cm3]", 500,\
 0.0,100000.0,500.0,0,512)
DeclareVariable("outdoorCoil.refFlow.states[10].T", "Temperature [K|degC]", 298,\
 1.0,10000.0,350.0,0,512)
DeclareVariable("outdoorCoil.refFlow.states[11].phase", "Phase of the fluid: 1 for 1-phase, 2 for two-phase, 0 for not known, e.g., interactive use [:#(type=Integer)]",\
 0, 0.0,2.0,0.0,0,644)
DeclareVariable("outdoorCoil.refFlow.states[11].p", "Pressure [Pa|bar]", \
1000000.0, 100000.0,4800000.0,4000000.0,0,512)
DeclareVariable("outdoorCoil.refFlow.states[11].h", "Specific enthalpy [J/kg]", \
400000.0, 100000.0,500000.0,500000.0,0,512)
DeclareVariable("outdoorCoil.refFlow.states[11].d", "Density [kg/m3|g/cm3]", 500,\
 0.0,100000.0,500.0,0,512)
DeclareVariable("outdoorCoil.refFlow.states[11].T", "Temperature [K|degC]", 298,\
 1.0,10000.0,350.0,0,512)
DeclareVariable("outdoorCoil.refFlow.states[12].phase", "Phase of the fluid: 1 for 1-phase, 2 for two-phase, 0 for not known, e.g., interactive use [:#(type=Integer)]",\
 0, 0.0,2.0,0.0,0,644)
DeclareVariable("outdoorCoil.refFlow.states[12].p", "Pressure [Pa|bar]", \
1000000.0, 100000.0,4800000.0,4000000.0,0,512)
DeclareVariable("outdoorCoil.refFlow.states[12].h", "Specific enthalpy [J/kg]", \
400000.0, 100000.0,500000.0,500000.0,0,512)
DeclareVariable("outdoorCoil.refFlow.states[12].d", "Density [kg/m3|g/cm3]", 500,\
 0.0,100000.0,500.0,0,512)
DeclareVariable("outdoorCoil.refFlow.states[12].T", "Temperature [K|degC]", 298,\
 1.0,10000.0,350.0,0,512)
DeclareVariable("outdoorCoil.refFlow.states[13].phase", "Phase of the fluid: 1 for 1-phase, 2 for two-phase, 0 for not known, e.g., interactive use [:#(type=Integer)]",\
 0, 0.0,2.0,0.0,0,644)
DeclareVariable("outdoorCoil.refFlow.states[13].p", "Pressure [Pa|bar]", \
1000000.0, 100000.0,4800000.0,4000000.0,0,512)
DeclareVariable("outdoorCoil.refFlow.states[13].h", "Specific enthalpy [J/kg]", \
400000.0, 100000.0,500000.0,500000.0,0,512)
DeclareVariable("outdoorCoil.refFlow.states[13].d", "Density [kg/m3|g/cm3]", 500,\
 0.0,100000.0,500.0,0,512)
DeclareVariable("outdoorCoil.refFlow.states[13].T", "Temperature [K|degC]", 298,\
 1.0,10000.0,350.0,0,512)
DeclareVariable("outdoorCoil.refFlow.states[14].phase", "Phase of the fluid: 1 for 1-phase, 2 for two-phase, 0 for not known, e.g., interactive use [:#(type=Integer)]",\
 0, 0.0,2.0,0.0,0,644)
DeclareVariable("outdoorCoil.refFlow.states[14].p", "Pressure [Pa|bar]", \
1000000.0, 100000.0,4800000.0,4000000.0,0,512)
DeclareVariable("outdoorCoil.refFlow.states[14].h", "Specific enthalpy [J/kg]", \
400000.0, 100000.0,500000.0,500000.0,0,512)
DeclareVariable("outdoorCoil.refFlow.states[14].d", "Density [kg/m3|g/cm3]", 500,\
 0.0,100000.0,500.0,0,512)
DeclareVariable("outdoorCoil.refFlow.states[14].T", "Temperature [K|degC]", 298,\
 1.0,10000.0,350.0,0,512)
DeclareVariable("outdoorCoil.refFlow.states[15].phase", "Phase of the fluid: 1 for 1-phase, 2 for two-phase, 0 for not known, e.g., interactive use [:#(type=Integer)]",\
 0, 0.0,2.0,0.0,0,644)
DeclareVariable("outdoorCoil.refFlow.states[15].p", "Pressure [Pa|bar]", \
1000000.0, 100000.0,4800000.0,4000000.0,0,512)
DeclareVariable("outdoorCoil.refFlow.states[15].h", "Specific enthalpy [J/kg]", \
400000.0, 100000.0,500000.0,500000.0,0,512)
DeclareVariable("outdoorCoil.refFlow.states[15].d", "Density [kg/m3|g/cm3]", 500,\
 0.0,100000.0,500.0,0,512)
DeclareVariable("outdoorCoil.refFlow.states[15].T", "Temperature [K|degC]", 298,\
 1.0,10000.0,350.0,0,512)
DeclareAlias2("outdoorCoil.refFlow.slipRatio.states[1].phase", "Phase of the fluid: 1 for 1-phase, 2 for two-phase, 0 for not known, e.g., interactive use [:#(type=Integer)]",\
 "outdoorCoil.refFlow.states[1].phase", 1, 5, 5920, 66)
DeclareAlias2("outdoorCoil.refFlow.slipRatio.states[1].p", "Pressure [Pa|bar]", \
"outdoorCoil.refFlow.states[1].p", 1, 5, 5921, 0)
DeclareAlias2("outdoorCoil.refFlow.slipRatio.states[1].h", "Specific enthalpy [J/kg]",\
 "outdoorCoil.refFlow.states[1].h", 1, 5, 5922, 0)
DeclareAlias2("outdoorCoil.refFlow.slipRatio.states[1].d", "Density [kg/m3|g/cm3]",\
 "outdoorCoil.refFlow.states[1].d", 1, 5, 5923, 0)
DeclareAlias2("outdoorCoil.refFlow.slipRatio.states[1].T", "Temperature [K|degC]",\
 "outdoorCoil.refFlow.states[1].T", 1, 5, 5924, 0)
DeclareAlias2("outdoorCoil.refFlow.slipRatio.states[2].phase", "Phase of the fluid: 1 for 1-phase, 2 for two-phase, 0 for not known, e.g., interactive use [:#(type=Integer)]",\
 "outdoorCoil.refFlow.states[2].phase", 1, 5, 5925, 66)
DeclareAlias2("outdoorCoil.refFlow.slipRatio.states[2].p", "Pressure [Pa|bar]", \
"outdoorCoil.refFlow.states[2].p", 1, 5, 5926, 0)
DeclareAlias2("outdoorCoil.refFlow.slipRatio.states[2].h", "Specific enthalpy [J/kg]",\
 "outdoorCoil.refFlow.states[2].h", 1, 5, 5927, 0)
DeclareAlias2("outdoorCoil.refFlow.slipRatio.states[2].d", "Density [kg/m3|g/cm3]",\
 "outdoorCoil.refFlow.states[2].d", 1, 5, 5928, 0)
DeclareAlias2("outdoorCoil.refFlow.slipRatio.states[2].T", "Temperature [K|degC]",\
 "outdoorCoil.refFlow.states[2].T", 1, 5, 5929, 0)
DeclareAlias2("outdoorCoil.refFlow.slipRatio.states[3].phase", "Phase of the fluid: 1 for 1-phase, 2 for two-phase, 0 for not known, e.g., interactive use [:#(type=Integer)]",\
 "outdoorCoil.refFlow.states[3].phase", 1, 5, 5930, 66)
DeclareAlias2("outdoorCoil.refFlow.slipRatio.states[3].p", "Pressure [Pa|bar]", \
"outdoorCoil.refFlow.states[3].p", 1, 5, 5931, 0)
DeclareAlias2("outdoorCoil.refFlow.slipRatio.states[3].h", "Specific enthalpy [J/kg]",\
 "outdoorCoil.refFlow.states[3].h", 1, 5, 5932, 0)
DeclareAlias2("outdoorCoil.refFlow.slipRatio.states[3].d", "Density [kg/m3|g/cm3]",\
 "outdoorCoil.refFlow.states[3].d", 1, 5, 5933, 0)
DeclareAlias2("outdoorCoil.refFlow.slipRatio.states[3].T", "Temperature [K|degC]",\
 "outdoorCoil.refFlow.states[3].T", 1, 5, 5934, 0)
DeclareAlias2("outdoorCoil.refFlow.slipRatio.states[4].phase", "Phase of the fluid: 1 for 1-phase, 2 for two-phase, 0 for not known, e.g., interactive use [:#(type=Integer)]",\
 "outdoorCoil.refFlow.states[4].phase", 1, 5, 5935, 66)
DeclareAlias2("outdoorCoil.refFlow.slipRatio.states[4].p", "Pressure [Pa|bar]", \
"outdoorCoil.refFlow.states[4].p", 1, 5, 5936, 0)
DeclareAlias2("outdoorCoil.refFlow.slipRatio.states[4].h", "Specific enthalpy [J/kg]",\
 "outdoorCoil.refFlow.states[4].h", 1, 5, 5937, 0)
DeclareAlias2("outdoorCoil.refFlow.slipRatio.states[4].d", "Density [kg/m3|g/cm3]",\
 "outdoorCoil.refFlow.states[4].d", 1, 5, 5938, 0)
DeclareAlias2("outdoorCoil.refFlow.slipRatio.states[4].T", "Temperature [K|degC]",\
 "outdoorCoil.refFlow.states[4].T", 1, 5, 5939, 0)
DeclareAlias2("outdoorCoil.refFlow.slipRatio.states[5].phase", "Phase of the fluid: 1 for 1-phase, 2 for two-phase, 0 for not known, e.g., interactive use [:#(type=Integer)]",\
 "outdoorCoil.refFlow.states[5].phase", 1, 5, 5940, 66)
DeclareAlias2("outdoorCoil.refFlow.slipRatio.states[5].p", "Pressure [Pa|bar]", \
"outdoorCoil.refFlow.states[5].p", 1, 5, 5941, 0)
DeclareAlias2("outdoorCoil.refFlow.slipRatio.states[5].h", "Specific enthalpy [J/kg]",\
 "outdoorCoil.refFlow.states[5].h", 1, 5, 5942, 0)
DeclareAlias2("outdoorCoil.refFlow.slipRatio.states[5].d", "Density [kg/m3|g/cm3]",\
 "outdoorCoil.refFlow.states[5].d", 1, 5, 5943, 0)
DeclareAlias2("outdoorCoil.refFlow.slipRatio.states[5].T", "Temperature [K|degC]",\
 "outdoorCoil.refFlow.states[5].T", 1, 5, 5944, 0)
DeclareAlias2("outdoorCoil.refFlow.slipRatio.states[6].phase", "Phase of the fluid: 1 for 1-phase, 2 for two-phase, 0 for not known, e.g., interactive use [:#(type=Integer)]",\
 "outdoorCoil.refFlow.states[6].phase", 1, 5, 5945, 66)
DeclareAlias2("outdoorCoil.refFlow.slipRatio.states[6].p", "Pressure [Pa|bar]", \
"outdoorCoil.refFlow.states[6].p", 1, 5, 5946, 0)
DeclareAlias2("outdoorCoil.refFlow.slipRatio.states[6].h", "Specific enthalpy [J/kg]",\
 "outdoorCoil.refFlow.states[6].h", 1, 5, 5947, 0)
DeclareAlias2("outdoorCoil.refFlow.slipRatio.states[6].d", "Density [kg/m3|g/cm3]",\
 "outdoorCoil.refFlow.states[6].d", 1, 5, 5948, 0)
DeclareAlias2("outdoorCoil.refFlow.slipRatio.states[6].T", "Temperature [K|degC]",\
 "outdoorCoil.refFlow.states[6].T", 1, 5, 5949, 0)
DeclareAlias2("outdoorCoil.refFlow.slipRatio.states[7].phase", "Phase of the fluid: 1 for 1-phase, 2 for two-phase, 0 for not known, e.g., interactive use [:#(type=Integer)]",\
 "outdoorCoil.refFlow.states[7].phase", 1, 5, 5950, 66)
DeclareAlias2("outdoorCoil.refFlow.slipRatio.states[7].p", "Pressure [Pa|bar]", \
"outdoorCoil.refFlow.states[7].p", 1, 5, 5951, 0)
DeclareAlias2("outdoorCoil.refFlow.slipRatio.states[7].h", "Specific enthalpy [J/kg]",\
 "outdoorCoil.refFlow.states[7].h", 1, 5, 5952, 0)
DeclareAlias2("outdoorCoil.refFlow.slipRatio.states[7].d", "Density [kg/m3|g/cm3]",\
 "outdoorCoil.refFlow.states[7].d", 1, 5, 5953, 0)
DeclareAlias2("outdoorCoil.refFlow.slipRatio.states[7].T", "Temperature [K|degC]",\
 "outdoorCoil.refFlow.states[7].T", 1, 5, 5954, 0)
DeclareAlias2("outdoorCoil.refFlow.slipRatio.states[8].phase", "Phase of the fluid: 1 for 1-phase, 2 for two-phase, 0 for not known, e.g., interactive use [:#(type=Integer)]",\
 "outdoorCoil.refFlow.states[8].phase", 1, 5, 5955, 66)
DeclareAlias2("outdoorCoil.refFlow.slipRatio.states[8].p", "Pressure [Pa|bar]", \
"outdoorCoil.refFlow.states[8].p", 1, 5, 5956, 0)
DeclareAlias2("outdoorCoil.refFlow.slipRatio.states[8].h", "Specific enthalpy [J/kg]",\
 "outdoorCoil.refFlow.states[8].h", 1, 5, 5957, 0)
DeclareAlias2("outdoorCoil.refFlow.slipRatio.states[8].d", "Density [kg/m3|g/cm3]",\
 "outdoorCoil.refFlow.states[8].d", 1, 5, 5958, 0)
DeclareAlias2("outdoorCoil.refFlow.slipRatio.states[8].T", "Temperature [K|degC]",\
 "outdoorCoil.refFlow.states[8].T", 1, 5, 5959, 0)
DeclareAlias2("outdoorCoil.refFlow.slipRatio.states[9].phase", "Phase of the fluid: 1 for 1-phase, 2 for two-phase, 0 for not known, e.g., interactive use [:#(type=Integer)]",\
 "outdoorCoil.refFlow.states[9].phase", 1, 5, 5960, 66)
DeclareAlias2("outdoorCoil.refFlow.slipRatio.states[9].p", "Pressure [Pa|bar]", \
"outdoorCoil.refFlow.states[9].p", 1, 5, 5961, 0)
DeclareAlias2("outdoorCoil.refFlow.slipRatio.states[9].h", "Specific enthalpy [J/kg]",\
 "outdoorCoil.refFlow.states[9].h", 1, 5, 5962, 0)
DeclareAlias2("outdoorCoil.refFlow.slipRatio.states[9].d", "Density [kg/m3|g/cm3]",\
 "outdoorCoil.refFlow.states[9].d", 1, 5, 5963, 0)
DeclareAlias2("outdoorCoil.refFlow.slipRatio.states[9].T", "Temperature [K|degC]",\
 "outdoorCoil.refFlow.states[9].T", 1, 5, 5964, 0)
DeclareAlias2("outdoorCoil.refFlow.slipRatio.states[10].phase", "Phase of the fluid: 1 for 1-phase, 2 for two-phase, 0 for not known, e.g., interactive use [:#(type=Integer)]",\
 "outdoorCoil.refFlow.states[10].phase", 1, 5, 5965, 66)
DeclareAlias2("outdoorCoil.refFlow.slipRatio.states[10].p", "Pressure [Pa|bar]",\
 "outdoorCoil.refFlow.states[10].p", 1, 5, 5966, 0)
DeclareAlias2("outdoorCoil.refFlow.slipRatio.states[10].h", "Specific enthalpy [J/kg]",\
 "outdoorCoil.refFlow.states[10].h", 1, 5, 5967, 0)
DeclareAlias2("outdoorCoil.refFlow.slipRatio.states[10].d", "Density [kg/m3|g/cm3]",\
 "outdoorCoil.refFlow.states[10].d", 1, 5, 5968, 0)
DeclareAlias2("outdoorCoil.refFlow.slipRatio.states[10].T", "Temperature [K|degC]",\
 "outdoorCoil.refFlow.states[10].T", 1, 5, 5969, 0)
DeclareAlias2("outdoorCoil.refFlow.slipRatio.states[11].phase", "Phase of the fluid: 1 for 1-phase, 2 for two-phase, 0 for not known, e.g., interactive use [:#(type=Integer)]",\
 "outdoorCoil.refFlow.states[11].phase", 1, 5, 5970, 66)
DeclareAlias2("outdoorCoil.refFlow.slipRatio.states[11].p", "Pressure [Pa|bar]",\
 "outdoorCoil.refFlow.states[11].p", 1, 5, 5971, 0)
DeclareAlias2("outdoorCoil.refFlow.slipRatio.states[11].h", "Specific enthalpy [J/kg]",\
 "outdoorCoil.refFlow.states[11].h", 1, 5, 5972, 0)
DeclareAlias2("outdoorCoil.refFlow.slipRatio.states[11].d", "Density [kg/m3|g/cm3]",\
 "outdoorCoil.refFlow.states[11].d", 1, 5, 5973, 0)
DeclareAlias2("outdoorCoil.refFlow.slipRatio.states[11].T", "Temperature [K|degC]",\
 "outdoorCoil.refFlow.states[11].T", 1, 5, 5974, 0)
DeclareAlias2("outdoorCoil.refFlow.slipRatio.states[12].phase", "Phase of the fluid: 1 for 1-phase, 2 for two-phase, 0 for not known, e.g., interactive use [:#(type=Integer)]",\
 "outdoorCoil.refFlow.states[12].phase", 1, 5, 5975, 66)
DeclareAlias2("outdoorCoil.refFlow.slipRatio.states[12].p", "Pressure [Pa|bar]",\
 "outdoorCoil.refFlow.states[12].p", 1, 5, 5976, 0)
DeclareAlias2("outdoorCoil.refFlow.slipRatio.states[12].h", "Specific enthalpy [J/kg]",\
 "outdoorCoil.refFlow.states[12].h", 1, 5, 5977, 0)
DeclareAlias2("outdoorCoil.refFlow.slipRatio.states[12].d", "Density [kg/m3|g/cm3]",\
 "outdoorCoil.refFlow.states[12].d", 1, 5, 5978, 0)
DeclareAlias2("outdoorCoil.refFlow.slipRatio.states[12].T", "Temperature [K|degC]",\
 "outdoorCoil.refFlow.states[12].T", 1, 5, 5979, 0)
DeclareAlias2("outdoorCoil.refFlow.slipRatio.states[13].phase", "Phase of the fluid: 1 for 1-phase, 2 for two-phase, 0 for not known, e.g., interactive use [:#(type=Integer)]",\
 "outdoorCoil.refFlow.states[13].phase", 1, 5, 5980, 66)
DeclareAlias2("outdoorCoil.refFlow.slipRatio.states[13].p", "Pressure [Pa|bar]",\
 "outdoorCoil.refFlow.states[13].p", 1, 5, 5981, 0)
DeclareAlias2("outdoorCoil.refFlow.slipRatio.states[13].h", "Specific enthalpy [J/kg]",\
 "outdoorCoil.refFlow.states[13].h", 1, 5, 5982, 0)
DeclareAlias2("outdoorCoil.refFlow.slipRatio.states[13].d", "Density [kg/m3|g/cm3]",\
 "outdoorCoil.refFlow.states[13].d", 1, 5, 5983, 0)
DeclareAlias2("outdoorCoil.refFlow.slipRatio.states[13].T", "Temperature [K|degC]",\
 "outdoorCoil.refFlow.states[13].T", 1, 5, 5984, 0)
DeclareAlias2("outdoorCoil.refFlow.slipRatio.states[14].phase", "Phase of the fluid: 1 for 1-phase, 2 for two-phase, 0 for not known, e.g., interactive use [:#(type=Integer)]",\
 "outdoorCoil.refFlow.states[14].phase", 1, 5, 5985, 66)
DeclareAlias2("outdoorCoil.refFlow.slipRatio.states[14].p", "Pressure [Pa|bar]",\
 "outdoorCoil.refFlow.states[14].p", 1, 5, 5986, 0)
DeclareAlias2("outdoorCoil.refFlow.slipRatio.states[14].h", "Specific enthalpy [J/kg]",\
 "outdoorCoil.refFlow.states[14].h", 1, 5, 5987, 0)
DeclareAlias2("outdoorCoil.refFlow.slipRatio.states[14].d", "Density [kg/m3|g/cm3]",\
 "outdoorCoil.refFlow.states[14].d", 1, 5, 5988, 0)
DeclareAlias2("outdoorCoil.refFlow.slipRatio.states[14].T", "Temperature [K|degC]",\
 "outdoorCoil.refFlow.states[14].T", 1, 5, 5989, 0)
DeclareAlias2("outdoorCoil.refFlow.slipRatio.states[15].phase", "Phase of the fluid: 1 for 1-phase, 2 for two-phase, 0 for not known, e.g., interactive use [:#(type=Integer)]",\
 "outdoorCoil.refFlow.states[15].phase", 1, 5, 5990, 66)
DeclareAlias2("outdoorCoil.refFlow.slipRatio.states[15].p", "Pressure [Pa|bar]",\
 "outdoorCoil.refFlow.states[15].p", 1, 5, 5991, 0)
DeclareAlias2("outdoorCoil.refFlow.slipRatio.states[15].h", "Specific enthalpy [J/kg]",\
 "outdoorCoil.refFlow.states[15].h", 1, 5, 5992, 0)
DeclareAlias2("outdoorCoil.refFlow.slipRatio.states[15].d", "Density [kg/m3|g/cm3]",\
 "outdoorCoil.refFlow.states[15].d", 1, 5, 5993, 0)
DeclareAlias2("outdoorCoil.refFlow.slipRatio.states[15].T", "Temperature [K|degC]",\
 "outdoorCoil.refFlow.states[15].T", 1, 5, 5994, 0)
DeclareVariable("indoorCoil.refFlow.states[1].phase", "Phase of the fluid: 1 for 1-phase, 2 for two-phase, 0 for not known, e.g., interactive use [:#(type=Integer)]",\
 0, 0.0,2.0,0.0,0,644)
DeclareVariable("indoorCoil.refFlow.states[1].p", "Pressure [Pa|bar]", 1000000.0,\
 100000.0,4800000.0,4000000.0,0,512)
DeclareVariable("indoorCoil.refFlow.states[1].h", "Specific enthalpy [J/kg]", \
400000.0, 100000.0,500000.0,500000.0,0,512)
DeclareVariable("indoorCoil.refFlow.states[1].d", "Density [kg/m3|g/cm3]", 500, \
0.0,100000.0,500.0,0,512)
DeclareVariable("indoorCoil.refFlow.states[1].T", "Temperature [K|degC]", 298, \
1.0,10000.0,350.0,0,512)
DeclareVariable("indoorCoil.refFlow.states[2].phase", "Phase of the fluid: 1 for 1-phase, 2 for two-phase, 0 for not known, e.g., interactive use [:#(type=Integer)]",\
 0, 0.0,2.0,0.0,0,644)
DeclareVariable("indoorCoil.refFlow.states[2].p", "Pressure [Pa|bar]", 1000000.0,\
 100000.0,4800000.0,4000000.0,0,512)
DeclareVariable("indoorCoil.refFlow.states[2].h", "Specific enthalpy [J/kg]", \
400000.0, 100000.0,500000.0,500000.0,0,512)
DeclareVariable("indoorCoil.refFlow.states[2].d", "Density [kg/m3|g/cm3]", 500, \
0.0,100000.0,500.0,0,512)
DeclareVariable("indoorCoil.refFlow.states[2].T", "Temperature [K|degC]", 298, \
1.0,10000.0,350.0,0,512)
DeclareVariable("indoorCoil.refFlow.states[3].phase", "Phase of the fluid: 1 for 1-phase, 2 for two-phase, 0 for not known, e.g., interactive use [:#(type=Integer)]",\
 0, 0.0,2.0,0.0,0,644)
DeclareVariable("indoorCoil.refFlow.states[3].p", "Pressure [Pa|bar]", 1000000.0,\
 100000.0,4800000.0,4000000.0,0,512)
DeclareVariable("indoorCoil.refFlow.states[3].h", "Specific enthalpy [J/kg]", \
400000.0, 100000.0,500000.0,500000.0,0,512)
DeclareVariable("indoorCoil.refFlow.states[3].d", "Density [kg/m3|g/cm3]", 500, \
0.0,100000.0,500.0,0,512)
DeclareVariable("indoorCoil.refFlow.states[3].T", "Temperature [K|degC]", 298, \
1.0,10000.0,350.0,0,512)
DeclareVariable("indoorCoil.refFlow.states[4].phase", "Phase of the fluid: 1 for 1-phase, 2 for two-phase, 0 for not known, e.g., interactive use [:#(type=Integer)]",\
 0, 0.0,2.0,0.0,0,644)
DeclareVariable("indoorCoil.refFlow.states[4].p", "Pressure [Pa|bar]", 1000000.0,\
 100000.0,4800000.0,4000000.0,0,512)
DeclareVariable("indoorCoil.refFlow.states[4].h", "Specific enthalpy [J/kg]", \
400000.0, 100000.0,500000.0,500000.0,0,512)
DeclareVariable("indoorCoil.refFlow.states[4].d", "Density [kg/m3|g/cm3]", 500, \
0.0,100000.0,500.0,0,512)
DeclareVariable("indoorCoil.refFlow.states[4].T", "Temperature [K|degC]", 298, \
1.0,10000.0,350.0,0,512)
DeclareVariable("indoorCoil.refFlow.states[5].phase", "Phase of the fluid: 1 for 1-phase, 2 for two-phase, 0 for not known, e.g., interactive use [:#(type=Integer)]",\
 0, 0.0,2.0,0.0,0,644)
DeclareVariable("indoorCoil.refFlow.states[5].p", "Pressure [Pa|bar]", 1000000.0,\
 100000.0,4800000.0,4000000.0,0,512)
DeclareVariable("indoorCoil.refFlow.states[5].h", "Specific enthalpy [J/kg]", \
400000.0, 100000.0,500000.0,500000.0,0,512)
DeclareVariable("indoorCoil.refFlow.states[5].d", "Density [kg/m3|g/cm3]", 500, \
0.0,100000.0,500.0,0,512)
DeclareVariable("indoorCoil.refFlow.states[5].T", "Temperature [K|degC]", 298, \
1.0,10000.0,350.0,0,512)
DeclareVariable("indoorCoil.refFlow.states[6].phase", "Phase of the fluid: 1 for 1-phase, 2 for two-phase, 0 for not known, e.g., interactive use [:#(type=Integer)]",\
 0, 0.0,2.0,0.0,0,644)
DeclareVariable("indoorCoil.refFlow.states[6].p", "Pressure [Pa|bar]", 1000000.0,\
 100000.0,4800000.0,4000000.0,0,512)
DeclareVariable("indoorCoil.refFlow.states[6].h", "Specific enthalpy [J/kg]", \
400000.0, 100000.0,500000.0,500000.0,0,512)
DeclareVariable("indoorCoil.refFlow.states[6].d", "Density [kg/m3|g/cm3]", 500, \
0.0,100000.0,500.0,0,512)
DeclareVariable("indoorCoil.refFlow.states[6].T", "Temperature [K|degC]", 298, \
1.0,10000.0,350.0,0,512)
DeclareVariable("indoorCoil.refFlow.states[7].phase", "Phase of the fluid: 1 for 1-phase, 2 for two-phase, 0 for not known, e.g., interactive use [:#(type=Integer)]",\
 0, 0.0,2.0,0.0,0,644)
DeclareVariable("indoorCoil.refFlow.states[7].p", "Pressure [Pa|bar]", 1000000.0,\
 100000.0,4800000.0,4000000.0,0,512)
DeclareVariable("indoorCoil.refFlow.states[7].h", "Specific enthalpy [J/kg]", \
400000.0, 100000.0,500000.0,500000.0,0,512)
DeclareVariable("indoorCoil.refFlow.states[7].d", "Density [kg/m3|g/cm3]", 500, \
0.0,100000.0,500.0,0,512)
DeclareVariable("indoorCoil.refFlow.states[7].T", "Temperature [K|degC]", 298, \
1.0,10000.0,350.0,0,512)
DeclareVariable("indoorCoil.refFlow.states[8].phase", "Phase of the fluid: 1 for 1-phase, 2 for two-phase, 0 for not known, e.g., interactive use [:#(type=Integer)]",\
 0, 0.0,2.0,0.0,0,644)
DeclareVariable("indoorCoil.refFlow.states[8].p", "Pressure [Pa|bar]", 1000000.0,\
 100000.0,4800000.0,4000000.0,0,512)
DeclareVariable("indoorCoil.refFlow.states[8].h", "Specific enthalpy [J/kg]", \
400000.0, 100000.0,500000.0,500000.0,0,512)
DeclareVariable("indoorCoil.refFlow.states[8].d", "Density [kg/m3|g/cm3]", 500, \
0.0,100000.0,500.0,0,512)
DeclareVariable("indoorCoil.refFlow.states[8].T", "Temperature [K|degC]", 298, \
1.0,10000.0,350.0,0,512)
DeclareVariable("indoorCoil.refFlow.states[9].phase", "Phase of the fluid: 1 for 1-phase, 2 for two-phase, 0 for not known, e.g., interactive use [:#(type=Integer)]",\
 0, 0.0,2.0,0.0,0,644)
DeclareVariable("indoorCoil.refFlow.states[9].p", "Pressure [Pa|bar]", 1000000.0,\
 100000.0,4800000.0,4000000.0,0,512)
DeclareVariable("indoorCoil.refFlow.states[9].h", "Specific enthalpy [J/kg]", \
400000.0, 100000.0,500000.0,500000.0,0,512)
DeclareVariable("indoorCoil.refFlow.states[9].d", "Density [kg/m3|g/cm3]", 500, \
0.0,100000.0,500.0,0,512)
DeclareVariable("indoorCoil.refFlow.states[9].T", "Temperature [K|degC]", 298, \
1.0,10000.0,350.0,0,512)
DeclareVariable("indoorCoil.refFlow.states[10].phase", "Phase of the fluid: 1 for 1-phase, 2 for two-phase, 0 for not known, e.g., interactive use [:#(type=Integer)]",\
 0, 0.0,2.0,0.0,0,644)
DeclareVariable("indoorCoil.refFlow.states[10].p", "Pressure [Pa|bar]", \
1000000.0, 100000.0,4800000.0,4000000.0,0,512)
DeclareVariable("indoorCoil.refFlow.states[10].h", "Specific enthalpy [J/kg]", \
400000.0, 100000.0,500000.0,500000.0,0,512)
DeclareVariable("indoorCoil.refFlow.states[10].d", "Density [kg/m3|g/cm3]", 500,\
 0.0,100000.0,500.0,0,512)
DeclareVariable("indoorCoil.refFlow.states[10].T", "Temperature [K|degC]", 298, \
1.0,10000.0,350.0,0,512)
DeclareVariable("indoorCoil.refFlow.states[11].phase", "Phase of the fluid: 1 for 1-phase, 2 for two-phase, 0 for not known, e.g., interactive use [:#(type=Integer)]",\
 0, 0.0,2.0,0.0,0,644)
DeclareVariable("indoorCoil.refFlow.states[11].p", "Pressure [Pa|bar]", \
1000000.0, 100000.0,4800000.0,4000000.0,0,512)
DeclareVariable("indoorCoil.refFlow.states[11].h", "Specific enthalpy [J/kg]", \
400000.0, 100000.0,500000.0,500000.0,0,512)
DeclareVariable("indoorCoil.refFlow.states[11].d", "Density [kg/m3|g/cm3]", 500,\
 0.0,100000.0,500.0,0,512)
DeclareVariable("indoorCoil.refFlow.states[11].T", "Temperature [K|degC]", 298, \
1.0,10000.0,350.0,0,512)
DeclareVariable("indoorCoil.refFlow.states[12].phase", "Phase of the fluid: 1 for 1-phase, 2 for two-phase, 0 for not known, e.g., interactive use [:#(type=Integer)]",\
 0, 0.0,2.0,0.0,0,644)
DeclareVariable("indoorCoil.refFlow.states[12].p", "Pressure [Pa|bar]", \
1000000.0, 100000.0,4800000.0,4000000.0,0,512)
DeclareVariable("indoorCoil.refFlow.states[12].h", "Specific enthalpy [J/kg]", \
400000.0, 100000.0,500000.0,500000.0,0,512)
DeclareVariable("indoorCoil.refFlow.states[12].d", "Density [kg/m3|g/cm3]", 500,\
 0.0,100000.0,500.0,0,512)
DeclareVariable("indoorCoil.refFlow.states[12].T", "Temperature [K|degC]", 298, \
1.0,10000.0,350.0,0,512)
DeclareVariable("indoorCoil.refFlow.states[13].phase", "Phase of the fluid: 1 for 1-phase, 2 for two-phase, 0 for not known, e.g., interactive use [:#(type=Integer)]",\
 0, 0.0,2.0,0.0,0,644)
DeclareVariable("indoorCoil.refFlow.states[13].p", "Pressure [Pa|bar]", \
1000000.0, 100000.0,4800000.0,4000000.0,0,512)
DeclareVariable("indoorCoil.refFlow.states[13].h", "Specific enthalpy [J/kg]", \
400000.0, 100000.0,500000.0,500000.0,0,512)
DeclareVariable("indoorCoil.refFlow.states[13].d", "Density [kg/m3|g/cm3]", 500,\
 0.0,100000.0,500.0,0,512)
DeclareVariable("indoorCoil.refFlow.states[13].T", "Temperature [K|degC]", 298, \
1.0,10000.0,350.0,0,512)
DeclareVariable("indoorCoil.refFlow.states[14].phase", "Phase of the fluid: 1 for 1-phase, 2 for two-phase, 0 for not known, e.g., interactive use [:#(type=Integer)]",\
 0, 0.0,2.0,0.0,0,644)
DeclareVariable("indoorCoil.refFlow.states[14].p", "Pressure [Pa|bar]", \
1000000.0, 100000.0,4800000.0,4000000.0,0,512)
DeclareVariable("indoorCoil.refFlow.states[14].h", "Specific enthalpy [J/kg]", \
400000.0, 100000.0,500000.0,500000.0,0,512)
DeclareVariable("indoorCoil.refFlow.states[14].d", "Density [kg/m3|g/cm3]", 500,\
 0.0,100000.0,500.0,0,512)
DeclareVariable("indoorCoil.refFlow.states[14].T", "Temperature [K|degC]", 298, \
1.0,10000.0,350.0,0,512)
DeclareVariable("indoorCoil.refFlow.states[15].phase", "Phase of the fluid: 1 for 1-phase, 2 for two-phase, 0 for not known, e.g., interactive use [:#(type=Integer)]",\
 0, 0.0,2.0,0.0,0,644)
DeclareVariable("indoorCoil.refFlow.states[15].p", "Pressure [Pa|bar]", \
1000000.0, 100000.0,4800000.0,4000000.0,0,512)
DeclareVariable("indoorCoil.refFlow.states[15].h", "Specific enthalpy [J/kg]", \
400000.0, 100000.0,500000.0,500000.0,0,512)
DeclareVariable("indoorCoil.refFlow.states[15].d", "Density [kg/m3|g/cm3]", 500,\
 0.0,100000.0,500.0,0,512)
DeclareVariable("indoorCoil.refFlow.states[15].T", "Temperature [K|degC]", 298, \
1.0,10000.0,350.0,0,512)
DeclareAlias2("indoorCoil.refFlow.slipRatio.states[1].phase", "Phase of the fluid: 1 for 1-phase, 2 for two-phase, 0 for not known, e.g., interactive use [:#(type=Integer)]",\
 "indoorCoil.refFlow.states[1].phase", 1, 5, 5995, 66)
DeclareAlias2("indoorCoil.refFlow.slipRatio.states[1].p", "Pressure [Pa|bar]", \
"indoorCoil.refFlow.states[1].p", 1, 5, 5996, 0)
DeclareAlias2("indoorCoil.refFlow.slipRatio.states[1].h", "Specific enthalpy [J/kg]",\
 "indoorCoil.refFlow.states[1].h", 1, 5, 5997, 0)
DeclareAlias2("indoorCoil.refFlow.slipRatio.states[1].d", "Density [kg/m3|g/cm3]",\
 "indoorCoil.refFlow.states[1].d", 1, 5, 5998, 0)
DeclareAlias2("indoorCoil.refFlow.slipRatio.states[1].T", "Temperature [K|degC]",\
 "indoorCoil.refFlow.states[1].T", 1, 5, 5999, 0)
DeclareAlias2("indoorCoil.refFlow.slipRatio.states[2].phase", "Phase of the fluid: 1 for 1-phase, 2 for two-phase, 0 for not known, e.g., interactive use [:#(type=Integer)]",\
 "indoorCoil.refFlow.states[2].phase", 1, 5, 6000, 66)
DeclareAlias2("indoorCoil.refFlow.slipRatio.states[2].p", "Pressure [Pa|bar]", \
"indoorCoil.refFlow.states[2].p", 1, 5, 6001, 0)
DeclareAlias2("indoorCoil.refFlow.slipRatio.states[2].h", "Specific enthalpy [J/kg]",\
 "indoorCoil.refFlow.states[2].h", 1, 5, 6002, 0)
DeclareAlias2("indoorCoil.refFlow.slipRatio.states[2].d", "Density [kg/m3|g/cm3]",\
 "indoorCoil.refFlow.states[2].d", 1, 5, 6003, 0)
DeclareAlias2("indoorCoil.refFlow.slipRatio.states[2].T", "Temperature [K|degC]",\
 "indoorCoil.refFlow.states[2].T", 1, 5, 6004, 0)
DeclareAlias2("indoorCoil.refFlow.slipRatio.states[3].phase", "Phase of the fluid: 1 for 1-phase, 2 for two-phase, 0 for not known, e.g., interactive use [:#(type=Integer)]",\
 "indoorCoil.refFlow.states[3].phase", 1, 5, 6005, 66)
DeclareAlias2("indoorCoil.refFlow.slipRatio.states[3].p", "Pressure [Pa|bar]", \
"indoorCoil.refFlow.states[3].p", 1, 5, 6006, 0)
DeclareAlias2("indoorCoil.refFlow.slipRatio.states[3].h", "Specific enthalpy [J/kg]",\
 "indoorCoil.refFlow.states[3].h", 1, 5, 6007, 0)
DeclareAlias2("indoorCoil.refFlow.slipRatio.states[3].d", "Density [kg/m3|g/cm3]",\
 "indoorCoil.refFlow.states[3].d", 1, 5, 6008, 0)
DeclareAlias2("indoorCoil.refFlow.slipRatio.states[3].T", "Temperature [K|degC]",\
 "indoorCoil.refFlow.states[3].T", 1, 5, 6009, 0)
DeclareAlias2("indoorCoil.refFlow.slipRatio.states[4].phase", "Phase of the fluid: 1 for 1-phase, 2 for two-phase, 0 for not known, e.g., interactive use [:#(type=Integer)]",\
 "indoorCoil.refFlow.states[4].phase", 1, 5, 6010, 66)
DeclareAlias2("indoorCoil.refFlow.slipRatio.states[4].p", "Pressure [Pa|bar]", \
"indoorCoil.refFlow.states[4].p", 1, 5, 6011, 0)
DeclareAlias2("indoorCoil.refFlow.slipRatio.states[4].h", "Specific enthalpy [J/kg]",\
 "indoorCoil.refFlow.states[4].h", 1, 5, 6012, 0)
DeclareAlias2("indoorCoil.refFlow.slipRatio.states[4].d", "Density [kg/m3|g/cm3]",\
 "indoorCoil.refFlow.states[4].d", 1, 5, 6013, 0)
DeclareAlias2("indoorCoil.refFlow.slipRatio.states[4].T", "Temperature [K|degC]",\
 "indoorCoil.refFlow.states[4].T", 1, 5, 6014, 0)
DeclareAlias2("indoorCoil.refFlow.slipRatio.states[5].phase", "Phase of the fluid: 1 for 1-phase, 2 for two-phase, 0 for not known, e.g., interactive use [:#(type=Integer)]",\
 "indoorCoil.refFlow.states[5].phase", 1, 5, 6015, 66)
DeclareAlias2("indoorCoil.refFlow.slipRatio.states[5].p", "Pressure [Pa|bar]", \
"indoorCoil.refFlow.states[5].p", 1, 5, 6016, 0)
DeclareAlias2("indoorCoil.refFlow.slipRatio.states[5].h", "Specific enthalpy [J/kg]",\
 "indoorCoil.refFlow.states[5].h", 1, 5, 6017, 0)
DeclareAlias2("indoorCoil.refFlow.slipRatio.states[5].d", "Density [kg/m3|g/cm3]",\
 "indoorCoil.refFlow.states[5].d", 1, 5, 6018, 0)
DeclareAlias2("indoorCoil.refFlow.slipRatio.states[5].T", "Temperature [K|degC]",\
 "indoorCoil.refFlow.states[5].T", 1, 5, 6019, 0)
DeclareAlias2("indoorCoil.refFlow.slipRatio.states[6].phase", "Phase of the fluid: 1 for 1-phase, 2 for two-phase, 0 for not known, e.g., interactive use [:#(type=Integer)]",\
 "indoorCoil.refFlow.states[6].phase", 1, 5, 6020, 66)
DeclareAlias2("indoorCoil.refFlow.slipRatio.states[6].p", "Pressure [Pa|bar]", \
"indoorCoil.refFlow.states[6].p", 1, 5, 6021, 0)
DeclareAlias2("indoorCoil.refFlow.slipRatio.states[6].h", "Specific enthalpy [J/kg]",\
 "indoorCoil.refFlow.states[6].h", 1, 5, 6022, 0)
DeclareAlias2("indoorCoil.refFlow.slipRatio.states[6].d", "Density [kg/m3|g/cm3]",\
 "indoorCoil.refFlow.states[6].d", 1, 5, 6023, 0)
DeclareAlias2("indoorCoil.refFlow.slipRatio.states[6].T", "Temperature [K|degC]",\
 "indoorCoil.refFlow.states[6].T", 1, 5, 6024, 0)
DeclareAlias2("indoorCoil.refFlow.slipRatio.states[7].phase", "Phase of the fluid: 1 for 1-phase, 2 for two-phase, 0 for not known, e.g., interactive use [:#(type=Integer)]",\
 "indoorCoil.refFlow.states[7].phase", 1, 5, 6025, 66)
DeclareAlias2("indoorCoil.refFlow.slipRatio.states[7].p", "Pressure [Pa|bar]", \
"indoorCoil.refFlow.states[7].p", 1, 5, 6026, 0)
DeclareAlias2("indoorCoil.refFlow.slipRatio.states[7].h", "Specific enthalpy [J/kg]",\
 "indoorCoil.refFlow.states[7].h", 1, 5, 6027, 0)
DeclareAlias2("indoorCoil.refFlow.slipRatio.states[7].d", "Density [kg/m3|g/cm3]",\
 "indoorCoil.refFlow.states[7].d", 1, 5, 6028, 0)
DeclareAlias2("indoorCoil.refFlow.slipRatio.states[7].T", "Temperature [K|degC]",\
 "indoorCoil.refFlow.states[7].T", 1, 5, 6029, 0)
DeclareAlias2("indoorCoil.refFlow.slipRatio.states[8].phase", "Phase of the fluid: 1 for 1-phase, 2 for two-phase, 0 for not known, e.g., interactive use [:#(type=Integer)]",\
 "indoorCoil.refFlow.states[8].phase", 1, 5, 6030, 66)
DeclareAlias2("indoorCoil.refFlow.slipRatio.states[8].p", "Pressure [Pa|bar]", \
"indoorCoil.refFlow.states[8].p", 1, 5, 6031, 0)
DeclareAlias2("indoorCoil.refFlow.slipRatio.states[8].h", "Specific enthalpy [J/kg]",\
 "indoorCoil.refFlow.states[8].h", 1, 5, 6032, 0)
DeclareAlias2("indoorCoil.refFlow.slipRatio.states[8].d", "Density [kg/m3|g/cm3]",\
 "indoorCoil.refFlow.states[8].d", 1, 5, 6033, 0)
DeclareAlias2("indoorCoil.refFlow.slipRatio.states[8].T", "Temperature [K|degC]",\
 "indoorCoil.refFlow.states[8].T", 1, 5, 6034, 0)
DeclareAlias2("indoorCoil.refFlow.slipRatio.states[9].phase", "Phase of the fluid: 1 for 1-phase, 2 for two-phase, 0 for not known, e.g., interactive use [:#(type=Integer)]",\
 "indoorCoil.refFlow.states[9].phase", 1, 5, 6035, 66)
DeclareAlias2("indoorCoil.refFlow.slipRatio.states[9].p", "Pressure [Pa|bar]", \
"indoorCoil.refFlow.states[9].p", 1, 5, 6036, 0)
DeclareAlias2("indoorCoil.refFlow.slipRatio.states[9].h", "Specific enthalpy [J/kg]",\
 "indoorCoil.refFlow.states[9].h", 1, 5, 6037, 0)
DeclareAlias2("indoorCoil.refFlow.slipRatio.states[9].d", "Density [kg/m3|g/cm3]",\
 "indoorCoil.refFlow.states[9].d", 1, 5, 6038, 0)
DeclareAlias2("indoorCoil.refFlow.slipRatio.states[9].T", "Temperature [K|degC]",\
 "indoorCoil.refFlow.states[9].T", 1, 5, 6039, 0)
DeclareAlias2("indoorCoil.refFlow.slipRatio.states[10].phase", "Phase of the fluid: 1 for 1-phase, 2 for two-phase, 0 for not known, e.g., interactive use [:#(type=Integer)]",\
 "indoorCoil.refFlow.states[10].phase", 1, 5, 6040, 66)
DeclareAlias2("indoorCoil.refFlow.slipRatio.states[10].p", "Pressure [Pa|bar]", \
"indoorCoil.refFlow.states[10].p", 1, 5, 6041, 0)
DeclareAlias2("indoorCoil.refFlow.slipRatio.states[10].h", "Specific enthalpy [J/kg]",\
 "indoorCoil.refFlow.states[10].h", 1, 5, 6042, 0)
DeclareAlias2("indoorCoil.refFlow.slipRatio.states[10].d", "Density [kg/m3|g/cm3]",\
 "indoorCoil.refFlow.states[10].d", 1, 5, 6043, 0)
DeclareAlias2("indoorCoil.refFlow.slipRatio.states[10].T", "Temperature [K|degC]",\
 "indoorCoil.refFlow.states[10].T", 1, 5, 6044, 0)
DeclareAlias2("indoorCoil.refFlow.slipRatio.states[11].phase", "Phase of the fluid: 1 for 1-phase, 2 for two-phase, 0 for not known, e.g., interactive use [:#(type=Integer)]",\
 "indoorCoil.refFlow.states[11].phase", 1, 5, 6045, 66)
DeclareAlias2("indoorCoil.refFlow.slipRatio.states[11].p", "Pressure [Pa|bar]", \
"indoorCoil.refFlow.states[11].p", 1, 5, 6046, 0)
DeclareAlias2("indoorCoil.refFlow.slipRatio.states[11].h", "Specific enthalpy [J/kg]",\
 "indoorCoil.refFlow.states[11].h", 1, 5, 6047, 0)
DeclareAlias2("indoorCoil.refFlow.slipRatio.states[11].d", "Density [kg/m3|g/cm3]",\
 "indoorCoil.refFlow.states[11].d", 1, 5, 6048, 0)
DeclareAlias2("indoorCoil.refFlow.slipRatio.states[11].T", "Temperature [K|degC]",\
 "indoorCoil.refFlow.states[11].T", 1, 5, 6049, 0)
DeclareAlias2("indoorCoil.refFlow.slipRatio.states[12].phase", "Phase of the fluid: 1 for 1-phase, 2 for two-phase, 0 for not known, e.g., interactive use [:#(type=Integer)]",\
 "indoorCoil.refFlow.states[12].phase", 1, 5, 6050, 66)
DeclareAlias2("indoorCoil.refFlow.slipRatio.states[12].p", "Pressure [Pa|bar]", \
"indoorCoil.refFlow.states[12].p", 1, 5, 6051, 0)
DeclareAlias2("indoorCoil.refFlow.slipRatio.states[12].h", "Specific enthalpy [J/kg]",\
 "indoorCoil.refFlow.states[12].h", 1, 5, 6052, 0)
DeclareAlias2("indoorCoil.refFlow.slipRatio.states[12].d", "Density [kg/m3|g/cm3]",\
 "indoorCoil.refFlow.states[12].d", 1, 5, 6053, 0)
DeclareAlias2("indoorCoil.refFlow.slipRatio.states[12].T", "Temperature [K|degC]",\
 "indoorCoil.refFlow.states[12].T", 1, 5, 6054, 0)
DeclareAlias2("indoorCoil.refFlow.slipRatio.states[13].phase", "Phase of the fluid: 1 for 1-phase, 2 for two-phase, 0 for not known, e.g., interactive use [:#(type=Integer)]",\
 "indoorCoil.refFlow.states[13].phase", 1, 5, 6055, 66)
DeclareAlias2("indoorCoil.refFlow.slipRatio.states[13].p", "Pressure [Pa|bar]", \
"indoorCoil.refFlow.states[13].p", 1, 5, 6056, 0)
DeclareAlias2("indoorCoil.refFlow.slipRatio.states[13].h", "Specific enthalpy [J/kg]",\
 "indoorCoil.refFlow.states[13].h", 1, 5, 6057, 0)
DeclareAlias2("indoorCoil.refFlow.slipRatio.states[13].d", "Density [kg/m3|g/cm3]",\
 "indoorCoil.refFlow.states[13].d", 1, 5, 6058, 0)
DeclareAlias2("indoorCoil.refFlow.slipRatio.states[13].T", "Temperature [K|degC]",\
 "indoorCoil.refFlow.states[13].T", 1, 5, 6059, 0)
DeclareAlias2("indoorCoil.refFlow.slipRatio.states[14].phase", "Phase of the fluid: 1 for 1-phase, 2 for two-phase, 0 for not known, e.g., interactive use [:#(type=Integer)]",\
 "indoorCoil.refFlow.states[14].phase", 1, 5, 6060, 66)
DeclareAlias2("indoorCoil.refFlow.slipRatio.states[14].p", "Pressure [Pa|bar]", \
"indoorCoil.refFlow.states[14].p", 1, 5, 6061, 0)
DeclareAlias2("indoorCoil.refFlow.slipRatio.states[14].h", "Specific enthalpy [J/kg]",\
 "indoorCoil.refFlow.states[14].h", 1, 5, 6062, 0)
DeclareAlias2("indoorCoil.refFlow.slipRatio.states[14].d", "Density [kg/m3|g/cm3]",\
 "indoorCoil.refFlow.states[14].d", 1, 5, 6063, 0)
DeclareAlias2("indoorCoil.refFlow.slipRatio.states[14].T", "Temperature [K|degC]",\
 "indoorCoil.refFlow.states[14].T", 1, 5, 6064, 0)
DeclareAlias2("indoorCoil.refFlow.slipRatio.states[15].phase", "Phase of the fluid: 1 for 1-phase, 2 for two-phase, 0 for not known, e.g., interactive use [:#(type=Integer)]",\
 "indoorCoil.refFlow.states[15].phase", 1, 5, 6065, 66)
DeclareAlias2("indoorCoil.refFlow.slipRatio.states[15].p", "Pressure [Pa|bar]", \
"indoorCoil.refFlow.states[15].p", 1, 5, 6066, 0)
DeclareAlias2("indoorCoil.refFlow.slipRatio.states[15].h", "Specific enthalpy [J/kg]",\
 "indoorCoil.refFlow.states[15].h", 1, 5, 6067, 0)
DeclareAlias2("indoorCoil.refFlow.slipRatio.states[15].d", "Density [kg/m3|g/cm3]",\
 "indoorCoil.refFlow.states[15].d", 1, 5, 6068, 0)
DeclareAlias2("indoorCoil.refFlow.slipRatio.states[15].T", "Temperature [K|degC]",\
 "indoorCoil.refFlow.states[15].T", 1, 5, 6069, 0)
DeclareVariable("liquidLine.refFlow.states[1].phase", "Phase of the fluid: 1 for 1-phase, 2 for two-phase, 0 for not known, e.g., interactive use [:#(type=Integer)]",\
 0, 0.0,2.0,0.0,0,644)
DeclareVariable("liquidLine.refFlow.states[1].p", "Pressure [Pa|bar]", 1000000.0,\
 100000.0,4800000.0,4000000.0,0,512)
DeclareVariable("liquidLine.refFlow.states[1].h", "Specific enthalpy [J/kg]", \
400000.0, 100000.0,500000.0,500000.0,0,512)
DeclareVariable("liquidLine.refFlow.states[1].d", "Density [kg/m3|g/cm3]", 500, \
0.0,100000.0,500.0,0,512)
DeclareVariable("liquidLine.refFlow.states[1].T", "Temperature [K|degC]", 298, \
1.0,10000.0,350.0,0,512)
DeclareVariable("liquidLine.refFlow.states[2].phase", "Phase of the fluid: 1 for 1-phase, 2 for two-phase, 0 for not known, e.g., interactive use [:#(type=Integer)]",\
 0, 0.0,2.0,0.0,0,644)
DeclareVariable("liquidLine.refFlow.states[2].p", "Pressure [Pa|bar]", 1000000.0,\
 100000.0,4800000.0,4000000.0,0,512)
DeclareVariable("liquidLine.refFlow.states[2].h", "Specific enthalpy [J/kg]", \
400000.0, 100000.0,500000.0,500000.0,0,512)
DeclareVariable("liquidLine.refFlow.states[2].d", "Density [kg/m3|g/cm3]", 500, \
0.0,100000.0,500.0,0,512)
DeclareVariable("liquidLine.refFlow.states[2].T", "Temperature [K|degC]", 298, \
1.0,10000.0,350.0,0,512)
DeclareVariable("liquidLine.refFlow.states[3].phase", "Phase of the fluid: 1 for 1-phase, 2 for two-phase, 0 for not known, e.g., interactive use [:#(type=Integer)]",\
 0, 0.0,2.0,0.0,0,644)
DeclareVariable("liquidLine.refFlow.states[3].p", "Pressure [Pa|bar]", 1000000.0,\
 100000.0,4800000.0,4000000.0,0,512)
DeclareVariable("liquidLine.refFlow.states[3].h", "Specific enthalpy [J/kg]", \
400000.0, 100000.0,500000.0,500000.0,0,512)
DeclareVariable("liquidLine.refFlow.states[3].d", "Density [kg/m3|g/cm3]", 500, \
0.0,100000.0,500.0,0,512)
DeclareVariable("liquidLine.refFlow.states[3].T", "Temperature [K|degC]", 298, \
1.0,10000.0,350.0,0,512)
DeclareVariable("liquidLine.refFlow.states[4].phase", "Phase of the fluid: 1 for 1-phase, 2 for two-phase, 0 for not known, e.g., interactive use [:#(type=Integer)]",\
 0, 0.0,2.0,0.0,0,644)
DeclareVariable("liquidLine.refFlow.states[4].p", "Pressure [Pa|bar]", 1000000.0,\
 100000.0,4800000.0,4000000.0,0,512)
DeclareVariable("liquidLine.refFlow.states[4].h", "Specific enthalpy [J/kg]", \
400000.0, 100000.0,500000.0,500000.0,0,512)
DeclareVariable("liquidLine.refFlow.states[4].d", "Density [kg/m3|g/cm3]", 500, \
0.0,100000.0,500.0,0,512)
DeclareVariable("liquidLine.refFlow.states[4].T", "Temperature [K|degC]", 298, \
1.0,10000.0,350.0,0,512)
DeclareVariable("liquidLine.refFlow.states[5].phase", "Phase of the fluid: 1 for 1-phase, 2 for two-phase, 0 for not known, e.g., interactive use [:#(type=Integer)]",\
 0, 0.0,2.0,0.0,0,644)
DeclareVariable("liquidLine.refFlow.states[5].p", "Pressure [Pa|bar]", 1000000.0,\
 100000.0,4800000.0,4000000.0,0,512)
DeclareVariable("liquidLine.refFlow.states[5].h", "Specific enthalpy [J/kg]", \
400000.0, 100000.0,500000.0,500000.0,0,512)
DeclareVariable("liquidLine.refFlow.states[5].d", "Density [kg/m3|g/cm3]", 500, \
0.0,100000.0,500.0,0,512)
DeclareVariable("liquidLine.refFlow.states[5].T", "Temperature [K|degC]", 298, \
1.0,10000.0,350.0,0,512)
DeclareAlias2("liquidLine.refFlow.slipRatio.states[1].phase", "Phase of the fluid: 1 for 1-phase, 2 for two-phase, 0 for not known, e.g., interactive use [:#(type=Integer)]",\
 "liquidLine.refFlow.states[1].phase", 1, 5, 6070, 66)
DeclareAlias2("liquidLine.refFlow.slipRatio.states[1].p", "Pressure [Pa|bar]", \
"liquidLine.refFlow.states[1].p", 1, 5, 6071, 0)
DeclareAlias2("liquidLine.refFlow.slipRatio.states[1].h", "Specific enthalpy [J/kg]",\
 "liquidLine.refFlow.states[1].h", 1, 5, 6072, 0)
DeclareAlias2("liquidLine.refFlow.slipRatio.states[1].d", "Density [kg/m3|g/cm3]",\
 "liquidLine.refFlow.states[1].d", 1, 5, 6073, 0)
DeclareAlias2("liquidLine.refFlow.slipRatio.states[1].T", "Temperature [K|degC]",\
 "liquidLine.refFlow.states[1].T", 1, 5, 6074, 0)
DeclareAlias2("liquidLine.refFlow.slipRatio.states[2].phase", "Phase of the fluid: 1 for 1-phase, 2 for two-phase, 0 for not known, e.g., interactive use [:#(type=Integer)]",\
 "liquidLine.refFlow.states[2].phase", 1, 5, 6075, 66)
DeclareAlias2("liquidLine.refFlow.slipRatio.states[2].p", "Pressure [Pa|bar]", \
"liquidLine.refFlow.states[2].p", 1, 5, 6076, 0)
DeclareAlias2("liquidLine.refFlow.slipRatio.states[2].h", "Specific enthalpy [J/kg]",\
 "liquidLine.refFlow.states[2].h", 1, 5, 6077, 0)
DeclareAlias2("liquidLine.refFlow.slipRatio.states[2].d", "Density [kg/m3|g/cm3]",\
 "liquidLine.refFlow.states[2].d", 1, 5, 6078, 0)
DeclareAlias2("liquidLine.refFlow.slipRatio.states[2].T", "Temperature [K|degC]",\
 "liquidLine.refFlow.states[2].T", 1, 5, 6079, 0)
DeclareAlias2("liquidLine.refFlow.slipRatio.states[3].phase", "Phase of the fluid: 1 for 1-phase, 2 for two-phase, 0 for not known, e.g., interactive use [:#(type=Integer)]",\
 "liquidLine.refFlow.states[3].phase", 1, 5, 6080, 66)
DeclareAlias2("liquidLine.refFlow.slipRatio.states[3].p", "Pressure [Pa|bar]", \
"liquidLine.refFlow.states[3].p", 1, 5, 6081, 0)
DeclareAlias2("liquidLine.refFlow.slipRatio.states[3].h", "Specific enthalpy [J/kg]",\
 "liquidLine.refFlow.states[3].h", 1, 5, 6082, 0)
DeclareAlias2("liquidLine.refFlow.slipRatio.states[3].d", "Density [kg/m3|g/cm3]",\
 "liquidLine.refFlow.states[3].d", 1, 5, 6083, 0)
DeclareAlias2("liquidLine.refFlow.slipRatio.states[3].T", "Temperature [K|degC]",\
 "liquidLine.refFlow.states[3].T", 1, 5, 6084, 0)
DeclareAlias2("liquidLine.refFlow.slipRatio.states[4].phase", "Phase of the fluid: 1 for 1-phase, 2 for two-phase, 0 for not known, e.g., interactive use [:#(type=Integer)]",\
 "liquidLine.refFlow.states[4].phase", 1, 5, 6085, 66)
DeclareAlias2("liquidLine.refFlow.slipRatio.states[4].p", "Pressure [Pa|bar]", \
"liquidLine.refFlow.states[4].p", 1, 5, 6086, 0)
DeclareAlias2("liquidLine.refFlow.slipRatio.states[4].h", "Specific enthalpy [J/kg]",\
 "liquidLine.refFlow.states[4].h", 1, 5, 6087, 0)
DeclareAlias2("liquidLine.refFlow.slipRatio.states[4].d", "Density [kg/m3|g/cm3]",\
 "liquidLine.refFlow.states[4].d", 1, 5, 6088, 0)
DeclareAlias2("liquidLine.refFlow.slipRatio.states[4].T", "Temperature [K|degC]",\
 "liquidLine.refFlow.states[4].T", 1, 5, 6089, 0)
DeclareAlias2("liquidLine.refFlow.slipRatio.states[5].phase", "Phase of the fluid: 1 for 1-phase, 2 for two-phase, 0 for not known, e.g., interactive use [:#(type=Integer)]",\
 "liquidLine.refFlow.states[5].phase", 1, 5, 6090, 66)
DeclareAlias2("liquidLine.refFlow.slipRatio.states[5].p", "Pressure [Pa|bar]", \
"liquidLine.refFlow.states[5].p", 1, 5, 6091, 0)
DeclareAlias2("liquidLine.refFlow.slipRatio.states[5].h", "Specific enthalpy [J/kg]",\
 "liquidLine.refFlow.states[5].h", 1, 5, 6092, 0)
DeclareAlias2("liquidLine.refFlow.slipRatio.states[5].d", "Density [kg/m3|g/cm3]",\
 "liquidLine.refFlow.states[5].d", 1, 5, 6093, 0)
DeclareAlias2("liquidLine.refFlow.slipRatio.states[5].T", "Temperature [K|degC]",\
 "liquidLine.refFlow.states[5].T", 1, 5, 6094, 0)
DeclareVariable("vaporLine.refFlow.states[1].phase", "Phase of the fluid: 1 for 1-phase, 2 for two-phase, 0 for not known, e.g., interactive use [:#(type=Integer)]",\
 0, 0.0,2.0,0.0,0,644)
DeclareVariable("vaporLine.refFlow.states[1].p", "Pressure [Pa|bar]", 1000000.0,\
 100000.0,4800000.0,4000000.0,0,512)
DeclareVariable("vaporLine.refFlow.states[1].h", "Specific enthalpy [J/kg]", \
400000.0, 100000.0,500000.0,500000.0,0,512)
DeclareVariable("vaporLine.refFlow.states[1].d", "Density [kg/m3|g/cm3]", 500, \
0.0,100000.0,500.0,0,512)
DeclareVariable("vaporLine.refFlow.states[1].T", "Temperature [K|degC]", 298, \
1.0,10000.0,350.0,0,512)
DeclareVariable("vaporLine.refFlow.states[2].phase", "Phase of the fluid: 1 for 1-phase, 2 for two-phase, 0 for not known, e.g., interactive use [:#(type=Integer)]",\
 0, 0.0,2.0,0.0,0,644)
DeclareVariable("vaporLine.refFlow.states[2].p", "Pressure [Pa|bar]", 1000000.0,\
 100000.0,4800000.0,4000000.0,0,512)
DeclareVariable("vaporLine.refFlow.states[2].h", "Specific enthalpy [J/kg]", \
400000.0, 100000.0,500000.0,500000.0,0,512)
DeclareVariable("vaporLine.refFlow.states[2].d", "Density [kg/m3|g/cm3]", 500, \
0.0,100000.0,500.0,0,512)
DeclareVariable("vaporLine.refFlow.states[2].T", "Temperature [K|degC]", 298, \
1.0,10000.0,350.0,0,512)
DeclareVariable("vaporLine.refFlow.states[3].phase", "Phase of the fluid: 1 for 1-phase, 2 for two-phase, 0 for not known, e.g., interactive use [:#(type=Integer)]",\
 0, 0.0,2.0,0.0,0,644)
DeclareVariable("vaporLine.refFlow.states[3].p", "Pressure [Pa|bar]", 1000000.0,\
 100000.0,4800000.0,4000000.0,0,512)
DeclareVariable("vaporLine.refFlow.states[3].h", "Specific enthalpy [J/kg]", \
400000.0, 100000.0,500000.0,500000.0,0,512)
DeclareVariable("vaporLine.refFlow.states[3].d", "Density [kg/m3|g/cm3]", 500, \
0.0,100000.0,500.0,0,512)
DeclareVariable("vaporLine.refFlow.states[3].T", "Temperature [K|degC]", 298, \
1.0,10000.0,350.0,0,512)
DeclareVariable("vaporLine.refFlow.states[4].phase", "Phase of the fluid: 1 for 1-phase, 2 for two-phase, 0 for not known, e.g., interactive use [:#(type=Integer)]",\
 0, 0.0,2.0,0.0,0,644)
DeclareVariable("vaporLine.refFlow.states[4].p", "Pressure [Pa|bar]", 1000000.0,\
 100000.0,4800000.0,4000000.0,0,512)
DeclareVariable("vaporLine.refFlow.states[4].h", "Specific enthalpy [J/kg]", \
400000.0, 100000.0,500000.0,500000.0,0,512)
DeclareVariable("vaporLine.refFlow.states[4].d", "Density [kg/m3|g/cm3]", 500, \
0.0,100000.0,500.0,0,512)
DeclareVariable("vaporLine.refFlow.states[4].T", "Temperature [K|degC]", 298, \
1.0,10000.0,350.0,0,512)
DeclareVariable("vaporLine.refFlow.states[5].phase", "Phase of the fluid: 1 for 1-phase, 2 for two-phase, 0 for not known, e.g., interactive use [:#(type=Integer)]",\
 0, 0.0,2.0,0.0,0,644)
DeclareVariable("vaporLine.refFlow.states[5].p", "Pressure [Pa|bar]", 1000000.0,\
 100000.0,4800000.0,4000000.0,0,512)
DeclareVariable("vaporLine.refFlow.states[5].h", "Specific enthalpy [J/kg]", \
400000.0, 100000.0,500000.0,500000.0,0,512)
DeclareVariable("vaporLine.refFlow.states[5].d", "Density [kg/m3|g/cm3]", 500, \
0.0,100000.0,500.0,0,512)
DeclareVariable("vaporLine.refFlow.states[5].T", "Temperature [K|degC]", 298, \
1.0,10000.0,350.0,0,512)
DeclareAlias2("vaporLine.refFlow.slipRatio.states[1].phase", "Phase of the fluid: 1 for 1-phase, 2 for two-phase, 0 for not known, e.g., interactive use [:#(type=Integer)]",\
 "vaporLine.refFlow.states[1].phase", 1, 5, 6095, 66)
DeclareAlias2("vaporLine.refFlow.slipRatio.states[1].p", "Pressure [Pa|bar]", \
"vaporLine.refFlow.states[1].p", 1, 5, 6096, 0)
DeclareAlias2("vaporLine.refFlow.slipRatio.states[1].h", "Specific enthalpy [J/kg]",\
 "vaporLine.refFlow.states[1].h", 1, 5, 6097, 0)
DeclareAlias2("vaporLine.refFlow.slipRatio.states[1].d", "Density [kg/m3|g/cm3]",\
 "vaporLine.refFlow.states[1].d", 1, 5, 6098, 0)
DeclareAlias2("vaporLine.refFlow.slipRatio.states[1].T", "Temperature [K|degC]",\
 "vaporLine.refFlow.states[1].T", 1, 5, 6099, 0)
DeclareAlias2("vaporLine.refFlow.slipRatio.states[2].phase", "Phase of the fluid: 1 for 1-phase, 2 for two-phase, 0 for not known, e.g., interactive use [:#(type=Integer)]",\
 "vaporLine.refFlow.states[2].phase", 1, 5, 6100, 66)
DeclareAlias2("vaporLine.refFlow.slipRatio.states[2].p", "Pressure [Pa|bar]", \
"vaporLine.refFlow.states[2].p", 1, 5, 6101, 0)
DeclareAlias2("vaporLine.refFlow.slipRatio.states[2].h", "Specific enthalpy [J/kg]",\
 "vaporLine.refFlow.states[2].h", 1, 5, 6102, 0)
DeclareAlias2("vaporLine.refFlow.slipRatio.states[2].d", "Density [kg/m3|g/cm3]",\
 "vaporLine.refFlow.states[2].d", 1, 5, 6103, 0)
DeclareAlias2("vaporLine.refFlow.slipRatio.states[2].T", "Temperature [K|degC]",\
 "vaporLine.refFlow.states[2].T", 1, 5, 6104, 0)
DeclareAlias2("vaporLine.refFlow.slipRatio.states[3].phase", "Phase of the fluid: 1 for 1-phase, 2 for two-phase, 0 for not known, e.g., interactive use [:#(type=Integer)]",\
 "vaporLine.refFlow.states[3].phase", 1, 5, 6105, 66)
DeclareAlias2("vaporLine.refFlow.slipRatio.states[3].p", "Pressure [Pa|bar]", \
"vaporLine.refFlow.states[3].p", 1, 5, 6106, 0)
DeclareAlias2("vaporLine.refFlow.slipRatio.states[3].h", "Specific enthalpy [J/kg]",\
 "vaporLine.refFlow.states[3].h", 1, 5, 6107, 0)
DeclareAlias2("vaporLine.refFlow.slipRatio.states[3].d", "Density [kg/m3|g/cm3]",\
 "vaporLine.refFlow.states[3].d", 1, 5, 6108, 0)
DeclareAlias2("vaporLine.refFlow.slipRatio.states[3].T", "Temperature [K|degC]",\
 "vaporLine.refFlow.states[3].T", 1, 5, 6109, 0)
DeclareAlias2("vaporLine.refFlow.slipRatio.states[4].phase", "Phase of the fluid: 1 for 1-phase, 2 for two-phase, 0 for not known, e.g., interactive use [:#(type=Integer)]",\
 "vaporLine.refFlow.states[4].phase", 1, 5, 6110, 66)
DeclareAlias2("vaporLine.refFlow.slipRatio.states[4].p", "Pressure [Pa|bar]", \
"vaporLine.refFlow.states[4].p", 1, 5, 6111, 0)
DeclareAlias2("vaporLine.refFlow.slipRatio.states[4].h", "Specific enthalpy [J/kg]",\
 "vaporLine.refFlow.states[4].h", 1, 5, 6112, 0)
DeclareAlias2("vaporLine.refFlow.slipRatio.states[4].d", "Density [kg/m3|g/cm3]",\
 "vaporLine.refFlow.states[4].d", 1, 5, 6113, 0)
DeclareAlias2("vaporLine.refFlow.slipRatio.states[4].T", "Temperature [K|degC]",\
 "vaporLine.refFlow.states[4].T", 1, 5, 6114, 0)
DeclareAlias2("vaporLine.refFlow.slipRatio.states[5].phase", "Phase of the fluid: 1 for 1-phase, 2 for two-phase, 0 for not known, e.g., interactive use [:#(type=Integer)]",\
 "vaporLine.refFlow.states[5].phase", 1, 5, 6115, 66)
DeclareAlias2("vaporLine.refFlow.slipRatio.states[5].p", "Pressure [Pa|bar]", \
"vaporLine.refFlow.states[5].p", 1, 5, 6116, 0)
DeclareAlias2("vaporLine.refFlow.slipRatio.states[5].h", "Specific enthalpy [J/kg]",\
 "vaporLine.refFlow.states[5].h", 1, 5, 6117, 0)
DeclareAlias2("vaporLine.refFlow.slipRatio.states[5].d", "Density [kg/m3|g/cm3]",\
 "vaporLine.refFlow.states[5].d", 1, 5, 6118, 0)
DeclareAlias2("vaporLine.refFlow.slipRatio.states[5].T", "Temperature [K|degC]",\
 "vaporLine.refFlow.states[5].T", 1, 5, 6119, 0)
EndNonAlias(5)
PreNonAliasNew(6)
