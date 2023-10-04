within DynamicVCC.Components.Types;
type ModelStructure = enumeration(
    av_vb "volume cells connected to port_a, port_b",
    av_b "volume cell connected to port_a, flow cell connected to port_b",
    a_vb "flow cell connected to port_a, volume cell connected to port_b",
    a_v_b "flow cells connected to port_a, port_b")
      "Model structures of 1-D flow";
