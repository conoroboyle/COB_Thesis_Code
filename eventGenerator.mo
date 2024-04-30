model eventGenerator

  // State event generator

  // A model for creating state events to trigger the pre() function allowing variables to be stored at desired timesteps/events.

  // The model requires a variable input [Real] which signals whether a heating/cooling load is being applied - 1=load-applied, 0=no-load//load-removed.
  // When the input (sig) passes a threshold value (0.5) the internal variables (on/off [Boolean]) will flip between true and false, this discontinuity will trigger a state event.
  // The state events allow the pre() function to output the value of any variable at the instant prior to the event occuring;
      // This is used to hold the previous local sensation values before a heating/cooling load is applied/removed to calculate how they have changed over time.

  // This model has been created to provide inputs to the UCB Overall Sensation model (UCB_OS.mo).
  // Developed by Conor O'Boyle (2024)


  
  // Model Inputs
  Modelica.Blocks.Interfaces.RealInput sig                         annotation (Placement(transformation(extent={{-120,-20},{-80,20}}))); // Input signal

  // Model Outputs
  Modelica.Blocks.Interfaces.BooleanOutput event[2] "1=on,2=off"   annotation (Placement(transformation(extent={{90,-10},{110,10}})));   // Event signals

  
  // Model Variables
  Boolean on;
  Boolean off;

  

equation

  
  on = if (sig > 0.5) then true else false;    // when the input signal is 1 then a heating/cooling application is applied (trigger event when crossing threshold)
  off = if (sig < 0.5) then true else false;   // when the input signal is 0 then a heating/cooling application is removed (trigger event when crossing threshold)

  event[1] = on;    // output event signal (on)
  event[2] = off;   // output event signal (off)
  //Note: 2 signals are required for the model to work due to its current setup. This will be reduced in future versions.


end eventGenerator;
