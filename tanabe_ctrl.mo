record tanabe_ctrl

  // Record holder for variables within the Tanabe 65MN physiological model.


  // This record holds information for use in the Tanabe_65MN.mo model; specifically the sweat, shiver, vasodilation, and vasoconstriction control constants.
  // For information on their use see Tanabe_65MN.mo.

  // Developed by Conor O'Boyle (2024)

    extends Modelica.Icons.Record;

    parameter Real sw "sweat [W/K]";
    parameter Real ch "shiver [W/K]";
    parameter Real dl "vasodilation [1/K]";
    parameter Real st "vasoconstriction [1/K]";

end tanabe_ctrl;
