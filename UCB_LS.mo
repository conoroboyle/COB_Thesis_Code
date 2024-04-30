model UCB_LS

  // UCB Local Comfort Model


  // A model for calculating the local thermal sensation levels in 16 body segments:
  // 1-head, 2-chest, 3-back, 4-pelvis, 5-lShoulder, 6-rShoulder, 7-lArm, 8-rArm,
  // 9-lHand, 10-rHand, 11-lThigh, 12-rThigh, 13-lLeg, 14-rLeg, 15-lFoot, 16-rFoot,
  // 17-core (used by physiological model - not a UCB model output)

  // The model requires inputs from the Tanabe_65MN (physiological) model.

  // Developed based on the work published in https://doi.org/10.1016/j.buildenv.2009.06.018
  // and https://escholarship.org/uc/item/11m0n1wt by Zhang et al.

  

  // Model Inputs
  Modelica.Blocks.Interfaces.RealInput T_skin[17]	annotation (Placement(transformation(extent={{-120,-20},{-80,20}})));	// Local Skin/Core Temperature(s)

  // Model Outputs
  Modelica.Blocks.Interfaces.RealOutput LS[16]		annotation (Placement(transformation(extent={{80,-20},{120,20}})));	// Local Sensation


  //Model Regression Coefficients
  Real c1[16,2] = {{0.40,1.30}, {0.35,0.60}, {0.30,0.70}, {0.20,0.40}, {0.30,0.40}, {0.30,0.40}, {0.30,0.70}, {0.30,0.70},
                   {0.20,0.45}, {0.20,0.45}, {0.20,0.30}, {0.20,0.30}, {0.30,0.40}, {0.30,0.40}, {0.25,0.25}, {0.25,0.25}};
  Real c2[16,2] = {{543,90}, {39,136}, {88,192}, {75,137}, {156,167}, {156,167}, {144,125}, {144,125},
                   {19,46}, {19,46}, {151,263}, {151,263}, {206,212}, {206,212}, {109,162}, {109,162}};
  Real C3[16] = {-2289, -2135, -4054, -5053, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
  Real k1[16,2] = {{0.20,0.20}, {0.10,0.10}, {0.10,0.10}, {0.15,0.15}, {0.10,0.10}, {0.10,0.10}, {0.10,0.10}, {0.10,0.10},
                   {0.15,0.15}, {0.15,0.15}, {0.10,0.10}, {0.10,0.10}, {0.10,0.10}, {0.10,0.10}, {0.15,0.15}, {0.15,0.15}};
  Real C1[16]; // used to determine which C1 variable is required
  Real C2[16]; // used to determine which C2 variable is required
  Real K1[16]; // used to determine which K1 variable is required

  //Model Parameters
  parameter Real T_set[16] = {34.90,34.46,34.48,34.64,34.11,34.11,33.79,33.79,34.13,34.13,33.77,33.77,33.02,33.02,32.85,32.85}; // local skin neutral setpoint temperature - taken from 65MN model under neutral PMV conditions

  //Model Variables
  Real T_s[16];	// local skin temperature
  Real T_c;	// core temperature
  Real Ts_bar;	// mean skin temperature
  Real Tset_bar; // mean neutral skin setpoint temperature
  Real ls[16];	// local sensation



equation


  //temperature conversions & averages
  T_c = T_skin[17] - 273.15;
  for i in 1:16 loop
    T_s[i] = T_skin[i] - 273.15;
  end for;
  Ts_bar = sum(T_s)/16;
  Tset_bar = sum(T_set)/16;


  //determining which regression coefficients to use based on overall sensation preference
  for i in 1:16 loop    
    if (T_s[i] < T_set[i]) then
      C1[i] = c1[i,1];
      K1[i] = k1[i,1];
    else
      C1[i] = c1[i,2];
      K1[i] = k1[i,2];
    end if;

    C2[i] = if (der(T_s[i])<0) then c2[i,1] else c2[i,2];
  end for;


    // Model
  for i in 1:16 loop

    ls[i] = 4*(( 2 / (1 + exp(-C1[i]*(T_s[i] - T_set[i]) - K1[i]*((T_s[i] - Ts_bar) - (T_set[i] - Tset_bar))))) - 1) + (C2[i])*der(T_s[i]) + (C3[i])*der(T_c);
    
    LS[i] = if (abs(ls[i])>4) then 4*(ls[i]/abs(ls[i])) else ls[i]; // limiting values to +-4

  end for;



end UCB_LS;
