model UCB_LC

  // UCB Local Comfort Model


  // A model for calculating the local thermal comfort levels in 16 body segments:
  // 1-head, 2-chest, 3-back, 4-pelvis, 5-lShoulder, 6-rShoulder, 7-lArm, 8-rArm,
  // 9-lHand, 10-rHand, 11-lThigh, 12-rThigh, 13-lLeg, 14-rLeg, 15-lFoot, 16-rFoot.

  // The model requires inputs from the UCB_LS (local sensation) and UCB_OS (overall sensation) sub-models.

  // Developed based on the work published in https://doi.org/10.1016/j.buildenv.2009.06.015
  // and https://escholarship.org/uc/item/11m0n1wt by Zhang et al.



  // Model Inputs
  Modelica.Blocks.Interfaces.RealInput LS[16]	annotation (Placement(transformation(extent={{-120,-20},{-80,20}})));	// Local Sensation
  Modelica.Blocks.Interfaces.RealInput OS	annotation (Placement(transformation(extent={{-120,40},{-80,80}})));	// Overall Sensation

  // Model Outputs
  Modelica.Blocks.Interfaces.RealOutput LC[16]	annotation (Placement(transformation(extent={{80,-20},{120,20}})));	// Local Comfort


  //Model Regression Coefficients
  Real c3[16,2] = {{0.35,0.35}, {-0.66,0.66}, {-0.45,0.45}, {0.59,0.00}, {0.30,0.35}, {0.30,0.35}, {-0.23,0.23}, {-0.23,0.23},
        	   {-0.80,0.80}, {-0.80,0.80}, {0.00,0.00}, {0.00,0.00}, {-0.20,0.61}, {-0.20,0.61}, {-0.91,0.40}, {-0.91,0.40}};
  Real c7[16,2] = {{0.28,0.40}, {0.39,0.90}, {0.96,0.00}, {0.50,0.00}, {0.00,0.00}, {0.00,0.00}, {0.00,0.71}, {0.00,0.71},
        	   {0.48,0.48}, {0.48,0.48}, {0.00,0.00}, {0.00,0.00}, {0.67,0.00}, {0.67,0.00}, {0.50,0.30}, {0.50,0.30}};
  Real C6[16] =   {2.17, 2.1, 2.1, 2.06, 2.14, 2.14, 2, 2, 1.98, 1.98, 1.98, 1.98, 2, 2, 2.13, 2.13};
  Real C8[16] =   {0.5, 0, 0, -0.51, -0.4, -0.4, -0.68, -0.68, 0, 0, 0, 0, 0, 0, 0, 0};
  Real n[16] = 	  {2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1.5, 1.5, 2, 2};

  Real C3[16]; // used to determine which c3 variable is required
  Real C7[16]; // used to determine which c7 variable is required
  //Note: changed c7[2,1] from 1.39 to 0.39, c7[7,2]&c7[8,2] from 1.71 to 0.71, and c7[13,1]&c7[14,1] from 1.67 to 0.67

  //Model Variables - used to simplify the code
  Real leftslope[16];
  Real rightslope[16];
  Real offset[16];
  Real maxcomfort[16];
  Real lc[16];



equation


  //determining which regression coefficients to use based on overall sensation preference
  for i in 1:16 loop
    if (OS < 0) then
      C7[i] = c7[i,1];
      C3[i] = c3[i,1];
    else
      C7[i] = c7[i,2];
      C3[i] = c3[i,2];
    end if;
  end for;


  // Model
  for i in 1:16 loop

    leftslope[i] =  (-4 - (C6[i] + C7[i]*abs(OS))) / ((abs(-4 + C3[i]*abs(OS) + C8[i]))^(n[i]));
    rightslope[i] = (-4 - (C6[i] + C7[i]*abs(OS))) / ((abs(4 + C3[i]*abs(OS) + C8[i]))^(n[i]));
    offset[i] =     (LS[i] + C3[i]*abs(OS) + C8[i]);
    maxcomfort[i] =  C6[i] + C7[i]*abs(OS);

    lc[i] = ((( (leftslope[i]-rightslope[i])/(exp(25*offset[i])+1))  + rightslope[i]) * ((abs(offset[i]))^(n[i]))) + maxcomfort[i];
    //Note: changed offset to be absolute value to correct model shape

    LC[i] = if (abs(lc[i])>4) then 4*(lc[i]/abs(lc[i])) else lc[i]; // limiting values to +-4

  end for;



end UCB_LC;
