model UCB_OC

  // UCB Overall Comfort Model


  // A model for calculating the local overall comfort level based on local thermal comfort levels in 16 body segments:
  // 1-head, 2-chest, 3-back, 4-pelvis, 5-lShoulder, 6-rShoulder, 7-lArm, 8-rArm,
  // 9-lHand, 10-rHand, 11-lThigh, 12-rThigh, 13-lLeg, 14-rLeg, 15-lFoot, 16-rFoot.

  // The model requires inputs from the UCB_LC (local comfort) model.

  // Developed based on the work published in https://doi.org/10.1016/j.buildenv.2009.06.020
  // and https://escholarship.org/uc/item/11m0n1wt by Zhang et al.

  
  
  // Model Inputs
  Modelica.Blocks.Interfaces.RealInput LC[16  annotation (Placement(transformation(extent={{-120,-20},{-80,20}}))); // Local Comfort

  // Model Outputs
  Modelica.Blocks.Interfaces.RealOutput OC    annotation (Placement(transformation(extent={{80,-20},{120,20}})));   // Overall Comfort

  
  //Model Parameters
  parameter Boolean transient(start=true);  // are conditions transient? (Y=true, N=false)
  parameter Boolean control(start=false);   // does the occupant have control over their comfort? (Y=true, N=false)

  //Model Variables
  Real LC_reduced[14];   // list of local comfort levels with removal of duplicate extremity votes (hands/feet)
  Real LC_rank[14];      // used to sort reduced variable list in descending order (most comfortable -> most uncomfortable)
  Integer LC_rankI[14];  // used to track the body segment index for sorted (reduced) variables
   Real LC_sort[16];     // used to sort all local comfort levels in descending order (most comfortable -> most uncomfortable)
  Integer LC_sortI[16];  // used to track the body segment index for sorted (all) variables


  
equation

  // Rule 1 is only concerned with minimum values. If minimums contain both hands (9,10) or both feet (15,16) then only the minimum comfort level should be counted
  LC_reduced[9] = min(LC[9],LC[10]);    // removing duplicate hand votes
  LC_reduced[14] = min(LC[15],LC[16]);  // removing duplicate feet votes

  //creating new array with only 1 hand and 1 foot vote
  for i in 1:8 loop
    LC_reduced[i] = LC[i];
  end for;
  for i in 1:4 loop
    LC_reduced[9+i] = LC[10+i];
  end for;

  //sort new reduced local comfort array in order of most comfortable to least comfortable
  (LC_rank, LC_rankI) = Modelica.Math.Vectors.sort(LC_r1, ascending=false);

  // Rule 2  is used if the environment is transient or the subject has control.
  // Under asymmetric loading conditions the hands||feet may be both the highest and lowest values. A second sort function will catch any asymmetric loads from being filtered out.
  //sort original local comfort array in order of most comfortable to least comfortable
  (LC_sort,LC_sortI) = Modelica.Math.Vectors.sort(LC, ascending=false);


  // Model

  if transient or control then
    //Rule 2 - 2 worst votes and 1 best vote averaged
    OC = (LC_rank[14] + LC_rank[13] + LC_sort[1]) / 3;
  else
    //Rule 1 - 2 worst votes averaged
    OC = (LC_rank[14] + LC_rank[13]) / 2;
  end if;


end UCB_OC;
