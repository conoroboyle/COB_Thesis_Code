model UCB_OS

  // UCB Overall Sensation Model


  // A model for calculating the overall thermal sensation level based on the local thermal sensation levels in 16 body segments:
  // 1-head, 2-chest, 3-back, 4-pelvis, 5-lShoulder, 6-rShoulder, 7-lArm, 8-rArm,
  // 9-lHand, 10-rHand, 11-lThigh, 12-rThigh, 13-lLeg, 14-rLeg, 15-lFoot, 16-rFoot.

  // The model requires inputs from the UCB_LS (local sensation) sub-model.

  // Developed based on the work published in https://doi.org/10.1016/j.buildenv.2009.06.020
  // and https://escholarship.org/uc/item/11m0n1wt by Zhang et al.



  // Model Inputs
  Modelica.Blocks.Interfaces.RealInput LS[16]	      annotation (Placement(transformation(extent={{-120,-20},{-80,20}})));	// Local Sensation
  Modelica.Blocks.Interfaces.BooleanInput event[2]  annotation (Placement(transformation(extent={{-120,40},{-80,80}})));  // Event trigger
  //Note: the opposite-sensation model requires knowledge of when local heating/cooling loads are applied and removed. The event trigger is used to supply this information;
  //      see event.mo for more details.

  // Model Outputs
  Modelica.Blocks.Interfaces.RealOutput OS	        annotation (Placement(transformation(extent={{80,-20},{120,20}})));	  // Overall Sensation


  //Model Regression Coefficients
  Real a[16,3] = {{0.54,0.50,0.69}, {0.91,0.57,0.97}, {0.91,0.46,0.75}, {0.94,0.32,0.75}, {0.43,0.28,0.37}, {0.43,0.28,0.37}, {0.37,0.38,0.30}, {0.37,0.38,0.30},
                  {0.25,0.00,0.33}, {0.25,0.00,0.33}, {0.81,0.30,0.82}, {0.81,0.30,0.82}, {0.70,0.29,0.80}, {0.70,0.29,0.80}, {0.50,0.00,0.44}, {0.50,0.00,0.44}};
  Real b[16,3] = {{-1.1,0.00,1.10}, {-1.14,0.00,1.14}, {-0.92,0.00,0.92}, {-0.64,0.00,0.64}, {-0.56,0.00,0.56}, {-0.56,0.00,0.56}, {-0.73,0.00,0.73}, {-0.73,0.00,0.73},
                  {0.00,0.00,0.00}, {0.00,0.00,0.00}, {-0.60,0.00,0.60}, {-0.60,0.00,0.60}, {-0.59,0.00,0.59}, {-0.59,0.00,0.59}, {0.00,0.00,0.00}, {0.00,0.00,0.00}};
  Real c[16,3] = {{-2,0,2}, {-2,0,2}, {-2,0,2}, {-2,0,2}, {-2,0,2}, {-2,0,2}, {-2,0,2}, {-2,0,2}, {-2,0,2}, {-2,0,2}, {-2,0,2}, {-2,0,2}, {-2,0,2}, {-2,0,2}, {-2,0,1}, {-2,0,1}};
  
  //Model Variables
  //Sensation Groups
  Real plus[16];   // array contains 1 if that body segment index has a positive local sensation (LS>=0, warm) and a 0 if not
  Real minus[16];  // array contains 1 if that body segment index has a negative local sensation (LS<0, cool) and a 0 if not
  Real n_plus;     // total number of body segments with positive local sensations (sum(plus))
  Real n_minus;    // total number of body segments with negative local sensations (sum(minus))

  //Sensation Model Checks
  Boolean No_Opposite_Sensation_Model;  // apply only no-opposite-sensation model
  Boolean LTSwarm;  // conditions met to trigger low thermal sensation (warm) sub-model
  Boolean LTScool;  // conditions met to trigger low thermal sensation (cool) sub-model
  Boolean HTSwarm;  // conditions met to trigger high thermal sensation (warm) sub-model
  Boolean HTScool;  // conditions met to trigger high thermal sensation (cool) sub-model

  //No-Opposite-Sensation Rankings
  Real LS_NOS[14];    // array used to remove the least impactful local sensations from the hands&&feet for the n-o-s model calculations
  Real rank[14];      // used to sort the reduced local sensations in descending order (warmest -> coolest)
  Integer ranki[14];  // used to track the body segment index of the sorted array


  //Low Thermal Sensation Votes
  Integer np[1] = {14}; //total number of potential sensation votes for counting interval
  Real interval;        // counting interval
  Real LTS[14];         // qualifying local thermal sensations
  Real nLTS[14];        // number of qualifying local thermal sensations

  //Opposite-Sensation Model
  Real OS_big;  // main contribution to the overall sensation value (from n-o-s model)
  Real OS_mod;  // modifying contribution to the overall sensation value (from opposite sensations)

  //Local sensations prior to heating/cooling
  Real dLS[16](start={0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0});    // array used to save the state of the local sensations at the last event trigger
  Real dLSon[16](start={0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0});  // array used to save the state of the local sensations when a heating/cooling load is applied (event trigger 1)
  Real dLSoff[16](start={0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0}); // array used to save the state of the local sensations when a heating/cooling load is removed (event trigger 2)
  Real deltaLS[16];      // change in local sensation since the last event trigger
  Real force[16];        // potential individual force from each local sensation
  Real iF[16];           // individual force from opposing local sensations (removes non-qualifying values from force)
  Integer j[16];         // used to track which regression coefficients to apply (i.e. opposite=warm or opposite=cool)
  Real iF_sort[16];      // used to sort the individual forces in descending order
  Integer iF_index[16];  // used to track the body segment index for sorted individual forces
  Boolean OSMwarm;       // modifying contribution if opposite sensations are warm
  Boolean OSMcool;       // modifying contribution if opposite sensations are cool



equation

  //Creating Sensation Groups
  for i in 1:16 loop
    plus[i] = if (LS[i]>0) then 1 else 0;    // positive sensations (LS>0, warm)
    minus[i] = if (LS[i]<0) then 1 else 0;   // negative sensations (LS<0, cool)
  end for;
  n_plus = sum(plus);    // total number of positive sensations
  n_minus = sum(minus);  // total number of negative sensations


  //No-Opposite-Sensation Model trigger
  //Check 1: are all local sensation all warm or all cool?
  if ((n_plus<=0) or (n_minus<=0)) then
    No_Opposite_Sensation_Model = true;  // trigger no-opposite-sensation model
  //Check 2: if opposing sensations are present do they meet the threshold to trigger the opposite-sensation model?
  elseif (LTSwarm and (rank[14]>=-1)) or (LTScool and (rank[1]<=1)) then
    No_Opposite_Sensation_Model = true;  // trigger no-opposite-sensation model
  else
    No_Opposite_Sensation_Model = false;  // trigger opposite-sensation model
  end if;


  //Removal of duplicate votes for extremeties (hands/feet), highest impact value only is taken for n-o-s model
  for i in 1:8 loop
    LS_NOS[i] = LS[i];      // fill array with local sensations in all other body parts
  end for;
  for i in 1:4 loop
    LS_NOS[9+i] = LS[10+i]; // fill array with local sensations in all other body parts
  end for;
  if (n_plus>=n_minus) then
    LS_NOS[9] = max(LS[9],LS[10]);    // if overall sensation will be warm use highest (warmest) hand vote
    LS_NOS[14] = max(LS[15],LS[16]);  // if overall sensation will be warm use highest (warmest) foot vote
  else
    LS_NOS[9] = min(LS[9],LS[10]);    // if overall sensation will be cool use lowest (coolest) hand vote
    LS_NOS[14] = min(LS[15],LS[16]);  // if overall sensation will be cool use lowest (coolest) foot vote
  end if;
  (rank, ranki) = Modelica.Math.Vectors.sort(LS_NOS, ascending=false);  // sort votes in descending order (warmest -> coolest)

  //Low Thermal Sensation sub-model check
  LTSwarm = if (n_plus>=n_minus) then true else false;  // if more warm votes (LS>0) than cool votes, use LTS (warm) sub-model
  LTScool = if (n_plus<n_minus) then true else false;   // if more cool votes (LS<0) than warm votes, use LTS (cool) sub-model
  //High Thermal Sensation sub-model check
  HTSwarm = if (rank[3] >= 2) then true else false;     // if third highest local sensation meets HTS threshold (LS>=+2), use HTS (warm) sub-model
  HTScool = if (rank[12] <= -2) then true else false;   // if third lowest local sensation meets HTS threshold (LS<=-2), use HTS (cool) sub-model
  //Note: triggering the HTS sub-model will superceed the LTS sub-model

  //Number of local sensation votes to be counted in LTS sub-model
  if LTSwarm then         //Low Thermal Sensation (warm) sub-model triggered
    for i in 1:14 loop 
      if (i<=3) then      //Automatically consider the 3 warmest local sensation values
        LTS[i] = rank[i];    //add local sensation to array
        nLTS[i] = 1;         //add 1 to the number of body segments which qualify
      elseif (rank[i] > (2 - (interval*(i-3)))) then //Check remaining local sensations against validity criteria
        LTS[i] = rank[i];    //add local sensation to array if they qualify
        nLTS[i] = 1;         //add 1 to the number of body segments which qualify
      else
        LTS[i] = 0;          //add 0 to array elements which do not qualify
        nLTS[i] = 0;         //add 0 to the number of body segments which qualify
      end if;
    end for;

  else                    //Low Thermal Sensation (cool) sub-model triggered
    for i in 1:14 loop
      if (i<=3) then      //Automatically consider the 3 coolest local sensation values
        LTS[i] = rank[15-i]; //add local sensation to array
        nLTS[i] = 1;         //add 1 to the number of body segments which qualify
      elseif (rank[15-i] > (-2 + (interval*(i-3)))) then  //Check remaining local sensations against validity criteria
        LTS[i] = rank[15-i]; //add local sensation to array if they qualify
        nLTS[i] = 1;         //add 1 to the number of body segments which qualify
      else
        LTS[i] = 0;          //add 0 to array elements which do not qualify
        nLTS[i] = 0;         //add 0 to the number of body segments which qualify
      end if;
    end for;

  end if;
  interval = 2/(np[1]-2);  // counting interval for validity criteria


  //Overall Sensation Model
  if No_Opposite_Sensation_Model then      // if n-o-s model has been triggered

    if (HTSwarm or HTScool) then           // check if HTS sub-model has been triggered
      OS = if HTSwarm then 0.5*rank[1] + 0.5*rank[3] else 0.38*rank[14] + 0.62*rank[12];  // overall sensation is based on complaint model
    else                                   // otherwise apply LTS sub-model
      OS = sum(LTS)/sum(nLTS);             // overall sensation is the average of the qualifying local thermal sensations
    end if;
    //tie up loose ends (Modelica requires these variables to have an output or an error will be generated)
    OS_big = 0;
    OS_mod = 0;
    OSMwarm = false;
    OSMcool = false;

  else                                     // if n-o-s model only conditions are not met then apply opposite sensation conditions
                                             // determine which group of sensations are larger (warm sensations or cool sensations)
    if (n_plus>=n_minus) then                // if the warm group is larger it becomes the 'big' group / main contributor

      if ((LS[2]<=-1) or (LS[3]<=-1) or (LS[4]<=-1)) then //check for dominant cooling in the chest, back, or pelvis
        OS = min(LS[2],min(LS[3],LS[4]));                 //overall sensation will be the lowest (most cool) value if this rule is triggered
        //tie up loose ends
        OS_big = 0;
        OS_mod = 0;
        OSMwarm = false;
        OSMcool = false;
      
      else                                    //otherwise calculate main contribution using n-o-s methods
        if HTSwarm then                       //check for high thermal sensations
          OS_big = 0.5*rank[1] + 0.5*rank[3];   //if HTS then the main contribution to overall sensation is determined by HTS (warm) sub-model
        else
          OS_big = sum(LTS)/sum(nLTS);          //otherwise the main contribution to overall sensation is determined by LTS (warm) sub-model
        end if;

        OS_mod = iF_sort[16] + 0.1*iF_sort[15]; //modifying contribution is determined by the most negative individual forces
        OS = OS_big + OS_mod;                   //overall sensation is the combination of the main and modifying contributions

        //tie up loose ends
        OSMwarm = true;
        OSMcool = false;
      end if;

    else                                    // if the cool group is larger it becomes the 'big' group / main contributor
      if HTScool then                       //check for high thermal sensations
        OS_big = 0.38*rank[14] + 0.62*rank[12]; //if HTS then the main contribution to overall sensation is determined by HTS (cool) sub-model
      else
        OS_big = sum(LTS)/sum(nLTS);            //otherwise the main contribution to overall sensation is determined by LTS (warm) sub-model
      end if;

      OS_mod = iF_sort[1] + 0.1*iF_sort[2];     //modifying contribution is determined by the most positive individual forces
      OS = OS_big + OS_mod;                     //overall sensation is the combination of the main and modifying contributions

      //tie up loose ends
      OSMwarm = false;
      OSMcool = true;
    end if;
  end if;

  //event operations
  when event[1] then
    for i in 1:16 loop
      dLSon[i] = pre(LS[i]);
    end for;
  end when;
  when event[2] then
    for i in 1:16 loop
      dLSoff[i] = pre(LS[i]);
    end for;
  end when;
  //calculating individual forces
  for i in 1 : 16 loop
    dLS[i] = if event[1] then dLSon[i] else dLSoff[i];
    deltaLS[i] = LS[i] - dLS[i];
    force[i] = a[i,j[i]]*(deltaLS[i]-c[i,j[i]]) + b[i,j[i]];
    iF[i] = if ( (OSMwarm and (LS[i]<=0)) or (OSMcool and (LS[i]>=0)))  then force[i] else 0; //should that have been a +-2??? dont think it affects the results
    j[i] = if deltaLS[i]<=-2 then 1 elseif deltaLS[i]>=2 then 3 else 2;
  end for;
  (iF_sort, iF_index) = Modelica.Math.Vectors.sort(iF, ascending=false);


end UCB_OS;
