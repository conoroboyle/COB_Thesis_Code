model Tanabe_65MN

  // Tanabe 65-MultiNode (65MN) Physiological Model

  // A model for calculating the local physiological and thermoregulatory response to different environmental conditions in 16 body segments (represented by i in array indices):
    // i=1-head, i=2-chest, i=3-back, i=4-pelvis, i=5-lShoulder, i=6-rShoulder, i=7-lArm, i=8-rArm,
    // i=9-lHand, i=10-rHand, i=11-lThigh, i=12-rThigh, i=13-lLeg, i=14-rLeg, i=15-lFoot, i=16-rFoot;
  // And across 4 different types of body tissue (represented by j in array indices):
    // j=1-viscera, j=2-muscle, j=3-fat, j=4-skin;
  // With an additional node representing the response of the core pool of central blood in the body (node 65: labeled as '_cb' in the equations, and as index i=17 in model outputs).

  // The model requires inputs from high-fidelity numerical analysis to provide the local heat transfer coefficients, local air temperatures, and local radiative heat flux loads
  // for each body segment, as well as the bulk mean radiant temperature of the environment
  // This model also requires the tanabe_ctrl.mo record block to be available.

  // Developed based on the work published in https://doi.org/10.1016/S0378-7788(02)00014-2 by Tanabe et al.

  
  // Model Inputs
  Modelica.Blocks.Interfaces.RealInput T_air[16]              annotation (Placement(transformation(extent={{-120,-20},{-80,20}})));    // Local Air Temperature
  Modelica.Blocks.Interfaces.RealInput HTC[16]                annotation (Placement(transformation(extent={{-120,40},{-80,80}})));     // Local Heat Transfer Coefficient
  Modelica.Blocks.Interfaces.RealInput RAD[16]                annotation (Placement(transformation(extent={{-120,-80},{-80,-40}})));   // Local Radiative Heat Flux
  Modelica.Blocks.Interfaces.RealInput MRT                    annotation (Placement(transformation(extent={{-120,-110},{-80,-70}})));  // Mean Radiant Temperature
  Modelica.Blocks.Interfaces.RealInput Tamb                   annotation (Placement(transformation(extent={{-120,70},{-80,110}})));    // Bulk Air Temperature

  // Model Outputs
  Modelica.Blocks.Interfaces.RealOutput T_skin[17](unit="K")  annotation (Placement(transformation(extent={{80,-20},{120,20}})));      // Local Skin & Core Temperature(s)
  

  //Model Parameters
  parameter Real met = 1.0; //metabolic activity level [met]
  parameter Real I_cl[16] = {0, 0, 0, 0.34, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};  //local clothing insulation levels [clo]

  //Model Constants
  Real A_Du[16] = {0.14, 0.175, 0.161, 0.221, 0.096, 0.096, 0.063, 0.063, 0.05, 0.05, 0.209, 0.209, 0.112, 0.112, 0.056, 0.056};    //local DuBois body surface area {i} [m2]

  Real c[16,4] = {{2.576,0.386,0.258,0.282}, {2.915,5.669,1.496,0.418}, {2.471,5.022,1.322,0.386}, {6.017,7.997,2.102,0.606},
                  {0.503,1.078,0.207,0.151}, {0.503,1.078,0.207,0.151}, {0.321,0.681,0.131,0.099}, {0.321,0.681,0.131,0.099},
                  {0.082,0.037,0.052,0.099}, {0.082,0.037,0.052,0.099}, {1.665,3.604,0.560,0.423}, {1.665,3.604,0.560,0.423},
                  {0.793,1.715,0.268,0.204}, {0.793,1.715,0.268,0.204}, {0.139,0.037,0.077,0.125}, {0.139,0.037,0.077,0.125}};    //heat capacity in each node {i,j} [W.h/K]
  Real C[16,4];  //array used for conversion from [h]ours -> [s]econds

  Real Q_b[16,4] = {{16.843,0.217,0.109,0.131}, {21.182,2.537,0.568,0.179}, {18.699,2.537,0.501,0.158}, {8.05,4.067,0.804,0.254},
                    {0.181,0.423,0.610,0.050}, {0.181,0.423,0.610,0.050}, {0.094,0.220,0.031,0.026}, {0.094,0.220,0.031,0.026},
                    {0.045,0.022,0.023,0.050}, {0.045,0.022,0.023,0.050}, {0.343,0.824,0.151,0.122}, {0.343,0.824,0.151,0.122},
                    {0.102,0.220,0.035,0.023}, {0.102,0.220,0.035,0.023}, {0.122,0.035,0.056,0.100}, {0.122,0.035,0.056,0.100}};    //basal metabolic rate in each node {i,j} [W]

  Real bfb[16,4] = {{45,0.87,0.34,2.24}, {77.85,7.66,1.34,1.8}, {76.34,7.66,1.34,1.35}, {18.19,12.28,2.16,2.08},
                     {0.320,1.280,0.160,0.860}, {0.320,1.280,0.160,0.860}, {0.160,0.670,0.085,0.450}, {0.160,0.670,0.085,0.450},
                     {0.091,0.078,0.042,0.910}, {0.091,0.078,0.042,0.910}, {0.364,0.855,0.150,0.380}, {0.364,0.855,0.150,0.380},
                     {0.071,0.070,0.019,0.110}, {0.071,0.070,0.019,0.110}, {0.049,0.010,0.019,0.450}, {0.049,0.010,0.019,0.450}};    //basal blood flow rate in each node {i,j} [L/h]
  Real BFB[16,4];  //array used for conversion from [h]ours -> [s]econds

  Real C_d[16,3] = {{1.601,13.224,16.008}, {0.616,2.100,9.164}, {0.594,2.018,8.700}, {0.379,1.276,5.104},
                    {0.441,2.946,7.308}, {0.441,2.946,7.308}, {0.244,2.227,7.888}, {0.244,2.227,7.888},
                    {2.181,6.484,5.858}, {2.181,6.484,5.858}, {2.401,4.536,30.160}, {2.401,4.536,30.160},
                    {1.891,2.656,7.540}, {1.891,2.656,7.540}, {8.12,10.266,8.178}, {8.12,10.266,8.178}};    //thermal conductance between node {i,j} and its neighbour {i,j+1} [W/K]

  Real T_set[16,4] = {{36.9,36.1,35.8,35.6}, {36.5,36.2,34.5,33.6}, {36.5,35.8,34.4,33.2}, {36.3,35.6,34.5,33.4},
                      {35.8,34.6,33.8,33.4}, {35.8,34.6,33.8,33.4}, {35.5,34.8,34.7,34.6}, {35.5,34.8,34.7,34.6},
                      {35.4,35.3,35.3,35.2}, {35.4,35.3,35.3,35.2}, {35.8,35.2,34.4,33.8}, {35.8,35.2,34.4,33.8},
                      {35.6,34.4,33.9,33.4}, {35.6,34.4,33.9,33.4}, {35.1,34.9,34.4,33.9}, {35.1,34.9,34.4,33.9}};

  Real SKINR[16] = {0.07, 0.149, 0.132, 0.212, 0.023, 0.023, 0.012, 0.012, 0.092, 0.092, 0.05, 0.05, 0.025, 0.025, 0.017, 0.017};    //distribution coefficient for skin sensors
  Real SKINS[16] = {0.081, 0.146, 0.129, 0.206, 0.051, 0.051, 0.026, 0.026, 0.016, 0.016, 0.073, 0.073, 0.036, 0.036, 0.018, 0.018}; //distribution coefficient for skin sweat
  Real SKINV[16] = {0.32, 0.098, 0.086, 0.138, 0.031, 0.031, 0.016, 0.016, 0.061, 0.061, 0.092, 0.092, 0.023, 0.023, 0.05, 0.05};    //distribution coefficient for skin vasodilation
  Real SKINC[16] = {0.022, 0.065, 0.065, 0.065, 0.022, 0.022, 0.022, 0.022, 0.152, 0.152, 0.022, 0.022, 0.022, 0.022, 0.152, 0.152}; //distribution coefficient for skin vasoconstriction
  Real Chilf[16] = {0.02, 0.258, 0.227, 0.365, 0.004, 0.004, 0.026, 0.026, 0, 0, 0.023, 0.023, 0.012, 0.012, 0, 0};                  //distribution coefficient for muscle heat production (shiver)
  Real Metf[16] = {0, 0.091, 0.08, 0.129, 0.026, 0.026, 0.014, 0.014, 0.005, 0.005, 0.201, 0.201, 0.099, 0.099, 0.005, 0.005};       //distribution coefficient for muscle heat production (work)
  
  Real SHIV;
  
  Real C_cb = 2.610*3600;
  Real arC = 1.067;
  Real RH = 70;
  Real T_d[16];
  Real p_d[16];
  Real T0;

  Real T_cb(start = 37.11); //t=60
//  Real T_cb(start = 36.23); //t=240
  Real res;
  Real WRMS;
  Real CLDS;
  //Real LR = 2.2; //old value
  Real LR = 16.5; //JOS3 value
  Real dl;
  Real DL;
  Real st;
  Real ST;
  Real RT = 10;
  Real Qb = 0.778; //basal metabolic output for whole body

  Real T_a[16];
  Real T_s[17];

  Real T_ave;
  Real E_v;
  Real T_ave2;

  Real B[16,4];
  Real BF[16,4];
  Real ch[16,4];
  Real C_h[16,4];
  Real Cld[16,4];
  Real clds[16];
  Real D[16,4];
  Real e[16];
  Real E[16];
  Real E_b[16];
  Real esw[16];
  Real E_sw[16];
  Real E_max[16];
  Real Err[16,4];
  Real f_cl[16];
  //Real f_cl[16] = {1.57, 1.04, 1.16, 1.05, 1.23, 1.23, 1.25, 1.25, 1, 1, 1.06, 1.06, 1.11, 1.11, 1.38, 1.38};
  //Real h_c[16] = {4.0, 4.2, 3.8, 4.2, 5.2, 5.2, 5, 5, 6, 6, 4.3, 4.3, 4.1, 4.1, 6, 6};
  //Real h_r[16] = {4.6, 3.7, 4.3, 4, 4.6, 4.6, 3.9, 3.9, 4.4, 4.4, 3.8, 3.8, 4.4, 4.4, 5.9, 5.9};
  Real h_c[16];
  Real h_r[16];
  Real h_e[16];
  Real h_t[16];
  Real h[16];
  Real i_cl[16];
  Real km[16];
  Real p_a[16];
  Real p_sks[16];
  Real Q[16,4];
  Real Q_t[16];
  Real RES[16];
  Real RATE[16,4];
  Real T[16,4]; //(start = {{36.9,36.1,35.8,35.6}, {36.5,36.2,34.5,33.6}, {36.5,35.8,34.4,33.2}, {36.3,35.6,34.5,33.4},
                //          {35.8,34.6,33.8,33.4}, {35.8,34.6,33.8,33.4}, {35.5,34.8,34.7,34.6}, {35.5,34.8,34.7,34.6},
                //          {35.4,35.3,35.3,35.2}, {35.4,35.3,35.3,35.2}, {35.8,35.2,34.4,33.8}, {35.8,35.2,34.4,33.8},
                //          {35.6,34.4,33.9,33.4}, {35.6,34.4,33.9,33.4}, {35.1,34.9,34.4,33.9}, {35.1,34.9,34.4,33.9}});
  Real t_0[16];
  Real w[16,4];
  Real W[16,4];
  Real Wrm[16,4];
  Real wrms[16];

  
  Real Mshiv[16];

  Real sumQ;
  Real sumB;
  Real sumW;
  Real sumCh;
  Real sumQt;
  Real sumE;
  Real sumEsw;

  Real T_o[16];

  tanabe_ctrl Core(sw=371.2,ch=0,dl=117,st=11.5)  annotation (Placement(transformation(extent={{-30,-70},{-10,-50}})));
  tanabe_ctrl Skin(sw=33.6,ch=0,dl=7.5,st=11.5)   annotation (Placement(transformation(extent={{10,-70},{30,-50}})));
  tanabe_ctrl P(sw=0,ch=24.4,dl=0,st=0)           annotation (Placement(transformation(extent={{50,-70},{70,-50}})));

// run 1 - start = 60
  Real T_init[16,4] = {{36.47, 34.6, 34.16, 33.75}, {36.27, 36.58, 34.65, 33.92}, {36.41, 36.56, 34.61, 33.85}, {36.58, 36.67, 33.68, 31.63},
                        {36.01, 35.61, 34.8, 34.36}, {36.01, 35.61, 34.8, 34.36}, {35.89, 35.47, 34.4, 34.08}, {35.88, 35.45, 34.38, 34.05},
                        {30.19, 29.9, 29.58, 29.18}, {29.92, 29.62, 29.29, 28.88}, {35.37, 35.04, 32, 31.51}, {35.29, 34.93, 31.73, 31.21},
                        {34.33, 34.12, 30.85, 29.67}, {34.28, 34.06, 30.71, 29.51}, {30.16, 30.09, 29.95, 29.73}, {29.74, 29.67, 29.52, 29.29}}; //for 17deg at 1.3clo at crash point

// run 2 - start = 240
//  Real T_init[16,4] = { {36.530077,34.5789055,34.0849915,33.588355}, {36.324846,36.593997,34.6445165,33.7975265}, {36.482216,36.574137,34.607044,33.778338}, {36.600231,36.6844985,33.670663,31.3624265},
//   {36.020325,35.6213725,34.755017,34.2572595}, {36.020325,35.6213725,34.755017,34.2572595}, {35.899527,35.473067,34.3341275,33.9589255}, {35.8894225,35.452642,34.279998,33.8824445},
//   {30.168994,29.8204325,29.4748175,29.01317}, {29.895504,29.5282315,29.1642955,28.671571}, {35.3677425,35.0413535,31.784742,31.1758225}, {35.2865585,34.931363,31.522199,30.886781},
//   {34.3215215,34.1180015,30.729245,29.380078}, {34.27057,34.057146,30.6005365,29.2456255}, {30.1308955,30.064905,29.9221805,29.7150395}, {29.7158605,29.650447,29.5033695,29.2883225}};

//  Real T_init[16,4] = {{36.9,36.1,35.8,35.6}, {36.5,36.2,34.5,33.6}, {36.5,35.8,34.4,33.2}, {36.3,35.6,34.5,33.4},
//                      {35.8,34.6,33.8,33.4}, {35.8,34.6,33.8,33.4}, {35.5,34.8,34.7,34.6}, {35.5,34.8,34.7,34.6},
//                      {35.4,35.3,35.3,35.2}, {35.4,35.3,35.3,35.2}, {35.8,35.2,34.4,33.8}, {35.8,35.2,34.4,33.8},
//                      {35.6,34.4,33.9,33.4}, {35.6,34.4,33.9,33.4}, {35.1,34.9,34.4,33.9}, {35.1,34.9,34.4,33.9}};

//  Real T_init[16,4] = {{36.3801,34.369923,33.840313,33.34496},{36.193737,36.55834,34.013107,33.04382},{36.327477,36.53142,34.004562,33.020077},{36.50928,36.656628,32.89653,30.320864},
//                        {35.686703,34.980026,33.681522,33.01845},{35.686703,34.980026,33.681522,33.01845},{35.77191,35.275238,33.824665,33.38565},{35.77191,35.275238,33.824665,33.38565},
//                        {29.107143,28.769588,28.413647,27.955654},{29.107143,28.769588,28.413647,27.955654},{35.15843,34.764465,31.245749,30.684666},{35.15843,34.764465,31.245749,30.684666},
//                        {34.359196,34.166264,30.854486,29.66598},{34.359196,34.166264,30.854486,29.66598},{28.433496,28.355068,28.183147,27.93326},{28.433496,28.355068,28.183147,27.93326}}; //for 17deg at 0.65clo


//  Real T_init[16,4] = {{36.183342,33.48974,32.49579,31.590584},{35.9551,36.516846,31.82774,30.043242},{36.154015,36.48825,31.854177,30.047497},{36.348495,36.667465,31.55319,28.082823},
//                       {34.79925,33.31224,30.536646,29.206947},{34.79925,33.31224,30.536646,29.206947},{35.500713,34.881542,31.851164,30.945454},{35.500713,34.881542,31.851164,30.945454},
//                       {26.288668,25.835714,25.358637,24.745087},{26.288668,25.835714,25.358637,24.745087},{34.968018,34.54971,30.271544,29.592018},{34.968018,34.54971,30.271544,29.592018},
//                       {34.21073,34.02398,30.017149,28.582853},{34.21073,34.02398,30.017149,28.582853},{22.732122,22.625874,22.367723,21.999928},{22.732122,22.625874,22.367723,21.999928}}; //for 17deg at 0.0clo



  
initial equation

  for i in 1:16 loop
    for j in 1:4 loop
      T[i,j] = T_init[i,j];
      //T_init[i,j] = T_air[i];
    end for;
  end for;

equation

  sumQ = sum(Q);
  sumB = sum(B);
  sumW = sum(W);
  sumCh = sum(C_h);
  sumQt = sum(Q_t);
  sumE = sum(E);
  sumEsw = sum(E_sw);

  //heat balance blood
  C_cb * der(T_cb) = sum(B);

  SHIV = 24.36*Cld[1,1]*CLDS;

  for i in 1:16 loop

    T_a[i] = T_air[i] - 273.15;
    h[i] = h_c[i] + h_r[i];
//    h_r[i] = if (RAD[i]>10) then abs(RAD[i] / abs(T_s[i] - MRT)) else 4; //T_s is the one in Celcius
//    h_r[i] = if (RAD[i]>10) then abs(RAD[i] / abs(T_skin[i] - MRT)) else 4; //T_skin is the one in Kelvin -- USE THIS ONE!!!
    h_r[i] = RAD[i]; //only use for paper yoke
//    h_r[i] = 0.96 * (5.670374e-8) * ((T_skin[i]^2)+(MRT^2))* (T_skin[i]^2+MRT);
    h_c[i] = if (HTC[i]>0) then HTC[i] else 2;
    T_d[i] = T_a[i] - ((100-RH)/5);
    p_d[i] = 0.61078 * exp((17.2694*(T_d[i]))/((T_d[i])+238.3));

    //heat balance
    C[i,1] * der(T[i,1]) = Q[i,1] - B[i,1] - D[i,1] - RES[i];
    C[i,2] * der(T[i,2]) = Q[i,2] - B[i,2] + D[i,1] - D[i,2];
    C[i,3] * der(T[i,3]) = Q[i,3] - B[i,3] + D[i,2] - D[i,3];
    C[i,4] * der(T[i,4]) = Q[i,4] - B[i,4] + D[i,3] - Q_t[i] - E[i];

    //respiration
    RES[i] = if (i==2) then res else 0;

    //evaporation
    e[i] = E_b[i] + E_sw[i];
    E[i] = if (e[i]>E_max[i]) then E_max[i] else e[i];
    E_b[i] = 0.06 * ((1-E_sw[i])/E_max[i]) * E_max[i];
    E_max[i] = h_e[i] * (p_sks[i] - p_d[i]) * A_Du[i];
    h_e[i] = (LR * i_cl[i])/((0.155*I_cl[i]) + (i_cl[i]/(h_c[i]*f_cl[i])));

///////////////////NEED VALUES///////////////////////////////
//    h_c[i] = 4; /////////////
//    I_cl[i] = 0.6; //////////////
    i_cl[i] = 0.4; ////////////
    f_cl[i] = 1 + 0.3*I_cl[i]; ///////////
//    h_r[i] = 4; /////////////
    t_0[i] = T_a[i]; /////////////
    T_o[i] = (h_r[i]*(MRT-273.15) + h_c[i]*T_a[i])/(h_r[i]+h_c[i]);
///////////////////NEED VALUES///////////////////////////////

    //sensible
//    Q_t[i] = h_t[i] * (T[i,4] - T_o[i]) * A_Du[i]; //technically correct but sim is using other eqn
    Q_t[i] = h_t[i] * (T[i,4] - t_0[i]) * A_Du[i]; //use this eqn
    1/(h_t[i]) = 0.155*I_cl[i] + (1/(h_c[i] + h_r[i]*f_cl[i]));

    //signal
    wrms[i] = SKINR[i] * Wrm[i,4];
    clds[i] = SKINR[i] * Cld[i,4];

    //local vasomotion factor
    km[i] = 2.0^(Err[i,4]/RT);

    //sweat
    esw[i] = (Core.sw*Err[1,1] + Skin.sw*(WRMS-CLDS) + P.sw*Wrm[1,1]*WRMS) * SKINS[i] * km[i];
    E_sw[i] = if (esw[i]>0) then esw[i] else 0;

//////////////CHECK THESE VALUES/////////////////////////////
//    LR[i] = 15.15 * ((T[i,4] + 273.2)/273.2);
    p_sks[i] = 0.61078 * exp((17.2694*(T[i,4]))/((T[i,4])+238.3));
    p_a[i] = 0.61078 * exp((17.2694*(T_a[i]))/((T_a[i])+238.3));

    Mshiv[i] = SHIV*Chilf[i];

      for j in 1:4 loop

        RATE[i,j] = 0*3600;
        BFB[i,j] = bfb[i,j];
        C[i,j] = c[i,j]*3600;

        //heat production
        Q[i,j] = Q_b[i,j] + W[i,j] + C_h[i,j];
        w[i,j] = if (j == 2) then 58.2 * (met - Qb) * sum(A_Du) * Metf[i] else 0;
        //w[i,j] = if (j == 2) then 58.2 * (MET - Qb) * sum(A_Du) * Metf[i] else 0;
        W[i,j] = if (w[i,j] > 0) then w[i,j] else 0;

        //blood flow
        B[i,j] = (arC) * BF[i,j] * (T[i,j] - T_cb);
        if (j==4) then
          BF[i,j] = ((BFB[i,j] + (SKINV[i]*DL)) / (1 + (SKINC[i]*ST))) * km[i];
        else
          BF[i,j] = BFB[i,j] + ((W[i,j] + C_h[i,j])/1.16);
        end if;

        //conduction
        D[i,j] = if (j<4) then C_d[i,j] * (T[i,j] - T[i,(j+1)]) else 0;

        //error
        Err[i,j] = (T[i,j] - T_set[i,j]) + RATE[i,j]*der(T[i,j]);
        Wrm[i,j] = if (Err[i,j]>0) then abs(Err[i,j]) else 0;
        Cld[i,j] = if (Err[i,j]<0) then abs(Err[i,j]) else 0;

        //shiver
        ch[i,j] = if (j==2) then (-Core.ch*Err[1,1] - Skin.ch*(WRMS-CLDS) + P.ch*Cld[1,1]*CLDS)*Chilf[i] else 0;
        C_h[i,j] = if (ch[i,j]>0) then ch[i,j] else 0;

      end for;

      T_s[i] = T[i,4];
      T_skin[i] = T_s[i] + 273.15;

      //f_cl[i] = (tcl1[i]-T_a[i])/(T_s[i]-T_a[i]);
      //Cv[i] + R[i] = (T_s[i]-tcl2[i])/0.155*I_cl[i];
      //Cv[i] = HTC[i]*(tcl2[i]-T_a[i]);
      //R[i] = (RAD[i]/(T_skin[i]-MRT))*(tcl2[i]-MRT+273.15);


  end for;

  //respiration loss
  res = (0.0014*(34 - T_a[1]) + 0.017*(5.867 - p_d[1])) * sum(Q);

  //warm & cold signal
  WRMS = sum(wrms);
  CLDS = sum(clds);

  //vasomotion signal
  dl = Core.dl * Err[1,1] + Skin.dl * (WRMS - CLDS) + P.dl*Wrm[1,1]*WRMS;
  st = -Core.st * Err[1,1] - Skin.st * (WRMS - CLDS) + P.st*Cld[1,1]*CLDS;
  DL = if (dl>0) then dl else 0;
  ST = if (st>0) then st else 0;

  T_s[17] = T_cb;
  T_skin[17] = T_s[17] + 273.15;

  T_ave = (T_s[1] + T_s[2] + T_s[3] + T_s[4] + T_s[11] + T_s[12] + T_s[13] + T_s[16] + T_s[5] + T_s[10])/10;
  T_ave2 = 0.07*T_s[1] + 0.175*T_s[2] + 0.175*T_s[3] + 0.07*T_s[5] + 0.07*T_s[7] + 0.05*T_s[9] + 0.19*T_s[11] + 0.20*T_s[13];
  E_v = 0.86 * sum(E) / 1.87;

  T0 = Tamb;


end Tanabe_65MN;
