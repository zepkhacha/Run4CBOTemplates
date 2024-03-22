void plotAll(){

    plotParam("slidingFits4/noRF_windowFits.root", "alpha_CBO", "dimensionless", -0.01, 0.01, "ALE PFC PLC");
    plotParam("slidingFits4/noRF_windowFits.root", "beta_CBO", "dimensionless", -0.01, 0.01, "ALE PFC PLC");
    plotParam("slidingFits4/noRF_windowFits.root", "A_CBO", "dimensionless", 0.0, 0.012, "ALE PFC PLC");
    plotParam("slidingFits4/noRF_windowFits.root", "phi_CBO", "[rad]", -3.2, 3.2, "ALE PFC PLC");
    plotParam("slidingFits4/noRF_windowFits.root", "w_CBO", "[rad/s]", 2.32, 2.34, "ALE PFC PLC");

    plotParam("slidingFits4/noRF_windowFits.root", "chisq", "[arb]", 0, 10000, "ALE PFC PLC");
    plotParam("slidingFits4/noRF_windowFits.root", "rchisq", "[arb]", 0, 50, "ALE PFC PLC");
    plotParam("slidingFits4/noRF_windowFits.root", "LM", "[arb]", -0.001, +0.001, "ALE PFC PLC");

    plotParam("slidingFits4/noRF_windowFits.root", "alpha_2CBO", "dimensionless", -0.002, 0.002, "ALE PFC PLC");
    plotParam("slidingFits4/noRF_windowFits.root", "beta_2CBO", "dimensionless", -0.002, 0.002, "ALE PFC PLC");
    plotParam("slidingFits4/noRF_windowFits.root", "A_2CBO", "dimensionless", 0.0, 0.002, "ALE PFC PLC");
    plotParam("slidingFits4/noRF_windowFits.root", "phi_2CBO", "[rad]", -3.2, 3.2, "ALE PFC PLC");

    plotParam("slidingFits4/noRF_windowFits.root", "alpha_y", "dimensionless", -0.002, 0.002, "ALE PFC PLC");
    plotParam("slidingFits4/noRF_windowFits.root", "beta_y", "dimensionless", -0.002, 0.002, "ALE PFC PLC");
    plotParam("slidingFits4/noRF_windowFits.root", "A_y", "dimensionless", 0.0, 0.002, "ALE PFC PLC");
    plotParam("slidingFits4/noRF_windowFits.root", "phi_y", "[rad]", -3.2, 3.2, "ALE PFC PLC");

    plotParam("slidingFits4/noRF_windowFits.root", "alpha_vw", "dimensionless", -0.002, 0.002, "ALE PFC PLC");
    plotParam("slidingFits4/noRF_windowFits.root", "beta_vw", "dimensionless", -0.002, 0.002, "ALE PFC PLC");
    plotParam("slidingFits4/noRF_windowFits.root", "A_vw", "dimensionless", 0.0, 0.002, "ALE PFC PLC");
    plotParam("slidingFits4/noRF_windowFits.root", "phi_vw", "[rad]", -3.2, 3.2, "ALE PFC PLC");

    plotParam("slidingFits4/noRF_windowFits.root", "alpha_xy", "dimensionless", -0.002, 0.002, "ALE PFC PLC");
    plotParam("slidingFits4/noRF_windowFits.root", "beta_xy", "dimensionless", -0.002, 0.002, "ALE PFC PLC");

    plotParam("slidingFits4/noRF_windowFits.root", "alpha_A0", "dimensionless", -0.02, 0.02, "ALE PFC PLC");
    plotParam("slidingFits4/noRF_windowFits.root", "beta_A0", "dimensionless", -0.02, 0.02, "ALE PFC PLC");
    plotParam("slidingFits4/noRF_windowFits.root", "A_A0", "dimensionless", 0.0, 0.02, "ALE PFC PLC");
    plotParam("slidingFits4/noRF_windowFits.root", "phi_A0", "[rad]", -3.2, 3.2, "ALE PFC PLC");
    
    plotParam("slidingFits4/noRF_windowFits.root", "alpha_phi", "dimensionless", -0.002, 0.002, "ALE PFC PLC");
    plotParam("slidingFits4/noRF_windowFits.root", "beta_phi", "dimensionless", -0.002, 0.002, "ALE PFC PLC");
    plotParam("slidingFits4/noRF_windowFits.root", "A_phi", "dimensionless", 0.0, 0.002, "ALE PFC PLC");
    plotParam("slidingFits4/noRF_windowFits.root", "phi_phi", "[rad]", -3.2, 3.2, "ALE PFC PLC");

}
