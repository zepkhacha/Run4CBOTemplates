#include <cmath>
#include <iostream>
#include "gm2util/blinders/Blinders.hh"

int desiredSeed = 0; // set via command-line argument
int desiredCalo = 0; // set via command-line argument

//histogram binning
int          startBin =   203;
const int    nBins    =  4357;
double       n        = 0.108;
double       wc_0     = 0.0;

// histogram declarations
static TH1D* wiggle = new TH1D("wiggle", "wiggle", nBins, 0.0, 650.0644);
static TH1D* lambda = new TH1D("lambda", "lambda", nBins, 0.0, 650.0644);
static TH1D* varcorr = new TH1D("varcorr", "varcorr", nBins, 0.0, 650.0644);
static TH1D* wiggleraw = new TH1D("wiggleraw", "wiggleraw", nBins, 0.0, 650.0644);

//fit type switches
bool pBinTau        = false; // momentum-binning switches & constraint
bool pBinPhi        = false; 
bool pBinCBO        = false;
bool includeMopTerm = false; 
bool constrainMop   = false; // mop constraint via chisq penalty

//path to full CBO templates
std::string templatePath;
std::string refFile; // to fix w_CBO ? not implemented yet
// to store alpha and beta TF1 templates
std::vector<TF1> alpha_CBO_TF1;
std::vector<TF1> beta_CBO_TF1;
std::vector<TF1> alpha_2CBO_TF1;
std::vector<TF1> beta_2CBO_TF1;
std::vector<TF1> alpha_y_TF1;
std::vector<TF1> beta_y_TF1;
std::vector<TF1> alpha_vw_TF1;
std::vector<TF1> beta_vw_TF1;
std::vector<TF1> alpha_A0_TF1;
std::vector<TF1> beta_A0_TF1;
std::vector<TF1> alpha_phi_TF1;
std::vector<TF1> beta_phi_TF1;

//bin-dependent variables
std::pair <double, double> magicTauConstraints (0.0,0.0);
std::pair <int,int> binRange (0,0); 
int magicBin = -10;
std::vector<double> binWeights;
std::vector<double> dp_p0;
std::vector<double> gammas;
std::vector<double> dPhi;
double weightedAvgTau;

//weights for calo-combination CBO term
std::vector<double> caloWeights(24, 0.0); // vector of weights, 0-based!!

// e-field 
double c_e_scale = 1.0; // set from 0 - 1
std::vector<double> c_e; //in ppm
std::vector<double> x_e;
double R_0;
double beta_0;

// blinding for fitting
static blinding::Blinders myblinders(blinding::Blinders::kOmega_a, "reunificationisblind");

// number of free parameters
int noFreeParams = 0;

double wcbo_i(double dp){
    double wc = (1. - (dp / (1.-n))) * wc_0;
    double nu_x = sqrt((1.-n) - ((n*(2.+n)) / (1.-n))*dp);
    return (1. - nu_x) * wc;
}

double cbo(double time, double amp){
    // For the endgame
    //return 1.0 + 2.927*exp(-time/79.05)/time + 2.239*exp(-time/6.94)/time;
    // For Run 2/3
    //return 1.0 + amp*exp(-time/24.4)/time;
    return 1.0;
}

double wy(double kappa, double wcbo, double time, double amp){
    double x = ((4*M_PI)/(0.1492*kappa*wcbo*cbo(time, amp))) - 1.0;
    return kappa*wcbo*cbo(time, amp)*sqrt(x);
}

double wvw(double kappa, double wcbo, double time, double amp){
    return ((2*M_PI)/0.1492) - 2.*wy(kappa, wcbo, time, amp);
}

double calcnu(double *dim, double *par){ // dim[0] = bin number
    double time = (dim[0]-0.5) * 0.1492 ;
    double blindr = myblinders.paramToFreq(par[4]);
    double phi = par[3];
    // common in all bins
    double nu = 0.0;

    int calcnuRangeLower = magicBin-1;
    int calcnuRangeUpper = magicBin+1;
    if (pBinTau or pBinPhi){
        calcnuRangeLower = binRange.first;
        calcnuRangeUpper = binRange.second;
    }

    // start iterating bin-by-bin
    // in kevin-style fits, b==magicBin has weight = 1.0, remainder have 0
    for (int b=calcnuRangeLower; b<=calcnuRangeUpper; b++){

        // first let's calculate all of our weighted-avg alpha and beta values
        double alpha_CBO = 0.0;
        double beta_CBO  = 0.0;
        double alpha_2CBO= 0.0;
        double beta_2CBO = 0.0;
        double alpha_y   = 0.0;
        double beta_y    = 0.0;
        double alpha_vw  = 0.0; 
        double beta_vw   = 0.0; 
        double alpha_A0  = 0.0;
        double beta_A0   = 0.0;
        double alpha_phi = 0.0;
        double beta_phi  = 0.0;

        for (int caloNum=0; caloNum<24; caloNum++){
            if (desiredCalo!=0 and caloNum!=desiredCalo-1){
                continue;
            }
            alpha_CBO  += caloWeights[caloNum]*(alpha_CBO_TF1 [caloNum]).Eval(time);
            beta_CBO   += caloWeights[caloNum]*(beta_CBO_TF1  [caloNum]).Eval(time);
            alpha_2CBO += caloWeights[caloNum]*(alpha_2CBO_TF1[caloNum]).Eval(time);
            beta_2CBO  += caloWeights[caloNum]*(beta_2CBO_TF1 [caloNum]).Eval(time);
            alpha_y    += caloWeights[caloNum]*(alpha_y_TF1   [caloNum]).Eval(time);
            beta_y     += caloWeights[caloNum]*(beta_y_TF1    [caloNum]).Eval(time);
            alpha_vw   += caloWeights[caloNum]*(alpha_vw_TF1  [caloNum]).Eval(time);
            beta_vw    += caloWeights[caloNum]*(beta_vw_TF1   [caloNum]).Eval(time);
            alpha_A0   += caloWeights[caloNum]*(alpha_A0_TF1  [caloNum]).Eval(time);
            beta_A0    += caloWeights[caloNum]*(beta_A0_TF1   [caloNum]).Eval(time);
            alpha_phi  += caloWeights[caloNum]*(alpha_phi_TF1 [caloNum]).Eval(time);
            beta_phi   += caloWeights[caloNum]*(beta_phi_TF1  [caloNum]).Eval(time);
        } // end loop over caloNum

        double b_dp = dp_p0[b];
        double binNu = 0.0;
        double phi_mod, A0_mod;

        // modulations on phi and asymmetry
        if (pBinCBO){
            // THIS ISNT CONFIGURED CORRECTLY FOR TEMPLATE FITS!
            //phi_mod = exp(-1.*(time)/par[5])*(par[20]*cos(par[6]*(wcbo_i(b_dp)/wcbo_i(0.0))*cbo(time, par[24])*time) 
            //        + par[21]*sin(par[6]*(wcbo_i(b_dp)/wcbo_i(0.0))*cbo(time, par[24])*time));
            //A0_mod      = exp(-1.*(time)/par[5])*(par[12]*cos(par[6]*(wcbo_i(b_dp)/wcbo_i(0.0))*time) 
            //        + par[13]*sin(par[6]*(wcbo_i(b_dp)/wcbo_i(0.0))*time));
        }else{
            phi_mod = par[20]*alpha_phi*cos(par[6]*cbo(time, par[24])*time) 
                + par[21]*beta_phi*sin(par[6]*cbo(time, par[24])*time);
            A0_mod  = par[12]*alpha_A0*cos(par[6]*cbo(time, par[24])*time) 
                + par[13]*beta_A0*sin(par[6]*cbo(time, par[24])*time);
        }   

        // 5-param term
        binNu    += pBinPhi ?    par[2]*cos((blindr*(1-c_e[b]))*time - ((phi + dPhi[b]) + phi_mod ))
            :                    par[2]*cos((blindr           )*time - ( phi            + phi_mod ));
        binNu    *= 1.0 + A0_mod;
        binNu    += 1.0;
        binNu    *= pBinTau ?    par[0]*exp(-1.*(time)/((gammas[b]/gammas[magicBin])*par[1]))
            : par[0]*exp(-1.*(time)/par[1]);

        // X & Y CBO
        double x, y;

        if (pBinCBO){
            // for now, there is NO calo-combination option for momentum-binned fits
            x = 1.0 +(exp(-1.*time/par[5])+ par[28])*(par[7] *cos(  par[6]*(wcbo_i(b_dp)/wcbo_i(0.0))*time) 
                    + par[8] *sin(  par[6]*cbo(time, par[24])*(wcbo_i(b_dp)/wcbo_i(0.0))*time))
                + exp(-2.*time/par[5]          )*(par[10]*cos(2*par[6]*(wcbo_i(b_dp)/wcbo_i(0.0))*time) 
                        + par[11]*sin(2*par[6]*cbo(time, par[24])*(wcbo_i(b_dp)/wcbo_i(0.0))*time));

            y = 1.0 + exp(-1.*time/par[17])*(par[14]*cos(wy (par[16], par[6]*(wcbo_i(b_dp)/wcbo_i(0.0)), time, par[24])*time) 
                    + par[15]*sin(wy (par[16], par[6]*(wcbo_i(b_dp)/wcbo_i(0.0)), time, par[24])*time))
                + exp(-2.*time/par[17])*(par[18]*cos(wvw(par[16], par[6]*(wcbo_i(b_dp)/wcbo_i(0.0)), time, par[24])*time) 
                        + par[19]*sin(wvw(par[16], par[6]*(wcbo_i(b_dp)/wcbo_i(0.0)), time, par[24])*time));

        } else{


            x = 1 + (par[7] *alpha_CBO *cos(  par[6]*cbo(time, par[24])*time)
                    +  par[8] * beta_CBO *sin(  par[6]*cbo(time, par[24])*time))
                + (par[10]*alpha_2CBO*cos(2*par[6]*cbo(time, par[24])*time) 
                        +  par[11]* beta_2CBO*sin(2*par[6]*cbo(time, par[24])*time));

            y = 1 + (par[14]*alpha_y *cos(wy (par[16], par[6], time, par[24])*time) 
                    +  par[15]*beta_y  *sin(wy (par[16], par[6], time, par[24])*time))
                + (par[18]*alpha_vw*cos(wvw(par[16], par[6], time, par[24])*time)
                        +  par[19]*beta_vw *sin(wvw(par[16], par[6], time, par[24])*time));


        }                                            
        double p = exp(-1.*time/par[25])*(par[22]*cos(time*par[26])+par[23]*sin(time*par[26]));
        binNu *= x*y + p;


        binNu *= (1.0 - par[9]*lambda->GetBinContent(int(dim[0])));
        nu += binWeights[b] * binNu;

    } // end loop over bins

    // if no mop term, p27=0; // if constrained mop, p27 ~ 0; // if floating mop, p27 ~ 2e-5
    nu *= exp(-time*par[27]);
    return nu;

} // end calcnu
