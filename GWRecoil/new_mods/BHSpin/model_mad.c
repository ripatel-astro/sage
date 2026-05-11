#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>
#include <assert.h>

#include <stdbool.h>

#include "core_allvars.h"
#include "core_proto.h"


//!! added declaration before implementation
double calcSpinupNT(double spin);

// Forward declaration for RK45 integrator
double integrate_spin_RK45(double mass_msun, double initial_spin, double eddington_ratio,
                          double dt_years, double allowedSpinError, 
                          double minimumFractionalTimeResolution, int *substeps_taken);


double r_isco(double spin) {
    double Z1 = 1.0 + pow(1.0 - pow(spin, 2), 1.0 / 3.0) * (pow(1.0 + spin, 1.0 / 3.0) + pow(1.0 - spin, 1.0 / 3.0));
    double Z2 = sqrt(3 * pow(spin, 2) + pow(Z1, 2));
    double r = 3 + Z2 - copysign(sqrt((3.0 - Z1) * (3.0 + Z1 + 2 * Z2)), spin);
    return r;
}

double radiativeEfficiencyThin(double spin) {
    return 1.0 - sqrt(1.0 - 2.0 / 3.0 / r_isco(spin));
}

double calcEddingtonAccretionRate(double mass, double spin) {

    // these are in SI units
    const double G_SI = 6.67430e-11;              // m³/(kg·s²) 
    const double m_p_SI = 1.6726219e-27;          // kg
    const double c_SI = 299792458.0;              // m/s
    const double sigma_T_SI = 6.6524587e-29;      // m²
    const double Msun_SI = 1.98847e30;            // kg
    const double yr_SI = 31557600.0;              // s

    double efficiency = radiativeEfficiencyThin(spin);
    double mass_kg = mass * Msun_SI;

    // Returns Eddington accretion rate in solar masses per year
    return 4 * M_PI * G_SI * mass_kg * m_p_SI / (efficiency * sigma_T_SI * c_SI) * yr_SI / Msun_SI;
}

double calcMaximumMagnetization(double spin) {
    return -20.2 * pow(spin, 3) - 14.9 * pow(spin, 2) + 34 * spin + 52.6;
}

double spinToHorizon(double a) {
    return 1.0 + sqrt(1.0 - pow(a, 2));
}

double Omega_H(double spin) {
    return fabs(spin) / (2 * spinToHorizon(spin));
}

double sHD_min(double spin) {
    return 0.86 - 2 * spin * 0.97;
}


double sHD(double spin, double eddingtonRatio) {
    double parameters = 0.01971853;
    double xi = parameters * eddingtonRatio;
    return (calcSpinupNT(spin) + sHD_min(spin) * xi) / (1 + xi);
}

double calcMagnetization(double spin, double eddingtonRatio) {

    if (eddingtonRatio > fEdd_ADAF_threshold) {
        double evolutionFactor = pow(eddingtonRatio / 1.88, 1.29);
        return calcMaximumMagnetization(spin) * evolutionFactor / (1.0 + evolutionFactor);
    } else {
        return calcMaximumMagnetization(spin);
    }
}

double jetEfficiencyFunction(double phi, double spin) {
    double kappa = 0.05;
    double horizonAngularVelocity = Omega_H(spin);
    return kappa / (4 * M_PI) * pow(phi, 2) * pow(horizonAngularVelocity, 2) * (1.0 + 1.38 * pow(horizonAngularVelocity, 2) - 9.2 * pow(horizonAngularVelocity, 4));
}

double etaEM(double spin, double eddingtonRatio) {
    return jetEfficiencyFunction(calcMagnetization(spin, eddingtonRatio), spin);
}

double k_ratio(double spin) {
    if (spin <= 0) {
        return 0.23;
    } else {
        return fmin(0.35, 0.1 + 0.5 * spin);
    }
}

double calcJetPower(double mass, double spin, double eddingtonRatio) {

    // Jet power: P = eta_EM * f_Edd * M_dot_Edd * c^2
    // Returns power in Watts (SI)
    
    const double c_SI = 299792458.0;              // m/s
    const double Msun_SI = 1.98847e30;            // kg
    const double yr_SI = 31557600.0;              // s
    
    double mdot_edd = calcEddingtonAccretionRate(mass, spin);  // Msun/yr
    
    // Convert Msun/yr to kg/s, then multiply by c^2. Jet power in units of Watts
    return etaEM(spin, eddingtonRatio) * eddingtonRatio * mdot_edd 
           * Msun_SI / yr_SI * c_SI * c_SI;

}

double calcMassEvolutionRate(double mass, double spin, double eddingtonRatio) {
    return calcEddingtonAccretionRate(mass, spin) * eddingtonRatio;
}

double calcSpinupNT(double spin) {
    double a = fmax(-1, fmin(1, spin));
    double r = r_isco(spin);
    double E = (pow(r, 2) - 2 * r + a * sqrt(r)) / (r * sqrt(pow(r, 2) - 3 * r + 2 * a * sqrt(r)));
    double l = sqrt(r) * (pow(r, 2) - 2 * a * sqrt(r) + pow(a, 2)) / (r * sqrt(pow(r, 2) - 3 * r + 2 * a * sqrt(r)));
    return l - 2 * a * E;
}

double calcSpinupMAD(double spin) {
    double fitCoefficients[6] = {0.44763054, -12.53175763, -7.79639428, 9.43649497, 5.70845519, -4.02630245};
    double result = 0;
    for (int i = 0; i < 6; ++i) {
        result += fitCoefficients[i] * pow(spin, i);
    }

    return result;
}

double calcSpinup(double spin, double eddingtonRatio) {

    double s;
    if (eddingtonRatio < fEdd_ADAF_threshold) {
        // You're a MAD ADAF
        return calcSpinupMAD(spin);
    } else {
        // Start with hydrodynamic spinup.
        s = sHD(spin, eddingtonRatio);
        if (spin != 0) {
            // Then, add electromagnetic spinup.
            s += -etaEM(spin, eddingtonRatio) * (1.0 / k_ratio(spin) / Omega_H(spin) - 2.0 * spin) * copysign(1, spin);
        }
        return s;
    }
}



double calcSpinEvolutionRate(double mass, double spin, double eddingtonRatio) {
    return calcSpinup(spin, eddingtonRatio) * eddingtonRatio * calcEddingtonAccretionRate(mass, spin) / mass;
}

//!! added stadard coefficients definition .. it makes no sense to put them in struct at each instantiation of the struct
const struct {
    double _A[6];
    double _B[6][5]; 
    double _C[6] ;
    double _CH[6];
    double _CT[6];	
}RKFCoefficients = {
	{0.0, 2.0 / 9.0, 1.0 / 3.0, 3.0 / 4.0, 1.0, 5.0 / 6.0},//_A
	{{0.0, 0.0, 0.0, 0.0, 0.0}, 
     {2.0 / 9.0, 0.0, 0.0, 0.0, 0.0}, 
     {1.0 / 12.0, 1.0 / 4.0, 0.0, 0.0, 0.0}, 
	 {69.0 / 128.0, -243.0 / 128.0, 135.0 / 64.0, 0.0, 0.0}, 
	 {-17.0 / 12.0, 27.0 / 4.0, -27.0 / 5.0, 16.0 / 15.0, 0.0}, 
	 {65.0 / 432.0, -5.0 / 16.0, 13.0 / 16.0, 4.0 / 27.0, 5.0 / 144.0}},//_B
    {1.0 / 9.0, 0.0, 9.0 / 20.0, 16.0 / 45.0, 1.0 / 12.0},//_C
    {47.0 / 450.0, 0.0, 12.0 / 25.0, 32.0 / 225.0, 1.0 / 30.0, 6.0 / 25.0}, //_CH
	{-1.0 / 150.0, 0.0, 3.0 / 100.0, -16.0 / 75.0, -1.0 / 20.0, 6.0 / 25.0} //_CT	
};


double integrate_spin_RK45(double mass_msun, double initial_spin, double eddington_ratio,
                          double dt_years, double allowedSpinError, double minimumFractionalTimeResolution,
                          int *substeps_taken)  // ADD THIS - return number of steps taken
{
    double spin = initial_spin;
    double time = 0.0;
    double substep = dt_years * 0.01;  // Initial substep (1% of full timestep)
    
    const int MAX_SUBSTEPS = 100;  // Hard limit as requested
    int substep_count = 0;
    

    // Integrate until we reach dt_years
    while(time < dt_years) {
        
        if(substep_count >= MAX_SUBSTEPS) {
            assert(substep_count < MAX_SUBSTEPS);  
        }
        
        // Don't overshoot the target time
        if(time + substep > dt_years) {
            substep = dt_years - time;
        }      
        //

        int step_accepted = 0;
        double proposed_spin;
        
        while(!step_accepted) {
            // RK45 substep calculations
            double k_spin[6];
            
            for(int k = 0; k < 6; k++) {
                double spin_pred = spin;
                for(int l = 0; l < k; l++) {
                    spin_pred += k_spin[l] * RKFCoefficients._B[k][l];
                }
                
                // Enforce spin bounds during prediction
                if(spin_pred > 0.998) spin_pred = 0.998;
                if(spin_pred < -0.998) spin_pred = -0.998;
                
                // Calculate da/dt at this prediction point
                double da_dt = calcSpinEvolutionRate(mass_msun, spin_pred, eddington_ratio);
                k_spin[k] = substep * da_dt;
            }
            
            // Proposed new spin (5th order)
            proposed_spin = spin;
            for(int k = 0; k < 6; k++) {
                proposed_spin += RKFCoefficients._CH[k] * k_spin[k];
            }
            
            // Enforce bounds
            if(proposed_spin > 0.998) proposed_spin = 0.998;
            if(proposed_spin < -0.998) proposed_spin = -0.998;
            
            // Estimate truncation error
            double truncation_error = 0.0;
            for(int k = 0; k < 6; k++) {
                truncation_error += RKFCoefficients._CT[k] * k_spin[k];
            }
            truncation_error = fabs(truncation_error);
            
            // Calculate optimal step size factor
            double stepsize_factor = 0.9 * pow(allowedSpinError / (truncation_error + 1e-30), 0.2);
            
            // Accept step if error is within tolerance
            if(stepsize_factor >= 1.0 || substep < dt_years * minimumFractionalTimeResolution) {
                step_accepted = 1;
                spin = proposed_spin;
                time += substep;
                substep_count++;  // Count accepted steps
            }
            
            // Adjust substep size for next iteration
            substep *= fmin(stepsize_factor, 2.0);
            
            // Ensure minimum resolution
            if(minimumFractionalTimeResolution > 0 && substep < dt_years * minimumFractionalTimeResolution) {
                substep = dt_years * minimumFractionalTimeResolution;
            }
        }
    }

    // Return number of substeps taken
    if(substeps_taken != NULL) {
        *substeps_taken = substep_count;
    }
    
    return spin;
}



// Simple wrapper function to update spin during accretion
void update_black_hole_spin(int p, double accreted_mass, double dt) {

    //static int debug_count = 0;
    double mass, spin, eddington_ratio, mdot_edd, mdot_actual, da_dt, new_spin;
    double magnetization, jet_power, jet_efficiency, thindisk_RadiativeEfficiency;
    static int total_calls = 0;
    static int total_substeps = 0;


    if(SpinEvolutionOn == 0) return;  // Skip if spin evolution is off
    
    mass = Gal[p].BlackHoleMass;
    spin = Gal[p].BlackHoleSpin;
    
    // Skip if BH is too small
    if(mass < 1e-10) {
        Gal[p].BlackHoleSpin = SpinSeedValue;
        Gal[p].BH_Magnetization = 0.0; 
        Gal[p].BH_JetPower = 0.0;
        Gal[p].Eddrate = 0.0;
        return;
    }
    
    // Calculate Eddington ratio  
    mass = mass*1e10 / Hubble_h;                                                 // Convert to solar masses
    mdot_edd = calcEddingtonAccretionRate(mass, spin);  
    mdot_actual = (accreted_mass*1e10/Hubble_h) / ( dt * UnitTime_in_s / SEC_PER_YEAR);
    eddington_ratio = mdot_actual / mdot_edd;
    
    // Protect against extreme values
    if(eddington_ratio > 1000.0) eddington_ratio = 1000.0;
    if(eddington_ratio < 1e-10) return;  // No significant accretion

    // RK45 parameters
    double allowedSpinError = 1e-6;
    double minimumFractionalTimeResolution = 0.0;
    double dt_years = dt * UnitTime_in_s / SEC_PER_YEAR;
    
    // Track substeps
    int substeps_taken = 0;
    
    // Integrate spin using RK45
    new_spin = integrate_spin_RK45(mass, spin, eddington_ratio, dt_years, 
                                   allowedSpinError, minimumFractionalTimeResolution,
                                   &substeps_taken);
    
    // Tracking number of substeps taken just in case
    total_calls++;
    total_substeps += substeps_taken;
        
    Gal[p].BlackHoleSpin = new_spin;

    magnetization = calcMagnetization(Gal[p].BlackHoleSpin, eddington_ratio);
    jet_power = calcJetPower(mass, Gal[p].BlackHoleSpin, eddington_ratio);  
    jet_efficiency = etaEM(Gal[p].BlackHoleSpin, eddington_ratio);
    thindisk_RadiativeEfficiency = radiativeEfficiencyThin(Gal[p].BlackHoleSpin);

    Gal[p].Eddrate = eddington_ratio;
    Gal[p].BH_Magnetization = magnetization; 
    Gal[p].BH_JetPower = jet_power;                                             // in Watts
    Gal[p].BH_JetEfficiency = jet_efficiency;
    Gal[p].BH_thindiskRadiativeEfficiency = thindisk_RadiativeEfficiency;

}


double calculate_merger_spin_aligned(double M1, double a1, double M2, double a2)
{
    // Final BH spin from aligned merger
    double q, nu, a_final, l_orb;
    double s4, s5, t0, t2, t3;
    
    // Ensure M1 >= M2 (swap if needed)
    if(M2 > M1) {
        double temp_M = M1, temp_a = a1;
        M1 = M2; a1 = a2;
        M2 = temp_M; a2 = temp_a;
    }
    
    // Handle edge cases
    if(M1 < 1e-10) return SpinSeedValue;  // Both essentially zero
    if(M2 < 1e-10) return a1;              // Only M1 has mass
    
    // Mass ratio q = M2/M1 <= 1
    q = M2 / M1;
    
    // Symmetric mass ratio nu = q/(1+q)^2
    nu = q / ((1.0 + q) * (1.0 + q));
    
    // Fitting coefficients from Rezzolla et al. (2008) Table 1
    // These come from fitting to numerical relativity simulations
    s4 = -0.129;
    s5 = -0.384;
    t0 = -2.686;
    t2 = -3.454;
    t3 = 2.353;
    
    // Specific orbital angular momentum at merger
    // This captures the contribution from the orbital motion
    l_orb = (s4 / (1.0 + q*q)) * (a1*a1 + a2*a2*q*q*q*q + 2.0*a1*a2*q*q) 
            + (s5*nu + t0 + 2.0) / (1.0 + q*q) 
            + t2*nu + t3*nu*nu;
    
    // Final spin formula (Rezzolla et al. 2008, Eq. 7)
    a_final = (a1 + a2*q*q + l_orb*q) / ((1.0 + q) * (1.0 + q));
    
    // Enforce physical limits: -0.998 < a < 0.998
    // (Using 0.998 instead of 1.0 to match your other spin limits)
    if(a_final > 0.998) a_final = 0.998;
    if(a_final < -0.998) a_final = -0.998;
    
    return a_final;
}


