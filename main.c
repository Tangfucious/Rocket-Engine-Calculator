#include <stdio.h>
#include <math.h>

//#define Mt 1.0
#define g 9.80665
#define Pi 3.1415926
#define Mt 1
#define Pe 101325
#define a 340.29

//LOX alcohol(75%) data

int main()
{
    double Fv = 9800*100; //VACUUM THRUST - N 5, 250
    double nc = 0.98;
    double nn = 0.98;
    double Is = 338; //SPECIFIC IMPULSE 358(RP1) 338(ethanol)
    double y = 1.26;//RP1
    double Qmc = 0; //MASS FLOW
    double R = 201.914;//GAS CONSTANT OF 518.28(CH4) 379.6(RP1)
    double Tc = 3230;//2773.15;//;//3083.15
    double MR = 1.43;//mix ratio of fuel and oxidizer. 2.327(RP1)
    double Pc = 5e+6; //UNIT IN Pa, THE PRESSURE OF THE COMBUSTION CHMABER
    double Lc = 1; //CHARACTERISTIC LENGTH OF THE

    // COMBUSTION CHAMBER
    double ec; //COMPRESSION RATIO OF THE ENGINE
    double ee; //EXPANSION RATIO OF THE ENGINE
    double p = 2.5, k = 1.5;

    //fuel and oxidizer
    double Fuel_Weight_PS;//mass of fuel runs per second
    double Oxi_Weight_PS;//mass of oxidizer runs per second

    //wall thickness
    double tw;//wall thick
    double S = 5.5158e+7;// work stress of copper

    //propellant tank and gas tank
    double Pp_oxi = 5e+6;//pressure of oxidizer tank
    double Pp_fuel = 5e+6;//pressure of fuel tank
    double Vp_oxi = 10*0.001;//volume of oxidizer tank
    double Vp_fuel = 5*0.001;//volume of fuel tank
    double P0_oxi;//initial gas tank pressure for oxidizer
    double P0_fuel;//initial gas tank pressure for fuel
    double V0_oxi = 5*0.001;//volume of gas tank for oxidizer, convert liters to m^3
    double V0_fuel = 5*0.001;//volume of gas tank for fuel, convert liters to m^3
    double Pg_oxi = 3.6e6;//final gas tank pressure for oxidizer
    double Pg_fuel = 3.6e6;//final gas tank pressure for fuel
    double m0_oxi;
    double m0_fuel;
    double R_gas = 296.7;
    double T0_gas = 63.17;

    //variables solve for
    double At;//AREA OF THE THROAT
    double Ae;//AREA AT EXIT
    double Ac;//COMBUSTION CHAMBER AREA
    double Dt;//THROAT DIAMETER
    double De;//EXIT DIAMETER
    double Dc;//COMBUSTION CHAMBER DIAMETER
    double R2;
    double R1;
    double R2_volume;
    double R1_volume;
    double Me;//MACH NUMBER AT EXIT
    double Rn = 0;
    double L_star =1;//characteristic length
    double Lc1;
    double Lc2;
    double h;
    double H;
    double Y;
    double ln = 6.5135;
    double Ln;
    double beta = 25;//exit angle
    beta = beta*(Pi/180);//convert to radian
    double R0;
    double ve;


    //Pc = Pe/pow(((pow(Me, 2)*(y+1))/2)+1, -(-(y-1)/y));
    //ve = sqrt(((2*y)/(y-1))*R*Tc*(1-pow(Pe/Pc, (y-1)/y)));
    //Me = ve/a;
    Me = sqrt((2 * (pow((Pe / Pc), -(y - 1) / y) - 1)) / (y - 1));
    ve = Me*a;
    Qmc = Fv/(ve*nc*nn);
    //Qmc = Fv/(Is*g);
    //solving for the throat
    At = (Qmc*pow(R*Tc, 1.0/2)) / (Pc*pow(y * pow(2 / (y + 1), (y + 1) / (y - 1)), 1.0/2));
    Dt = sqrt((4 * At) / Pi);
    //Vc = (R * Tc) / Pc;

    //solving for the combustion chamber;
    ec = 8 * pow(Dt*100, -0.6) + 1.25;
    Ac = At * ec;
    Dc = sqrt((4 * Ac) / Pi);
    //solving for the exit area
    Ae = At * (1 / Me) * pow((2 / (y + 1)) * (1 + ((y - 1) / 2) * pow(Me, 2)), (y + 1) / (2 * (y - 1)));
    De = sqrt((4 * Ae) / Pi);
    ee = Ae / At;
    Lc2 = (Dt / 2) * sqrt(pow((k + p * sqrt(ec)), 2) - pow(((p - 1) * sqrt(ec) + 1 + k), 2));
    h = Lc2 * (k / (k + p * sqrt(ec)));//throat
    H = Lc2 - h;//conversion
    Y = k * (Dt / 2) + (Dt / 2) - sqrt(pow(k, 2) * pow((Dt / 2), 2) - pow(h, 2));
    R2 =(pow(H, 2)+pow(((Dc/2)-Y), 2))/(2*((Dc/2)-Y));
    R1 = (pow(h, 2)+pow((Y-(Dt/2)), 2))/(2*((Y-(Dt/2))));

    double Rc = Dc/2;
    double Rt = Dt/2;

    R2_volume = Pi * (-asin(H/fabs(R2)) * pow(R2, 2) * (R2-Rc) - H * sqrt(pow(R2, 2)-pow(H, 2)) * (R2-Rc) - pow(H, 3)/3 + H * (2*pow(R2, 2)-2*Rc*R2+pow(Rc, 2)));
    R1_volume = Pi * (-asin(h/fabs(R1)) * pow(R1, 2) * (R1+Rt) - h * sqrt(pow(R1, 2)-pow(h, 2)) * (R1+Rt) - pow(h, 3)/3 + h * (2*pow(R1, 2)+2*Rt*R1+pow(Rt, 2)));
    double nozzle_volume = R1_volume+R2_volume;
    printf("R2 Volume: %lf\n", R2_volume);
    printf("R2 Volume: %lf\n", R1_volume);
    printf("Volume: %lf\n", nozzle_volume);

    //double nozzle_volume = 0.00042725272;//m^3, from solidworks. 1 ton: 0.00051712985
    //double nozzle_volume = R1_volume+R2_volume;
    Lc1 = (L_star*At-nozzle_volume)/(Pi*pow(Dc/2, 2));
    R0 = ((pow(ln, 2) + pow((1.5 - (De / Dt)), 2) - 1) / (2 * (1 - ln * (sin(beta* Pi/180)) - ((1.5 - (De / (2*Dt))) * (cos(beta* Pi/180))))))*Dt;
    Ln = ln * Dt;
    //solving oxidizer mass flow rate and fuel mass flow rate
    Fuel_Weight_PS = Qmc/(1+MR);
    Oxi_Weight_PS = Fuel_Weight_PS*MR;
    //wall thickness
    tw = (Pc*Dc)/(2*S);

    //injector
    //suppose cd is one for all of them
    //we have eight oxidizer injector and 1 fuel injector
    //oxidizer injector
    double Cd = 0.6;//0.735
    double rough_oxi = 1141;
    double rough_fuel = 789;
    double Pp = 1.5e+6;//2.5 Mpa, converted to Pa
    double A_oxi = (Oxi_Weight_PS/4)/(Cd*sqrt(2*rough_oxi*(Pp-Pc)));
    double D_oxi = 2*sqrt(A_oxi/Pi);
    double A_fuel = (Fuel_Weight_PS)/(Cd*sqrt(2*rough_fuel*(Pp-Pc)));
    double D_fuel = 2*sqrt(A_fuel/Pi);

    //volume of propellant tanks
    //fuel tank
    //gas tank for oxidizer tank
    P0_oxi = ((Pp_oxi*Vp_oxi)/V0_oxi)+Pg_oxi;
    m0_oxi = (Pp_oxi*Vp_oxi/(R_gas*T0_gas))*(1/(1-(Pg_oxi/P0_oxi)));

    //has tank for fuel tank
    P0_fuel = ((Pp_fuel*Vp_fuel)/V0_fuel)+Pg_fuel;
    m0_fuel = (Pp_fuel*Vp_fuel/(R_gas*T0_gas))*(1/(1-(Pg_fuel/P0_fuel)));

    //regenerative cooling
    //wall temperature should be under the melting point, take 800 degree celcius
    double Mole = 37.767;//mole of the compound
    double c_star;
    double area[60];//cross-sectional area of each section
    double temperature[60];
    double mach[60];
    double T_cool[60];
    double ug[60];
    double Taw[60];//adiabatic all temperature
    double alpha[60];//the correction factor
    double Tcool[60];//temperature of the coolant
    double viscosity[60];
    double v_H2O[60];//viscosity of H20
    double v_CO2[60];//viscosity of CO2
    double C_H2O = 1064;//Sutherland's constant of water
    double C_CO2 = 222;//Sutherland's constant of CO2
    double v0_H2O = 1.12*pow(10, -5);//base viscosity of H2O
    double v0_CO2 = 1.370*pow(10, -5);//base viscosity of CO2
    double T0_H2O = 350;//base temperature of H2O
    double T0_CO2 = 273;//based temperature of CO2
    double H2O_mfrac = 0.6;//mole fraction
    double CO2_mfrac = 0.4;//mole fraction
    double wall_temp = 1073.15;//800 degree celcius
    double Pr;
    double T;
    double Tt;
    double Te;
    double M;
    double Tw;//inner wall temperature
    double convec_TS;//convection temperature step
    double div_TS;//expansion

    c_star = Pc*At/Qmc;
    Pr = y/(y+(10.4/Mole));
    Tt = Tc*pow(1+(y-1)/2, -1);
    Te = Tt/(1+(y-1)/2*pow(Me,2));
    convec_TS = (Tc-Tt)/20;
    div_TS = (Tt-Te)/20;
    //speed in combustion chamber
    T = Tc-convec_TS;
    M = sqrt(((Tc/T)-1)/((y-1)/2));
    for(int i = 0; i<20; i++)
    {
        temperature[i] = Tc;
        mach[i] = M;
    }
    T = Tc;
    for(int i = 20; i<40; i++){
        T = T-convec_TS;
        M = sqrt(((Tc/T)-1)/((y-1)/2));
        temperature[i] = T;
        mach[i] = M;
    }
    T = Tt;
    for(int i = 40; i<60; i++){
        T = T-div_TS;
        M = sqrt(((Tc/T)-1)/((y-1)/2));
        temperature[i] = T;
        mach[i] = M;
    }
    for(int i = 0; i<60; i++)
    {
        //printf("T:%lf\n", temperature[i]);
        //printf("M:%lf\n", mach[i]);
    }
    for(int i = 0; i<60; i++)
    {
        Taw[i] = (Tc*(1+pow(Pr, 1.0/3)*((y-1)/2)*pow(mach[i], 2)))/(1+((y-1)/2)*pow(mach[i], 2));
        //printf("M: %lf\n", mach[i]);
        //printf("Taw: %lf\n", Taw[i]);
    }
    for(int i = 0; i<60; i++)
    {
        alpha[i] = pow((0.5*(wall_temp/Tc)*(1+((y-1)/2)*pow(mach[i], 2))+0.5), -0.86)*
                    pow(1+((y-1)/2)*pow(mach[i], 2), -0.12);
        //printf("alpah: %lf\n", alpha[i]);
    }
    for(int i = 0; i<60; i++)
    {
        v_H2O[i] = (v0_H2O*pow(temperature[i]/T0_H2O, 3.0/2)*(T0_H2O+C_H2O))/(C_H2O+temperature[i]);
        v_CO2[i] = (v0_CO2*pow(temperature[i]/T0_CO2, 3.0/2)*(T0_CO2+C_CO2))/(C_CO2+temperature[i]);
//        printf("v_H2O: %lf\n", v_H2O[i]);
//        printf("v_CO2: %lf\n", v_CO2[i]);
    }
    for(int i = 0; i<60; i++)
    {
        viscosity[i] = v_H2O[i]*H2O_mfrac+ v_CO2[i]*CO2_mfrac;
        //printf("viscosity: %lf\n", viscosity[i]);
    }
    for(int i = 0; i<60; i++)
    {
        ug[i] = (0.026/(pow(Dt, 0.2)))*((pow(viscosity[i], 0.2)*y)/(pow(Pr, 0.6)))
                *pow(Pc/c_star, 0.8)*(At/area[i], 0.8)*alpha[i];
    }
//    for(int i = 0; i<60; i++){
//        q[i] = ug[i]+(Taw[i]-Tw[i])
//    }
//    for(int i = 1; i<60; i++){
//        Tc =
//    }


    //cooling system
    double mc;//mass flow rate of film
    double P = 2*Pi*(Dc/2);//diameter of combustion chamber
    double film_length = Lc1;//this is the total length of the rocket engine, 0.54291
    double film_Qc = 213;//film latent heat evaporation
    double film_n = 0.5;//cooling - efficiency of the film
    double film_Tt = 54.36;//initial temperature
    double film_Ts = 90.19;//boiling temperature
    double film_hg;//convection constant of film
    double film_Tg = Tc;
    double film_cp = 0.918;//J/K
    double G = (Pi*pow((Dc/2), 2));//mass flow rate/area
    double L = Lc1;//total length of the chamber
    //film cooling
    film_hg = ((3.075*film_cp*pow(G,0.8))/pow(Dc, 0.2))*(1+pow((Dc/L), 0.7));
    mc = ((film_length/film_n)*P*film_hg*(film_Tg-film_Ts))/(film_cp*(film_Ts-film_Tt)+film_Qc);
    double film_percentage = (mc/Qmc)*100;
    double A_film = (mc/26)/(Cd*sqrt(2*rough_oxi*(Pp-Pc)));
    double D_film = 2*sqrt(A_film/Pi);
    double injector_height = 2*sqrt(((Oxi_Weight_PS)/(Cd*sqrt(2*rough_oxi*(Pp-Pc))))/Pi)
            +2*sqrt(((mc)/(Cd*sqrt(2*rough_oxi*(Pp-Pc))))/Pi);

    printf("Qmc: %lf\n", Qmc);
    printf("Pc: %lf\n", Pc);
    printf("R: %lf\n", R);
    printf("ve: %lf\n", ve);
    printf("Me: %lf\n", Me);
    printf("At: %lf\n", At);
    printf("Dt: %lf\n", Dt);
    printf("Ae: %lf\n", Ae);
    printf("De: %lf\n", De);
    printf("Dc: %lf\n", Dc);
    printf("R2: %lf\n", R2);
    printf("R1: %lf\n", R1);
    printf("Rn: %lf\n", Rn);
    printf("Lc: %lf\n", Lc);
    printf("Lc1: %lf\n", Lc1);
    printf("Lc2: %lf\n", Lc2);
    printf("H: %lf\n", H);
    printf("h: %lf\n", h);
    printf("Y: %lf\n", Y);
    printf("Lc: %lf\n", Lc);
    printf("R0: %lf\n", R0);
    printf("Ln: %lf\n", Ln);
    printf("ec: %lf\n", ec);
    printf("ee: %lf\n", ee);
    printf("R2 Volume: %lf\n", R2_volume);
    printf("R1 Volume: %lf\n", R1_volume);
    printf("Nozzle Volume: %lf\n", nozzle_volume);
    printf("Fuel Per Second: %lf\n", Fuel_Weight_PS);
    printf("Oxidizer Per Second: %lf\n", Oxi_Weight_PS);
    printf("Pressure Oxidizer: %lf\n", P0_oxi);
    printf("Pressure Fuel: %lf\n", P0_fuel);
    printf("Mass Oxidizer Pressurizer: %lf\n", m0_oxi);
    printf("Mass Fuel Pressurizer: %lf\n", m0_fuel);
    printf("Film FLow Amount: %lf\n", mc);
    printf("Film Percentage: %lf\n", film_percentage);
    printf("hg: %lf\n", film_hg);
    printf("tw: %lf\n", tw);
    printf("Area Oxidizer: %lf\n", A_oxi);
    printf("Diameter Oxidizer: %lf\n", D_oxi);
    printf("Area Fuel: %lf\n", A_fuel);
    printf("Diameter Fuel: %lf\n", D_fuel);
    printf("Injector Height: %lf\n", injector_height);
    printf("Film Injector Diameter: %lf\n", D_film);
    return 0;
}
