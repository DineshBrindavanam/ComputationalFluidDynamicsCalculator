#include <iostream>
#include<stdio.h>
#include<math.h>
#include<string.h>

using namespace std;

class fluid{
    public:
        fluid();
        double P_s, rho_s, T, u_avg, mu, T_s, R;
};
    fluid::fluid(){
    P_s = 101325, rho_s = 1.226, T = 300, u_avg = 1, mu = 1.785e-5, T_s = 300, R = 287;
    }

class geometry{
    public:
        geometry();
        double d, l, D_h, a, b;
};
    geometry::geometry(){
        d = 1, l = 1, D_h = 1, a = 1, b = 1;
    }

class mesh{
    public:
        mesh();
        double init_ht, er, nl, tot_ht, lay_ht, lst_cell_ht;
};
    mesh::mesh(){
        init_ht = 0, er = 1.1, nl =10, tot_ht = 0, lay_ht = 0, lst_cell_ht = 0;
    }

int i, t_no, nl;
double rho_s, P_s, R, T_s, ReNum, u, c_f, tau, mu, Y_Plus, y, x, er, init_ht, tot_ht, lay_ht, Ti, c_mu=0.09, Tls, mtv, tvr, tvr_KandL, tvr_KandMU, ke, e, e_KandL, e_KandMU, tv, tv_KandL, tv_KandMU, om;

double RenoldsNumber_Internal(){
    fluid air;
    geometry geo;
    cout << "\nEnter the static pressure of the fluid (N/m^2) (Default = 101325): ";
    cin >> air.P_s;
    cout << "\nEnter the static temperature of the fluid (Kelvin) (Default = 300): ";
    cin >> air.T_s;
    cout << "\nEnter the average velocity of the fluid (m/s) (Default = 1): ";
    cin >> air.u_avg;
    cout << "\nEnter the characteristic diameter (m) (Default = 1): ";
    cin >> geo.d;
    cout << "\nEnter the dynamic viscosity of the fluid (kg/m-s) (Default = 1.875e-5): ";
    cin >> air.mu;
    ReNum = ((air.P_s / (air.R * air.T_s)) * air.u_avg * geo.d) / air.mu;
    return ReNum;
}

double RenoldsNumber_External(){
    fluid air;
    geometry geo;
    cout << "\nEnter the static pressure of the fluid (N/m^2) (Default = 101325): ";
    cin >> air.P_s;
    cout << "\nEnter the static temperature of the fluid (Kelvin) (Default = 300): ";
    cin >> air.T_s;
    cout << "\nEnter the average velocity of the fluid (m/s) (Default = 1): ";
    cin >> air.u_avg;
    cout << "\nEnter the characteristic length (m) (Default = 1): ";
    cin >> geo.l;
    cout << "\nEnter the dynamic viscosity of the fluid (kg/m-s) (Default = 1.875e-5): ";
    cin >> air.mu;
    ReNum = ((air.P_s / (air.R * air.T_s)) * air.u_avg * geo.l) / air.mu;
    return ReNum;
}

double HydraulicDiameter_Circular(){
geometry geo;
cout << "\nEnter the Diameter of the cross section (m): ";
cin >> geo.D_h;
return geo.D_h/4;
}

double HydraulicDiameter_Square(){
geometry geo;
cout << "\nEnter the Length of the side of the cross section (m): ";
cin >> geo.a;
geo.D_h = geo.a;
return geo.D_h;
}

double HydraulicDiameter_Rectangular(){
geometry geo;
cout << "\nEnter the Length of the one side of the cross section (m): ";
cin >> geo.a;
cout << "\nEnter the Length of the other side of the cross section (m): ";
cin >> geo.b;
geo.D_h = (2 * geo.a * geo.b) / (geo.a +geo.b);
return geo.D_h;
}

double Y_plusCalc_Internal(){
fluid air;
geometry geo;
ReNum = RenoldsNumber_Internal();
c_f = 0.079 * pow(ReNum,-0.25);
tau = c_f * (0.5 * rho_s * pow(u,2));
double u_star = pow(tau / rho_s, 1/2);
cout << "\nRecommended values of Y+ are as follows:\n|\tTurbulence model\t|";
cout << "|\tS-A model\t|\t1 to 5\t\n|";
cout << "|\tk-e (standard wall function) model\t|\t1 to 5\t|\n";
cout << "|\tk-e (enhanced wall function) model\t|\t1 to 5\t|\n";
cout << "|\tk-w (low-re model) model\t|\t1 to 5\t|\n|";
cout << "|\tS-A model\t|\t1 to 5\t|\n";
cout << "Enter the required Y+ value: ";
cin >> Y_Plus;
y = (Y_Plus * air.mu) / (air.rho_s * u_star);
return y;
}

double Y_plusCalc_External(){
fluid air;
geometry geo;
ReNum = RenoldsNumber_External();
c_f = 0.058 * pow(ReNum,-0.2);
tau = c_f * (0.5 * rho_s * pow(u,2));
double u_star = pow(tau / rho_s, 1/2);
cout << "Enter the required Y+ value: ";
cin >> Y_Plus;
y = (Y_Plus * air.mu) / (air.rho_s * u_star);
return y;
}

int main(){
    fluid air;
    geometry geo;
    mesh msh;
    cout << "\t\t\t\t\t\tWelcome\n";
    cout << "\t\t\t\tComputational Fluid Dynamics (CFD) Calculator\n\n";
    cout << "1. Reynolds Number Calculator \n2. Y+ Calculator\n3. Boundary Layer Thickness Calculator"
         << "\n4. Boundary Layer Mesh Calculator\n5. Turbulence Property Calculator\n";
    cout << "\nEnter your choice for calculation: ";
    cin >> t_no;
    switch (t_no){
        case 1:
            cout << "\nYou have choose Reynolds Number Calculator\n";
            cout << "Enter the types of flow (1 - Internal Flow) (2 - External Flow): ";
            cin >> t_no;
            switch(t_no){
            case 1:
                cout << "You have chosen Internal Flow. \n";
                    ReNum = RenoldsNumber_Internal();
                    cout << "\nFlow Reynolds number = " << ReNum << endl;
                    if(ReNum <= 2300){
                        cout << "Flow is Laminar.\n";
                        }
                    else if(ReNum < 4000){
                        cout << "Flow is Transition.\n";
                        }
                    else{
                        cout << "Flow is Turbulent.\n";
                        }
                break;
            case 2:
                cout << "You have chosen External Flow. \n";
                    ReNum = RenoldsNumber_External();
                    cout << "\nFlow Reynolds number = " << ReNum << endl;
                    if(ReNum <= 2e5){
                        cout << "Flow is Laminar.\n";
                        }
                    else if(ReNum < 5e5){
                        cout << "Flow is Transition.\n";
                        }
                    else{
                        cout << "Flow is Turbulent.\n";
                        }
                break;
            default:
                cout << "You have chosen incorrect option. \n";
            }
        break;
        case 2:
            cout << "You have choose Y+ calculator.\n";
            cout << "\nEnter the types of flow (1 - Internal Flow) (2 - External Flow): ";
            cin >> t_no;
            switch(t_no){
            case 1:
                cout << "You have chosen Internal Flow. \n";
                y = Y_plusCalc_Internal();
                cout << "The First cell height (Y) is " << y << " m\n";
                break;
            case 2:
                cout << "You have chosen External Flow. \n";
                y = Y_plusCalc_External();
                cout << "The First cell height (Y) is " << y << " m\n";
                break;
            default:
                cout << "You have chosen incorrect option. \n";
            }
        break;
        case 3:
            cout << "You have choose Boundary Layer Thickness calculator. \n";
            cout << "Enter the types of flow (1 - Internal Flow) (2 - External Flow): ";
            cin >> t_no;
                switch (t_no){
                        case 1:
                            cout << "You have chosen Internal Flow. \n";
                                ReNum = RenoldsNumber_Internal();
                                cout << "\nFlow Reynolds number = " << ReNum << endl;
                                if(ReNum <= 2300){
                                    cout << "Flow is Laminar.\n";
                                    cout << "Laminar boundary layer calculator.\n" << endl;
                                    cout << "Total length up to which boundary layer thickness to be determined (in meters): ";
                                    cin >> x;
                                    cout << "\n\nBoundary layer thickness (Delta): " << (5 * x)/(pow(ReNum,0.5)) << " meters" << endl;
                                    cout << "\nDisplacement thickness (Delta_star): " << (1.72 * x)/(pow(ReNum,0.5)) << " meters" << endl;
                                    cout << "\nMomentum thickness (Theta_star): " << (0.664 * x)/(pow(ReNum,0.5)) << " meters" << endl;
                                    cout << "\nShape Factor (H): " << ((1.72 * x)/(pow(ReNum,0.5)) / (0.664 * x)/(pow(ReNum,0.5))) << endl;
                                    cout << "\nLocal skin friction coefficient (c_f): " << (0.664)/(pow(ReNum,0.5)) << endl;
                                    cout << "\nDrag coefficient (c_d): " << (1.328)/(pow(ReNum,0.5)) << endl;
                                    }
                                else{
                                    cout << "Flow is Turbulent.\n";
                                    cout << "Turbulent boundary layer calculator.\n" << endl;
                                    cout << "Total length up to which boundary layer thickness to be determined (in meters): ";
                                    cin >> x;
                                    cout << "\n\nBoundary layer thickness (Delta): " << (0.16 * x)/(pow(ReNum,1/7)) << " meters" << endl;
                                    cout << "\nDisplacement thickness (Delta_star): " << (0.02 * x)/(pow(ReNum,1/7)) << " meters" << endl;
                                    cout << "\nMomentum thickness (Theta_star): " << (0.016 * x)/(pow(ReNum,1/7)) << " meters" << endl;
                                    cout << "\nShape Factor (H): " << ((0.02 * x)/(pow(ReNum,1/7)) / (0.016 * x)/(pow(ReNum,1/7))) << endl;
                                    cout << "\nLocal skin friction coefficient (c_f): " << (0.027)/(pow(ReNum,1/7)) << endl;
                                    cout << "\nDrag coefficient (c_d): " << (0.031)/(pow(ReNum,1/7)) << endl;
                                    }
                            break;
                        case 2:
                            cout << "You have chosen External Flow. \n";
                                ReNum = RenoldsNumber_External();
                                cout << "\nFlow Reynolds number = " << ReNum << endl;
                                if(ReNum <= 2e5){
                                    cout << "Flow is Laminar.\n";
                                    cout << "Laminar boundary layer calculator.\n" << endl;
                                    cout << "Total length up to which boundary layer thickness to be determined (in meters): ";
                                    cin >> x;
                                    cout << "\n\nBoundary layer thickness (Delta): " << (5 * x)/(pow(ReNum,0.5)) << " meters" << endl;
                                    cout << "\nDisplacement thickness (Delta_star): " << (1.72 * x)/(pow(ReNum,0.5)) << " meters" << endl;
                                    cout << "\nMomentum thickness (Theta_star): " << (0.664 * x)/(pow(ReNum,0.5)) << " meters" << endl;
                                    cout << "\nShape Factor (H): " << ((1.72 * x)/(pow(ReNum,0.5)) / (0.664 * x)/(pow(ReNum,0.5))) << endl;
                                    cout << "\nLocal skin friction coefficient (c_f): " << (0.664)/(pow(ReNum,0.5)) << endl;
                                    cout << "\nDrag coefficient (c_d): " << (1.328)/(pow(ReNum,0.5)) << endl;

                                    }
                                else{
                                    cout << "Flow is Turbulent.\n";
                                    cout << "Turbulent boundary layer calculator.\n" << endl;
                                    cout << "Total length up to which boundary layer thickness to be determined (in meters): ";
                                    cin >> x;
                                    cout << "\n\nBoundary layer thickness (Delta): " << (0.16 * x)/(pow(ReNum,1/7)) << " meters" << endl;
                                    cout << "\nDisplacement thickness (Delta_star): " << (0.02 * x)/(pow(ReNum,1/7)) << " meters" << endl;
                                    cout << "\nMomentum thickness (Theta_star): " << (0.016 * x)/(pow(ReNum,1/7)) << " meters" << endl;
                                    cout << "\nShape Factor (H): " << ((0.02 * x)/(pow(ReNum,1/7)) / (0.016 * x)/(pow(ReNum,1/7))) << endl;
                                    cout << "\nLocal skin friction coefficient (c_f): " << (0.027)/(pow(ReNum,1/7)) << endl;
                                    cout << "\nDrag coefficient (c_d): " << (0.031)/(pow(ReNum,1/7)) << endl;
                                    }
                            break;
                        default:
                            cout << "You have chosen incorrect option. \n";
                 }
        break;
        case 4:
            cout << "\nYou have choose Boundary Layer Mesh calculator. \n";
            cout << "\nEnter the types of expansion (1 - Linear) (2 - Exponential): ";
            cin >> t_no;
            switch(t_no){
                case 1:
                    cout << "\nYou have choose Linear expansion Boundary layer calculator\n";
                    cout << "\nEnter the first cell height: ";
                    cin >> msh.init_ht;
                    cout << "\nEnter the expansion ratio (Default = 1.2): ";
                    cin >> msh.er;
                    cout << "\nEnter the number of boundary layer required: ";
                    cin >> msh.nl;
                    msh.tot_ht = msh.nl * msh.init_ht * ( ( ( msh.nl - 1) * ( msh.er - 1) + 2 ) / 2);
                    msh.lst_cell_ht = msh.init_ht * (1 + (msh.nl - 1) * (msh.er - 1));
                    cout << "\nTotal thickness of the Mesh Boundary layer is " << msh.tot_ht << " m and"
                         << " The last cell height is " << msh.lst_cell_ht << " m.\n";
                break;
                case 2:
                    cout << "\nYou have choose Exponential expansion Boundary layer calculator\n";
                    cout << "\nEnter the first cell height: ";
                    cin >> msh.init_ht;
                    cout << "\nEnter the expansion ratio ( Default = 1.2 ): ";
                    cin >> msh.er;
                    cout << "\nEnter the number of boundary layer required: ";
                    cin >> msh.nl;
                    msh.tot_ht = msh.init_ht * ((1 - pow(msh.er, msh.nl) / (1 - msh.er)));
                    msh.lst_cell_ht =  msh.init_ht * pow(msh.er, (msh.nl - 1));
                    cout << "\nTotal thickness of the Mesh Boundary layer is " << msh.tot_ht << " m and"
                         << " The last cell height is " << msh.lst_cell_ht << " m.\n";
                break;
                default:
                cout << "You have chosen incorrect option.\n";
            }
        break;
        case 5:
            cout << "\nYou have choose Turbulence property calculator. \n";
            cout << "\nChoose your turbulence model (1. Spalart-Allmaras) (2. k-epsilon) (3. k-omega): ";
            cin >> t_no;
            switch (t_no){
            case 1:
                cout << "You have choose Spalart-Allmaras Turbulent Model";
                cout << "\n\nEnter the types of flow (1 - Internal Flow) (2 - External Flow): ";
                cin >> t_no;
                switch (t_no){
                        case 1:
                            cout << "You have chosen Internal Flow. \n";
                                cout << "\nEnter the static pressure of the fluid (N/m^2) (Default = 101325): ";
                                cin >> air.P_s;
                                cout << "\nEnter the static temperature of the fluid (Kelvin) (Default = 300): ";
                                cin >> air.T_s;
                                cout << "\nEnter the average velocity of the fluid (m/s) (Default = 1): ";
                                cin >> air.u_avg;
                                cout << "\nEnter the characteristic diameter (m) (Default = 1): ";
                                cin >> geo.d;
                                cout << "\nEnter the dynamic viscosity of the fluid (kg/m-s) (Default = 1.875e-5): ";
                                cin >> air.mu;
                                ReNum = ((air.P_s / (air.R * air.T_s)) * air.u_avg * geo.d) / air.mu;
                                cout << "\nEnter the type of Cross Section (1 - Circular) (2 - Square) (3 - Rectangular): ";
                                cin >> t_no;
                                switch (t_no){
                                    case 1:
                                        cout << "\nYou have choose Circular cross section";
                                        geo.D_h = HydraulicDiameter_Circular();
                                    break;
                                    case 2:
                                        cout << "\nYou have choose Square cross section";
                                        geo.D_h = HydraulicDiameter_Square();
                                    break;
                                    case 3:
                                        cout << "\nYou have choose Rectangular cross section";
                                        geo.D_h = HydraulicDiameter_Rectangular();
                                    break;
                                    default:
                                        cout << "You have chosen incorrect option.\n";
                                }
                                Ti = 0.16 * pow(ReNum,-0.125);
                                Tls = (0.07 * geo.D_h) / pow(c_mu, 0.75);
                                mtv = (c_mu * pow(1.5,0.5) * air.u_avg * Ti * Tls);
                                tvr = mtv / air.mu;
                                cout << "\nFlow Reynolds number = " << ReNum << endl;
                                cout << "\nHydraulic Diameter (m) = " << geo.D_h << endl;
                                cout << "\nModified Turbulent Viscosity = " << mtv << endl;
                                cout << "\nTurbulent Intensity = " << Ti << endl;
                                cout << "\nTurbulent Length Scale (m) = " << Tls << endl;
                                cout << "\nTurbulent Viscosity Ratio = " << tvr << endl;
                        break;
                        case 2:
                            cout << "You have chosen External Flow. \n";
                                cout << "\nEnter the static pressure of the fluid (N/m^2) (Default = 101325): ";
                                cin >> air.P_s;
                                cout << "\nEnter the static temperature of the fluid (Kelvin) (Default = 300): ";
                                cin >> air.T_s;
                                cout << "\nEnter the average velocity of the fluid (m/s) (Default = 1): ";
                                cin >> air.u_avg;
                                cout << "\nEnter the characteristic Length (m) (Default = 1): ";
                                cin >> geo.D_h;
                                cout << "\nEnter the dynamic viscosity of the fluid (kg/m-s) (Default = 1.875e-5): ";
                                cin >> air.mu;
                                ReNum = ((air.P_s / (air.R * air.T_s)) * air.u_avg * geo.d) / air.mu;
                                Ti = 0.16 * pow(ReNum,-0.125);
                                Tls = (0.07 * geo.D_h) / pow(c_mu, 0.75);
                                mtv = (c_mu * pow(1.5,0.5) * air.u_avg * Ti * Tls);
                                tvr = mtv / air.mu;
                                cout << "\nFlow Reynolds number = " << ReNum << endl;
                                cout << "\nHydraulic Diameter or Length (m) = " << geo.D_h << endl;
                                cout << "\nModified Turbulent Viscosity = " << mtv << endl;
                                cout << "\nTurbulent Intensity = " << Ti << endl;
                                cout << "\nTurbulent Length Scale (m) = " << Tls << endl;
                                cout << "\nTurbulent Viscosity Ratio = " << tvr << endl;
                        break;
                        default:
                            cout << "You have chosen incorrect option. \n";
                }
            break;
            case 2:
                cout << "You have choose K-epsilon Turbulent Model";
                cout << "\n\nEnter the types of flow (1 - Internal Flow) (2 - External Flow): ";
                cin >> t_no;
                switch (t_no){
                        case 1:
                            cout << "You have chosen Internal Flow. \n";
                            cout << "\nEnter the static pressure of the fluid (N/m^2) (Default = 101325): ";
                            cin >> air.P_s;
                            cout << "\nEnter the static temperature of the fluid (Kelvin) (Default = 300): ";
                            cin >> air.T_s;
                            cout << "\nEnter the average velocity of the fluid (m/s) (Default = 1): ";
                            cin >> air.u_avg;
                            cout << "\nEnter the characteristic diameter (m) (Default = 1): ";
                            cin >> geo.d;
                            cout << "\nEnter the dynamic viscosity of the fluid (kg/m-s) (Default = 1.875e-5): ";
                            cin >> air.mu;
                            ReNum = ((air.P_s / (air.R * air.T_s)) * air.u_avg * geo.d) / air.mu;
                            cout << "\n\nEnter the type of Cross Section (1 - Circular) (2 - Square) (3 - Rectangular): ";
                            cin >> t_no;
                            switch (t_no){
                                case 1:
                                    cout << "\nYou have choose Circular cross section";
                                    geo.D_h = HydraulicDiameter_Circular();
                                break;
                                case 2:
                                    cout << "\nYou have choose Square cross section";
                                    geo.D_h = HydraulicDiameter_Square();
                                break;
                                case 3:
                                    cout << "\nYou have choose Rectangular cross section";
                                    geo.D_h = HydraulicDiameter_Rectangular();
                                break;
                                default:
                                    cout << "You have chosen incorrect option.\n";
                            }
                            Ti = 0.16 * pow(ReNum,-0.125);
                            ke = 1.5 * pow(air.u_avg * Ti , 2);
                            Tls = (0.07 * geo.D_h) / pow(c_mu, 0.75);
                            e = pow(ke,1.5) / Tls;
                            tv = ((air.P_s / (air.R * air.T_s) * 0.09 *ke *ke ) / e );
                            tvr = tv / air.mu;
                            cout << "\nFlow Reynolds number = " << ReNum << endl;
                            cout << "\nHydraulic Diameter (m) = " << geo.D_h << endl;
                            cout << "\nTurbulent Viscosity (kg/m-s) = " << tv << endl;
                            cout << "\nTurbulent Intensity = " << Ti << endl;
                            cout << "\nTurbulent Length Scale (m) = " << Tls << endl;
                            cout << "\nTurbulent Viscosity Ratio = " << tvr << endl;
                            cout << "\nTurbulent Kinetic Energy (m2/s2) = " << ke << endl;
                            cout << "\nTurbulent Dissipation Rate (m2/s3) = " << e << endl;
                        break;
                        case 2:
                            cout << "You have chosen External Flow. \n";
                            cout << "\nEnter the static pressure of the fluid (N/m^2) (Default = 101325): ";
                            cin >> air.P_s;
                            cout << "\nEnter the static temperature of the fluid (Kelvin) (Default = 300): ";
                            cin >> air.T_s;
                            cout << "\nEnter the average velocity of the fluid (m/s) (Default = 1): ";
                            cin >> air.u_avg;
                            cout << "\nEnter the characteristic length (m) (Default = 1): ";
                            cin >> geo.d;
                            cout << "\nEnter the dynamic viscosity of the fluid (kg/m-s) (Default = 1.875e-5): ";
                            cin >> air.mu;
                            ReNum = ((air.P_s / (air.R * air.T_s)) * air.u_avg * geo.d) / air.mu;
                            Ti = 0.16 * pow(ReNum,-0.125);
                            ke = 1.5 * pow(air.u_avg * Ti , 2);
                            Tls = (0.07 * geo.D_h) / pow(c_mu, 0.75);
                            e = pow(ke,1.5) / Tls;
                            tv = ((air.P_s / (air.R * air.T_s) * 0.09 *ke *ke ) / e );
                            tvr = tv / air.mu;
                            cout << "\nFlow Reynolds number = " << ReNum << endl;
                            cout << "\nHydraulic Diameter (m) = " << geo.D_h << endl;
                            cout << "\nTurbulent Viscosity (kg/m-s) = " << tv << endl;
                            cout << "\nTurbulent Intensity = " << Ti << endl;
                            cout << "\nTurbulent Length Scale (m) = " << Tls << endl;
                            cout << "\nTurbulent Viscosity Ratio = " << tvr << endl;
                            cout << "\nTurbulent Kinetic Energy (m2/s2) = " << ke << endl;
                            cout << "\nTurbulent Dissipation Rate (m2/s3) = " << e << endl;
                        break;
                        default:
                             cout << "You have chosen incorrect option.\n";
                        }
            break;
            case 3:
                cout << "You have choose K-omega Turbulent Model";
                cout << "\n\nEnter the types of flow (1 - Internal Flow) (2 - External Flow): ";
                cin >> t_no;
                switch (t_no){
                        case 1:
                            cout << "You have chosen Internal Flow. \n";
                            cout << "\nEnter the static pressure of the fluid (N/m^2) (Default = 101325): ";
                            cin >> air.P_s;
                            cout << "\nEnter the static temperature of the fluid (Kelvin) (Default = 300): ";
                            cin >> air.T_s;
                            cout << "\nEnter the average velocity of the fluid (m/s) (Default = 1): ";
                            cin >> air.u_avg;
                            cout << "\nEnter the characteristic diameter (m) (Default = 1): ";
                            cin >> geo.d;
                            cout << "\nEnter the dynamic viscosity of the fluid (kg/m-s) (Default = 1.875e-5): ";
                            cin >> air.mu;
                            ReNum = ((air.P_s / (air.R * air.T_s)) * air.u_avg * geo.d) / air.mu;
                            cout << "\n\nEnter the type of Cross Section (1 - Circular) (2 - Square) (3 - Rectangular): ";
                            cin >> t_no;
                            switch (t_no){
                                case 1:
                                    cout << "\nYou have choose Circular cross section";
                                    geo.D_h = HydraulicDiameter_Circular();
                                break;
                                case 2:
                                    cout << "\nYou have choose Square cross section";
                                    geo.D_h = HydraulicDiameter_Square();
                                break;
                                case 3:
                                    cout << "\nYou have choose Rectangular cross section";
                                    geo.D_h = HydraulicDiameter_Rectangular();
                                break;
                                default:
                                    cout << "You have chosen incorrect option.\n";
                            }
                            Ti = 0.16 * pow(ReNum,-0.125);
                            ke = 1.5 * pow(air.u_avg * Ti , 2);
                            Tls = (0.07 * geo.D_h) / pow(c_mu, 0.75);
                            e = pow(ke,1.5) / Tls;
                            tv = ((air.P_s / (air.R * air.T_s) * 0.09 *ke *ke ) / e );
                            tvr = tv / air.mu;
                            om = (pow(ke,0.5) / (c_mu *Tls));
                            cout << "\nFlow Reynolds number = " << ReNum << endl;
                            cout << "\nHydraulic Diameter (m) = " << geo.D_h << endl;
                            cout << "\nTurbulent Viscosity = " << tv << endl;
                            cout << "\nTurbulent Intensity = " << Ti << endl;
                            cout << "\nTurbulent Length Scale (m) = " << Tls << endl;
                            cout << "\nTurbulent Viscosity Ratio = " << tvr << endl;
                            cout << "\nTurbulent Kinetic Energy (m2/s2) = " << ke << endl;
                            cout << "\nSpecific Dissipation Rate (m2/s2) = " << tvr << endl;
                        break;
                        case 2:
                            cout << "You have chosen External Flow. \n";
                            cout << "\nEnter the static pressure of the fluid (N/m^2) (Default = 101325): ";
                            cin >> air.P_s;
                            cout << "\nEnter the static temperature of the fluid (Kelvin) (Default = 300): ";
                            cin >> air.T_s;
                            cout << "\nEnter the average velocity of the fluid (m/s) (Default = 1): ";
                            cin >> air.u_avg;
                            cout << "\nEnter the characteristic length (m) (Default = 1): ";
                            cin >> geo.D_h;
                            cout << "\nEnter the dynamic viscosity of the fluid (kg/m-s) (Default = 1.875e-5): ";
                            cin >> air.mu;
                            ReNum = ((air.P_s / (air.R * air.T_s)) * air.u_avg * geo.d) / air.mu;
                            Ti = 0.16 * pow(ReNum,-0.125);
                            ke = 1.5 * pow(air.u_avg * Ti , 2);
                            Tls = (0.07 * geo.D_h) / pow(c_mu, 0.75);
                            e = pow(ke,1.5) / Tls;
                            tv = ((air.P_s / (air.R * air.T_s) * 0.09 *ke *ke ) / e );
                            tvr = tv / air.mu;
                            om = (pow(ke,0.5) / (c_mu *Tls));
                            cout << "\nFlow Reynolds number = " << ReNum << endl;
                            cout << "\nHydraulic Diameter (m) = " << geo.D_h << endl;
                            cout << "\nTurbulent Viscosity = " << tv << endl;
                            cout << "\nTurbulent Intensity = " << Ti << endl;
                            cout << "\nTurbulent Length Scale (m) = " << Tls << endl;
                            cout << "\nTurbulent Viscosity Ratio = " << tvr << endl;
                            cout << "\nTurbulent Kinetic Energy (m2/s2) = " << ke << endl;
                            cout << "\nSpecific Dissipation Rate (m2/s2) = " << tvr << endl;
                        break;
                        default:
                             cout << "You have chosen incorrect option.\n";
                        }
            break;
            default:
            cout << "You have chosen incorrect option.\n";
            }
        break;
        default:
            cout << "You have chosen incorrect option.\n";
    }
return 0;
}
