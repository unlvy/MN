#include <iostream>
#include <cmath>
#include <vector>

double eps = 1.;
double delta = 0.1;
const int nx = 150;
const int ny = 100;

double V1 = 10.;
double V2 = 0.;

double xmax = delta*nx;
double ymax = delta*ny;

double sigmax = 0.1 * xmax;
double sigmay = 0.1 * ymax;

double TOL = std::pow(10,-8);

double ro1(double x, double y)
{
    double result = -1.*std::pow((x-0.35*xmax),2)/(sigmax*sigmax) - 1.*std::pow((y-0.5*ymax),2)/(sigmay*sigmay);
    return std::exp(result);
}

double ro2(double x, double y)
{
    double result = -1.*std::pow((x-0.65*xmax),2)/(sigmax*sigmax) - 1.*std::pow((y-0.5*ymax),2)/(sigmay*sigmay);
    return -std::exp(result);
}

int main()
{
    std::vector<double> omegaG = {0.6/*, 1.0*/};

    for (double omega : omegaG)
    {
        std::string vs = "dane/g_v" + std::to_string(omega) + ".dat";
        std::string ss = "dane/g_s" + std::to_string(omega) + ".dat";
        FILE *fV = fopen ("wynik.txt","w");
//        FILE *fS = fopen (ss.c_str(),"w");

        double Vn [nx+5][nx+5] = {0.};
        double P [nx+5][ny+5] = {0.};
        double Vs [nx+5][ny+5] = {0.};

        for(int i=0; i<nx+1; i++)
        {
            Vs[i][0] = V1;
            Vn[i][0] = V1;
            for(int j=0; j<ny+1; j++)
            {
                P[i][j] = ro1(i*delta, j*delta) + ro2(i*delta, j*delta);
            }
        }

        int count = 0;
        double S[50000] = {0.};
        S[0] = 0.;

        do
        {
            count++;
            for(int i=1; i<nx; i++)
            {
                for(int j = 1; j<ny; j++)
                {
                    Vn[i][j] = 1/4. * (Vs[i+1][j] + Vs[i-1][j] + Vs[i][j+1] + Vs[i][j-1] + (delta*delta)/eps * P[i][j]);
                }
            }
            for(int j=1; j<ny+1; j++)
            {
                Vn[0][j] = Vn[1][j];
                Vn[nx][j] = Vn[nx-1][j];
            }
            for(int i=0; i<nx+1; i++)
            {
                for(int j=0; j<ny+1; j++)
                {
                    Vs[i][j]=(1.-omega) * Vs[i][j] + omega * Vn[i][j];
                }
            }
            for(int i=0; i<nx; i++)
            {
                for(int j=0; j<ny; j++)
                {
                    S[count] += (delta*delta) * (0.5 * std::pow((Vn[i+1][j] - Vn[i][j])/delta, 2) + 0.5*std::pow((Vn[i][j+1]-Vn[i][j])/delta, 2 ) - (P[i][j]*Vn[i][j]) );
                }
            }
            //fprintf(fS,"%d %f\n",count-1,S[count]);
        }while(std::abs((S[count] - S[count-1]) / S[count-1]) > TOL);

        double Verr[nx+6][nx+6] = {0.};

        for(int i=1; i<nx; i++)
        {
            for(int j=1; j<ny; j++)
            {

                Verr[i][j] = (Vn[i+1][j] - 2.*Vn[i][j] + Vn[i-1][j]) / (delta*delta) + (Vn[i][j+1] - 2.*Vn[i][j] + Vn[i][j-1])/(delta*delta) + P[i][j]/eps;
            }
        }

        for(int i=0; i<=nx; i++)
        {
            for(int j=0; j<=ny; j++)
            {
                //fprintf(fV,"%f %f %f %f\n",i*delta, j*delta, Verr[i][j], Vn[i][j]);
                fprintf(fV, "%lf ", P[i][j]);
            }
            fprintf(fV, "\n");
        }
       fclose (fV);
        //fclose (fS);
    }

//    std::vector<double> omegaL = {1.0, 1.4, 1.8, 1.9};
//
//    for (double omega: omegaL)
//    {
//        std::string s = "dane/l_s" + std::to_string(omega) + ".dat";
//        FILE *f = fopen(s.c_str(),"w");
//        double V [nx+5][ny+5] = {0.};
//        double P [nx+5][ny+5] = {0.};
//        for(int i=0; i<nx+1; i++)
//        {
//            for(int j=0; j<ny+1; j++)
//            {
//                V[i][0] = V1;
//                P[i][j] = ro1(i*delta, j*delta) + ro2(i*delta, j*delta);
//            }
//        }
//
//        int count = 0;
//        double S[100000] = {0.};
//        S[0] = 0.;
//
//        do
//        {
//            count++;
//
//            for(int i=1; i<nx; i++)
//            {
//                for(int j=1; j<ny; j++)
//                {
//                    V[i][j] = (1-omega) * V[i][j] + (omega*0.25)*(V[i+1][j] + V[i-1][j] + V[i][j+1] + V[i][j-1] + (delta*delta)/eps * P[i][j]);
//                }
//            }
//
//            for(int j=1; j<ny; j++)
//            {
//                V[0][j] = V[1][j];
//                V[nx][j] = V[nx-1][j];
//            }
//
//            for(int i=0; i<nx; i++)
//            {
//                for(int j=0; j<ny; j++)
//                {
//                    S[count] += (delta*delta) * (0.5 * std::pow((V[i+1][j] - V[i][j])/delta, 2) + 0.5*std::pow((V[i][j+1]-V[i][j])/delta, 2 ) - (P[i][j]*V[i][j]) );
//                }
//            }
//
//            fprintf(f,"%d %f\n",count-1,S[count]);
//
//        }while(std::abs((S[count] - S[count-1]) / S[count-1]) > TOL);
//
//    }

    return 0;
}