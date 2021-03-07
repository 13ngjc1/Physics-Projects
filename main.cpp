#include <iostream>
#include <cstdlib>
#include <fstream>
#include <ctime>
#include <cmath>

#define lightspeed 2997924.58
#define Ar 39.938
#define u 9.314e8
#define Kb 8.6173324e-5
#define sigma 3.40 //length in angstrom
#define epsilon 1.03e-2 //energy in eV
#define PI 3.14159265

//declaration of constants
double const temp = 300; //in K
int const pAxis = 5;
int const particleCount = pAxis * pAxis * pAxis;
double const timestep = 0.01;
double const separation = 5; //in angstrom
int const sep = 5; //in angstrom, for int calculations
int const boxSize = sep*(pAxis+1);

//heat bath correction
double const tau = 0.9;
double const tempBath = 250; //in K

//mass
double mass = (Ar*u/pow(lightspeed,2));

using namespace std;

//initialization
void initialPosition(double[], double[], double[]);
void initialVelocity(double[], double[], double[]);
double BoxMullerTransform();

//acceleration function/force calculation function
void forcecalculation(double[], double[], double[], double[], double[], double[]);

int main()
{
    //initialize position and velocity arrays
    double x[particleCount], y[particleCount], z[particleCount];
    double Vx[particleCount] = {0}, Vy[particleCount] = {0}, Vz[particleCount] = {0};

    //acceleration arrays
    double ax[particleCount] = {0}, ay[particleCount] = {0}, az[particleCount] = {0};
    double ax2[particleCount] = {0}, ay2[particleCount] = {0}, az2[particleCount] = {0};

    //autocorrelation data struct
    //[i][j], where i = i-th particle and j = time array
    double xA[particleCount][501], yA[particleCount][501], zA[particleCount][501];

    //file pointers
    ofstream vmd;
    ofstream EnergyRec;
    ofstream autocor; //autocorrelation
    vmd.open("coord.xyz");
    EnergyRec.open("Energy Record.csv");
    EnergyRec<<"KE"<<","<<"PE"<<","<<"Total Energy"<<","<<"Temperature"<<","<<endl;
    autocor.open("autocorrelation.csv");

    //initialize velocity and positions
    initialPosition(x, y, z);
    initialVelocity(Vx, Vy, Vz);

    //heat bath correction factor
    double heatBath;

    //adjust CoM
    double VxTotal, VyTotal, VzTotal;
    for(int i=0; i<particleCount; i++)
    {
        VxTotal += Vx[i];
        VyTotal += Vy[i];
        VzTotal += Vz[i];
    }

    //cout<<VxTotal<<" "<<VyTotal<<" "<<VzTotal<<" "<<endl;

    for(int j=0; j<particleCount; j++)
    {
        Vx[j] = Vx[j] - VxTotal/particleCount;
        Vy[j] = Vy[j] - VyTotal/particleCount;
        Vz[j] = Vz[j] - VzTotal/particleCount;

        //cout<<Vx[j]<<" "<<Vy[j]<<" "<<Vz[j]<<endl;
    }

    double potentialShift = 4*1.03e-2*(pow((3.4/(2.5*sigma)),12) - pow((3.4/(2.5*sigma)),6));

    //for calculating energy later
    double KE, PE, TotalEnergy, T;

    //run force calculation once to ensure acceleration array is zero
    //and to ensure there is acceleration at time loop
    forcecalculation(x, y, z, ax, ay, az);

    int count = 0;
    int timecount;
    for(timecount = 0; timecount<5000; timecount++) //time evolution
    {
         //time evolution
        for(int p=0; p<particleCount; p++)
        {
            x[p] += Vx[p]*timestep+0.5*ax[p]*timestep*timestep;
            y[p] += Vy[p]*timestep+0.5*ay[p]*timestep*timestep;
            z[p] += Vz[p]*timestep+0.5*az[p]*timestep*timestep;

            //boundary check
            if(x[p] < 0)
            {
                x[p] = -x[p];
                Vx[p] = -Vx[p];
            }
            if(x[p] > boxSize)
            {
                x[p] = boxSize - (x[p] - boxSize);
                Vx[p] = -Vx[p];
            }
            if(y[p] < 0)
			{
				y[p] = -y[p];
				Vy[p] = -Vy[p];
			}
			if (y[p] > boxSize)
			{
				y[p] = boxSize-(y[p] - boxSize);
				Vy[p] = -Vy[p];
			}
			if(z[p] < 0)
			{
				z[p] = -z[p];
				Vz[p] = -Vz[p];
			}
			if(z[p] > boxSize)
			{
				z[p] = boxSize-(z[p]-boxSize);
				Vz[p] = -Vz[p];
			}

        }

        //two acceleration arrays needed for later calculations
        for(int i=0; i<particleCount; i++)
        {
            ax2[i] = ax[i];
            ay2[i] = ay[i];
            az2[i] = az[i];
        }

        forcecalculation(x, y, z, ax, ay, az);

        //update velocity
        for(int m=0; m<particleCount; m++)
        {
            Vx[m] += 0.5*(ax[m]+ax2[m])*timestep;
            Vy[m] += 0.5*(ay[m]+ay2[m])*timestep;
            Vz[m] += 0.5*(az[m]+az2[m])*timestep;
        }

        //calculate KE
        KE = 0;
        double v_2 = 0;
        for(int vcount = 0; vcount < particleCount; vcount++)
        {
            v_2 += pow(Vx[vcount], 2) + pow(Vy[vcount], 2) + pow(Vz[vcount], 2);
        }
        KE = (v_2*mass)/2;

        //calculating temperature
        T = (v_2*mass)/(3*Kb*particleCount);

        //heat bath correction factor
        heatBath = 0;
        heatBath = sqrt(1+(timestep/tau)*((tempBath/T)-1));
        //if(timecount % 100 == 0)
            //cout<<heatBath<<endl;

        //Heat bath correction
        for(int c=0; c<particleCount; c++)
        {
            Vx[c] = heatBath*Vx[c];
            Vy[c] = heatBath*Vy[c];
            Vz[c] = heatBath*Vz[c];
        }


        //calculating PE
        PE = 0;
        for(int k=0; k<particleCount; k++)
        {
            for(int l=k+1; l<particleCount-1; l++)
            {
                double d_2 = pow( (x[k]-x[l]), 2)+pow((y[k]-y[l]),2)+pow((z[k]-z[l]),2);
                double dist = pow((sigma*sigma/d_2), 3);

                if(d_2 <= 2.5*2.5*sigma*sigma)
                {
                    double P = 4*epsilon*(dist*dist-dist)-potentialShift;
                    PE += P;
                }
            }
        }

        //calculating total energy
        TotalEnergy = 0;
        TotalEnergy = KE + PE;

        //outputting to csv file
        EnergyRec<<KE<<","<<PE<<","<<TotalEnergy<<","<<T<<endl;

        //outputting every 100 frames
        if(timecount%100 == 0)
        {
            vmd<<particleCount<<endl;
            vmd<<"Frame "<<count+1<<endl;

            for(int pcount=0; pcount<particleCount; pcount++)
            {
                vmd<<"Ar "<<x[pcount]<<" "<<y[pcount]<<" "<<z[pcount]<<endl;
            }
            count++;
        }

        //autocorrelation position update - storing the positions
        if(timecount %10 == 0)
        {
            for(int k=0; k<particleCount; k++)
            {
                xA[k][timecount/11] = x[k];
                yA[k][timecount/11] = y[k];
                zA[k][timecount/11] = z[k];
            }
        }
    }

    double MSDTemp = 0;
    double MSD= 0;
    //MSD calc of autocorrelation
    for(int i=1; i<501; i++)
    {
        int n=0;
        MSDTemp = 0;
        for(int j=0; j<particleCount; j++)
        {
            for(int k=0; k<501-i; i++)
            {
                double x2 = xA[j][i+k] - xA[j][k];
                double y2 = yA[j][i+k] - yA[j][k];
                double z2 = zA[j][i+k] - zA[j][k];

                MSDTemp += x2*x2 + y2*y2 + z2*z2;
                n++;
            }
        }
        MSD = MSDTemp/n;
        autocor<<i<<MSD<<endl;
        cout<<i<<MSD<<endl;
    }
}


void initialPosition(double xPos[], double yPos[], double zPos[])
{
    int pIndex = 0;
    for(int i=0; i<pAxis; i++)
    {
        for(int j=0; j<pAxis; j++)
        {
            for(int k=0; k<pAxis; k++)
            {
                xPos[pIndex] = separation*i + (sep/2);
                yPos[pIndex] = separation*j + (sep/2);
                zPos[pIndex] = separation*k + (sep/2);

                //cout << xPos[pIndex] <<" "<< yPos[pIndex]<<" " << zPos[pIndex] << endl;
                pIndex++;
            }
        }
    }
    cout<<"Positions initialized"<<endl;
}

void initialVelocity(double xv[], double yv[], double zv[])
{
     for(int indexcount=0; indexcount<particleCount; indexcount++)
     {
         xv[indexcount] = BoxMullerTransform();
         yv[indexcount] = BoxMullerTransform();
         zv[indexcount] = BoxMullerTransform();
         //cout << xv[indexcount] <<" "<< yv[indexcount]<<" " << zv[indexcount] << endl;
     }
     cout<<"Velocities initialized"<<endl;
}

double BoxMullerTransform()
{
    double y1 = (double)rand()/(double)RAND_MAX;
	double y2 = (double)rand()/(double)RAND_MAX;
	double x1 = sqrt(-2*log(y1))*cos(2*PI*y2) * sqrt(Kb*temp/(Ar*u)) * lightspeed;

    return x1;
}

void forcecalculation(double px[], double py[], double pz[], double ax[], double ay[], double az[])
{
    for(int i=0; i<particleCount; i++)
    {
        ax[i] = 0;
        ay[i] = 0;
        az[i] = 0;
    }
    for(int j=0; j<particleCount-1; j++)
    {
        for(int k=j+1; k<particleCount; k++)
        {
            double d_2 = pow((px[k]-px[j]) ,2) + pow( (py[k]-py[j]) ,2) + pow( (pz[k]-pz[j]), 2);
		    double D = (sigma*sigma/d_2)*(sigma*sigma/d_2)*(sigma*sigma/d_2);

			if(d_2 <= 2.5*2.5*sigma*sigma)
			{
				double fx = 24*epsilon*(px[j]-px[k])*(2*D*D-D)*(1/d_2);
                double fy = 24*epsilon*(py[j]-py[k])*(2*D*D-D)*(1/d_2);
                double fz = 24*epsilon*(pz[j]-pz[k])*(2*D*D-D)*(1/d_2);

				ax[k] -= fx/mass;
                ay[k] -= fy/mass;
				az[k] -= fz/mass;

				ax[j] += fx/mass;
				ay[j] += fy/mass;
				az[j] += fz/mass;
            }
        }
    }
}
