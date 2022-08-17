#include <windows.h>
#include <stdio.h>
#include <GL/glut.h>
#include <math.h>
#include <time.h>
#define PI 3.1416

/*Delay function */

void delay(float number_of_seconds)
{
    // Converting time into milli_seconds
    float milli_seconds = 1000 * number_of_seconds;

    // Storing start time
    clock_t start_time = clock();

    // looping till required time is not achieved
    while (clock() < start_time + milli_seconds);
}

/* Random number generator function (Mean 0 and variance 1) */

double RNG()
{
	static double U, V;
	static int phase = 0;
	double Z;

	if(phase == 0) {
		U = (rand() + 1.) / (RAND_MAX + 2.);
		V = rand() / (RAND_MAX + 1.);
		Z = sqrt(-2 * log(U)) * sin(2 * PI * V);
	  }
    else
		Z = sqrt(-2 * log(U)) * cos(2 * PI * V);

	phase = 1 - phase;

	return Z;
}

void display(void)
{
/* Simulation of 1D Brownian particles in a heat bath at temperature T and a confining double well potential */
/* Important parameters (Constants) */

float m=pow(10,-14);        /* Mass of each Brownian particle(kg)*/
float R=pow(10,-6);         /* Radius of a Brownian particle (m)*/
float eta=0.001;            /* Viscosity of the fluid bath (at 25 deg C)*/
float T=30;                 /* Temperature of the fluid bath (25 deg C) */
int Np=200;                 /* Number of particles*/
float c=6*PI*eta*R;         /* Damping constant of the fluid bath */
float kb=1.38*pow(10,-23);  /* Boltzmann constant */
float dt=pow(10,-7);        /* Time step value (chosen to ensure convergence)*/
int Nm=10000;               /* Total number of time steps */
float ks1=2*pow(10,-3);     /* Stiffness of the wells */
float ks2=pow(10,-3);       /* Stiffness of springs connecting particles */
float delta=5*pow(10,-9);   /* Distance between right and left wells */
float a=1.5*pow(10,-8);     /* Distance between two left/right wells */
float F0=2*pow(10,-12);     /* Applied force on each particle */


float v[Np];                /* Velocities of all the particles */
float x[Np];                /* Positions of all the particles */
float F[Np];                /* Force on all the particle */
float noise[Np];            /* Noise on all particles */
int N=0;                    /* Initial time step number */


/* Initialize positions of all particles */
/* Half particles in right well (transformed) and other half in left (untransformed) */

for (int i=0;i<Np/2;i+=1)
{
        x[i]=delta+F0/ks1;
      // x[i]=delta;
        v[i]=0;
}

for (int i=Np/2;i<Np;i+=1)
{
       //x[i]=0;
        x[i]=F0/ks1;
        v[i]=0;
}


/* Loop over all time steps */

while(N<Nm)
{

/*  Force due to the confining potential */

for (int i=0;i<Np;i++)
    {
     if (x[i]<delta/2)
     {
         F[i]=-ks1*x[i];
     }

    else
    {
        F[i]=-ks1*(x[i]-delta);
    }

    }

/* Apply force on the particles */

  for (int i=0;i<Np;i++)
    {
          F[i]+=F0;
    }

/* Spring force due to neighbors */

 /* Particle on extreme left */

 F[0]+=ks2*(x[1]-x[0]);

 /* Particle on extreme right */

 F[Np-1]+=ks2*(x[Np-2]-x[Np-1]);


 /* Particles in between */

  for (int i=1;i<Np-1;i++)
    {
        F[i]+= ks2*(x[i+1]-2*x[i]+x[i-1]);
    }

/* Noise on each particle

for (int i=0;i<Np;i++)
    {
        F[i]+=(sqrt(2*kb*T*c)/sqrt(dt))*RNG();

    }

/* Damping (Dissipation) on the particle

  for (int i=0;i<Np;i++)
    {
          F[i]+=-c*v[i];
    }

/* Update velocities and positions with appropriate forces */

for (int k=0;k<Np;k++)
{
v[k]=v[k]+(F[k]*dt/m);
x[k]=x[k]+v[k]*dt;
}


/* Draw all the particles

float th;
glClear(GL_COLOR_BUFFER_BIT);

glColor3f(0.0,0.0,0.0);


  for (int k=0;k<Np;k++)
     {
       glBegin(GL_POLYGON);

        for(int i=0; i<360; i++)
        {
        th = i*PI/180;
        glVertex3f(pow(10,7)*(x[k]+a*k)+pow(10,3.8)*R*cos(th), pow(10,3.8)*(R*sin(th)),0.0);
        }

       glEnd();

     }

/* Draw the potential wells



/* Left wells
glColor3f(0.0, 0.0, 0.0);

for (int j=0;j<Np;j++)
{

 glBegin(GL_LINES);

    for(int i=-500; i<250; i++)
        {
             float xt=i*pow(10,-11)+a*j;
             glVertex3f(pow(10,7)*(xt),pow(10,18.7)*(0.5*ks1*pow(xt-a*j,2)),0.0);
        }
 glEnd();

}

/* Right wells
glColor3f(0.0, 0.0, 0.0);

for (int j=0;j<Np;j++)
{

 glBegin(GL_LINES);

    for(int i=-250; i<500; i++)
        {
             float xt=delta+i*pow(10,-11)+a*j;
             glVertex3f(pow(10,7)*(xt),pow(10,18.7)*(0.5*ks1*pow(xt-delta-a*j,2)),0.0);
        }
 glEnd();

}*/


/* Draw all particles */

float th;
glClear(GL_COLOR_BUFFER_BIT);

glColor3f(0.0,0.0,0.0);


  for (int k=0;k<Np;k++)
     {
       glBegin(GL_POLYGON);

        for(int i=0; i<360; i++)
        {
        th = i*PI/180;
        glVertex3f(0.01*k+pow(10,3.8)*R*cos(th),pow(10,7.5)*x[k]+pow(10,3.8)*(R*sin(th)),0.0);
        }

       glEnd();

     }




       glFlush();



   N+=1;


}

}



 int main(int argc, char **argv)
{

 glutInit(&argc, argv);
 glutInitDisplayMode ( GLUT_SINGLE | GLUT_RGB | GLUT_DEPTH);

 glutInitWindowPosition(0,0);
 glutInitWindowSize(2000,2000);
 glutCreateWindow ("square");


 glClearColor(1.0, 1.0, 1.0, 0.0);           // black background
 glMatrixMode(GL_MODELVIEW);                 // setup viewing projection
 glLoadIdentity();                           // start with identity matrix
 glOrtho(-0.1, 2.2, -0.1, 0.5, 0, 1.0);      // setup a 10x10x2 viewing world

 glutDisplayFunc(display);
 glutMainLoop();

 return 0;
}
