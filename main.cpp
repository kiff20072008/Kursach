#include "Forsythe.h"
using namespace std;




Float K_global;
Float L_global;
#define EPS  1.0e-6
#define M  1.
#define  G  9.81

Float myFunctIntegr1(Float x)
{
    return 0.5819767*(exp(1)-exp(0));
}

void MyFunc(Float t, Float *y, Float *dy)
{   
	dy[0] = y[1];
	dy[1] = -K_global*y[0]/M - G*(1 - cos(y[2])) + (L_global + y[0])*y[3]*y[3];
	dy[2] = y[3];
	dy[3] = -G*sin(y[2])/(L_global + y[0]) - 2*y[1]*y[3]/(L_global + y[0]);
}

Float averageQuadSumm(Float k_local)
{
	
	K_global = k_local;
	Float y0[] = { 0, 0, 0, 4 };//начальные значения
	Float t = 0;
	Float x[] = { 0, 0.303, -0.465, 0.592, -0.409, 0.164, 0.180 };//экспериментальные значения удлинения пружины

	unsigned char work[6 * (4 * sizeof(Float)) + sizeof(rkf_inside)];
	rkf_inside *p;

	rkf myRKF;
	myRKF.f = MyFunc;
	myRKF.Y = y0;
	myRKF.t = t;
	myRKF.tout = 0;
	myRKF.ae = EPS;
	myRKF.re = EPS;
	myRKF.neqn = 4;
	myRKF.flag = 1;
	myRKF.work = work;
	p = (rkf_inside *)myRKF.work;
	Float tout = 0;

	Float summ = 0; //среднеквадратичный критерий

	cout << "t " <<  "             x diff  "  << "              x2"  << "   x diff - x2 " <<endl;
	for (int i = 0; i < 7; i++)
	{
		rkf45(&myRKF);
		cout << setw(3) << myRKF.tout << "   "  << setw(15) << myRKF.Y[0] << "   " << setw(15) << x[i] << "   " << abs(myRKF.Y[0]-x[i]) << endl;
		summ += ((myRKF.Y[0] - x[i]) * (myRKF.Y[0] - x[i]));
		myRKF.tout += 0.4;
	}
	cout << "summ = " << summ << "   K = " << k_local << endl;
	return summ;
}


int main(void)
{
	Float errest, flag;
	int nofun;
	Float res = EPS;
	L_global =
		Quanc8(myFunctIntegr1, 0, 1, EPS, res, &errest, &nofun, &flag);//находим Л;
    Float K = FMin(averageQuadSumm, 36, 46, EPS);
    cout <<setprecision(10) <<"T optimal = " << K << " L =  " << L_global << endl;
	system("pause");
	return 0;
}
