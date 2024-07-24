#include <iostream>
#include <fstream>
#include <sstream>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_multiroots.h>
#include <gsl/gsl_interp.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_sf.h>
#include <list>
#include <vector>
#include <cmath>

#define numberGridPoints 99 
#define addGridPoints 10
//numberGridPoints is #of points when doing ONE grid from -inf to +inf. since we have on both 

using namespace std;


struct parameters{
double peclet1;
double peclet2;
double delta1;
double delta2;
double slope1;	// in entire code, suffix 1 means that the quantity is considered on the right branch.
double slope2;
};

struct integralParameters{
double sigma;
double xObserv;
double yObserv;
gsl_spline* splineG1loc;
gsl_interp_accel* accelG1loc;
gsl_spline* splineG2loc;
gsl_interp_accel* accelG2loc;
gsl_spline* kappaRHP_splineG1;
gsl_spline* kappaLHP_splineG2;
gsl_interp_accel* kappaRHP_accelG1;
gsl_interp_accel* kappaLHP_accelG2;
};
////////////////////////////////// testing variables
//bool shifted(0);
bool testedIntegrandAtxObs(0);
/////////////////////////////////// variables of standard types  ///////////////////////////
const  double PI=4.*atan(1.0);
double alpha(1e5),gammaD(1e5),pG1(1e5),pG2(1e5),deltaG1(1e5),deltaG2(1e5),slopeG1(1e5),slopeG2(1e5),initialDisc(1e5),A_N_L_end(1e5),B_N_L_end(1e5),A_N_R_end(1e5),B_N_R_end(1e5),C_N_L_end(1e5),C_N_R_end(1e5),tailFactor(1e5),angle_two_phi(100.),angle_delta(100.);
int BC_switch(10000);

//////////////////////////////////// gsl specific variables, arrays, lists, vectors, other STL type instances   
size_t iterations(0);
int errorCount(0);
const gsl_interp_type* splineT = gsl_interp_cspline;
double xGrid[numberGridPoints+2*addGridPoints],xGridG1[(numberGridPoints+1)/2+2*addGridPoints],xGridG2[(numberGridPoints+1)/2+2*addGridPoints];
double prolonguedYGridG1[(numberGridPoints+1)/2+2*addGridPoints],prolonguedYGridG2[(numberGridPoints+1)/2+2*addGridPoints]; // the xGrids and yGrids remain untouched - the corresponding value at x=0 ist set manually
vector <double> Discretization(numberGridPoints+2*addGridPoints,10.);
vector <double> A_N_L((numberGridPoints+1)/2+addGridPoints,10.);
vector <double> B_N_L((numberGridPoints+1)/2+addGridPoints,10.);
vector <double> A_N_R((numberGridPoints+1)/2+addGridPoints,10.);
vector <double> B_N_R((numberGridPoints+1)/2+addGridPoints,10.);
vector <double> C_N_L((numberGridPoints+1)/2+addGridPoints,10.);
vector <double> C_N_R((numberGridPoints+1)/2+addGridPoints,10.);
vector <double> locCurvature(numberGridPoints+1,10.);
vector <double> curvature(numberGridPoints+1,10.);
vector <double> HMK_curvature_RHP((numberGridPoints+1)/2,10.);
vector <double> HMK_curvature_LHP((numberGridPoints+1)/2,10.);
vector <double> integrand(numberGridPoints+1,10.);
////////////////////////////////////////// functions which do not return value ////////////////////////////
void setNumericalParameters(double& Length,double& initialDiscretization, double& error_tol, string& gridSwitch, double& p1_o_p2/*, const gsl_vector* x*/);
void initializeGridAndSolution(double& discretization, string& gridSwitch,string& loadValue, gsl_vector* x, double& p1_o_p2);
void printState (gsl_multiroot_fdfsolver* solver/*gsl_multiroot_fsolver* solver*/);
void printChangingParams(gsl_vector* x_solution);
void defineGrid(double& discretization, string& gridSwitch);
void loadSolution();
void initializeSolutionVector(gsl_vector* x,double& p1_o_p2);
void setFileNames(string& namePath,  string& YofXFileName, string& shapeCheckFileName, string& fileName_gnuThermal, int goodConvergeSwitch,string& fileName_IntegralsCheck,string& gridSwitch);
void displayInitialGuess(gsl_vector* x);
void testIvantsovRelationWithLocalSplines(gsl_vector* x);
void defineFittingParabola(const gsl_vector* x);
void displayFittingParabolaAndDisc();
void updateSplines(const gsl_vector* x, gsl_spline* localSplineG1, gsl_spline* localSplineG2, gsl_interp_accel* localAccelG1, gsl_interp_accel* localAccelG2, gsl_spline* kappaRHP_splineG1, gsl_spline* kappaLHP_splineG2, gsl_interp_accel* kappaRHP_accelG1, gsl_interp_accel* kappaLHP_accelG2);
void calc_curv_HMK_RHP(vector<double>& s,vector<double>& x,vector<double>&
y, vector<double>& nx, vector<double>& ny, vector<double>& curv);
void calc_sep_RHP(vector<double>& s,vector<double>& x,vector<double>&
y,vector<double>& xm,vector<double>& ym,vector<double>& sp);
void calc_curv_HMK_LHP(vector<double>& s,vector<double>& x,vector<double>&
y, vector<double>& nx, vector<double>& ny, vector<double>& curv);
void calc_sep_LHP(vector<double>& s,vector<double>& x,vector<double>&
y,vector<double>& xm,vector<double>& ym,vector<double>& sp);
float sgn(float val);


////////////////////////////////////////// functions which do return value ////////////////////////////
int set_equationSystem_Meiron(const gsl_vector* x, void* params, gsl_vector* equationVector);
int set_equationSystem_Meiron_fdf(const gsl_vector* x, void* params, gsl_vector* equationVector, gsl_matrix* jacobian);
int set_jacobian_equationSystem_Meiron(const gsl_vector* x, void* params, gsl_matrix* jacobian);
double exactDiffusionIntegrand_conservationLaw(double x, void* params);
double exactDiffusionIntegrand_locEq(double x, void* params);
 
string IntToString ( int number);
string DoubleToString ( double doubleNumber);

int main(int argc, char** argv)
{
	if ((argv[1]) == NULL)
	{	
		cout << " -------------------------------------------------------- " << endl;
		cout << " fix if load or noload of old result: " << endl;
		cout << " [program name] load " << endl;
		cout << " [program name] noload " << endl; 
		cout << " -------------------------------------------------------- " << endl;
		exit(1);
	}
	string loadSwitch(argv[1]);
	cout << " -------------------------------------------------------- " << endl;
	cout << " call: " << argv[0] << " " << argv[1] << endl;
	cout << " -------------------------------------------------------- " << endl;
///////////////////////////////////////////////////////////////////////////////////////////////
	double discretization(0.),error_tol(0.0001),systemLength(0.),p1_o_p2(0.);
	string gridSwitch;
	int status(0)/*,goodConvergeSwitch(0)*/;
		
	
	gsl_vector* x = gsl_vector_calloc(numberGridPoints);
	setNumericalParameters(systemLength,discretization,error_tol,gridSwitch,p1_o_p2/*,x*/);
	struct parameters meiron_Parameters = {pG1,pG2,deltaG1,deltaG2,slopeG1,slopeG2};
	initializeGridAndSolution(discretization,gridSwitch,loadSwitch,x,p1_o_p2);
	
	defineFittingParabola(x);
	displayFittingParabolaAndDisc();
////////////////////////////////////// UP TO HERE EVERYTHING FINE. //////////////////////////////
//// also tested in setEq -delta-function*1/2 of integrand if proper -Delta/2 comes out for p=0.01,0.1,1.0 - works 

	//exit(1);
////////////////////////////////// asymm Meiron ///////////////////////////////////////////////////////////////
  	const gsl_multiroot_fdfsolver_type* fdfsolverT = gsl_multiroot_fdfsolver_gnewton;/*hybridsj*/
  	gsl_multiroot_fdfsolver* solver1= gsl_multiroot_fdfsolver_alloc(fdfsolverT,numberGridPoints) ;

 	gsl_multiroot_function_fdf model_Meiron = {&set_equationSystem_Meiron,&set_jacobian_equationSystem_Meiron,&set_equationSystem_Meiron_fdf, numberGridPoints, &meiron_Parameters};

	gsl_multiroot_fdfsolver_set(solver1,&model_Meiron,x);
	
	printState(solver1);

//////////////////////////////////////////////////////////////////////////////////////////////////
	do
	{
		iterations++;
		status = gsl_multiroot_fdfsolver_iterate(solver1);
	
		printState(solver1);

		printChangingParams((solver1->x));
		
 

		if (status) {break;}

		status = gsl_multiroot_test_residual(solver1->f, error_tol);
	}
	while (status == GSL_CONTINUE && iterations < 1000);
	
	cout << "status= " << gsl_strerror(status) << endl;

	for (int i=0; i<numberGridPoints;i++)
	{
		cout << "x,root(x),y_Iv: " <<  ( (i<(numberGridPoints-1)/2) ? (xGridG2[addGridPoints+i]) :(xGridG1[addGridPoints+i-(numberGridPoints+1)/2]) ) << " " << gsl_vector_get(solver1->x,i) << " " << ( (i<(numberGridPoints-1)/2) ? -0.5*pow((xGridG2[addGridPoints+i]),2.) :-0.5*pow((xGridG1[addGridPoints+i-(numberGridPoints+1)/2]),2.0)  )  << endl;
		
	}  

	
 	gsl_multiroot_fdfsolver_free(solver1);
	gsl_vector_free(x);

	cout << " hier samma !" << endl;
////////////////////////////////////////////////////////////////////////////////////////////////////

return 0;
}	

///////////////////////////////////////////////////////////////////////
void initializeGridAndSolution(double& discretization, string& gridSwitch, string& loadValue, gsl_vector* x, double& p1_o_p2)
{
	if (loadValue == "load")
	{
		loadSolution();
	}
	else if (loadValue == "noload")
	{
		defineGrid(discretization, gridSwitch);
		initializeSolutionVector(x, p1_o_p2);
	}
	else 
	{
		cout << " -------------------------------------------------------- " << endl;
		cout << " please fix if you want to: " << endl;
		cout << " -------------------------------------------------------- " << endl;
		cout << " load old solutions: enter [program name] load " << endl;
		cout << " start from asymptotical Ivantsov solution: enter [program name] noload " << endl;
		cout << " -------------------------------------------------------- " << endl;
		cout << " now please call program again with described flag values" << endl;
		exit(1);
	}
}
/////////////////////////////////////////////////////////////////////////////////
void loadSolution()
{
	
}
//////////////////////////////////////////////////////////////////////////////////
void defineGrid(double& discretization, string& gridSwitch)
{
	const int halfNumberOfPoints((numberGridPoints+1)/2+addGridPoints);
	list<double> rightGrid;
	 
	int counter(0);
	// Discretization is vector <double> with n+2*add elmnts
/////////////////////////////////////////////////////////////////////////////////////////////////////////	
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////	
	if (gridSwitch=="tanh")
	{
		double x1(1e5),x2(1e5),h1(1e5),h2(1e5),xValue(1e5);
		cout << " (ap+n+1)/2= " << halfNumberOfPoints << "  " << endl;  
		
		rightGrid.push_back(0.);
		for (int xGridIndex=1; xGridIndex < halfNumberOfPoints; ++xGridIndex)
		{
			h1 = 5./((numberGridPoints+addGridPoints-1)/2);
			h2 = (((numberGridPoints+addGridPoints-1)/2)*discretization)/(((numberGridPoints+addGridPoints-1)/2)*(10.*tanh(0.)+11.));
			x1 = xGridIndex*h1;
			x2 = 10.*tanh(x1-5.)+11.;
			xValue = 2.*xGridIndex*h2*x2;
			xValue = 0.0075 + pow(xValue,1.5);
			rightGrid.push_back(xValue);
			
		}
		list<double> entireGrid(rightGrid);
	
		for (list<double>::iterator iter = rightGrid.begin(); iter != rightGrid.end(); ++iter)
		{
			if (iter != rightGrid.begin())
			{
				entireGrid.push_front(-(*iter));
			}
		}
	
		for ( list<double>::iterator iter = entireGrid.begin(); iter != entireGrid.end(); ++iter )
		{
			int xGridIndex(counter);
			if (xGridIndex < (numberGridPoints+2*addGridPoints))
			{
				xGrid[xGridIndex] = *iter; 
			}
			++counter;
		}	 
	
	}
///////////////////////////////////////   Discretization   //////////////////////////////////////////////////////////////////////
	for (int gridIndex = 0; gridIndex < (numberGridPoints+1)/2+addGridPoints; gridIndex++)
	{
		Discretization.at(gridIndex) = xGrid[gridIndex+1]-xGrid[gridIndex];
	}
	Discretization.at((numberGridPoints-1)/2) = Discretization.at((numberGridPoints-1)/2-1);
	for (int gridIndex = (numberGridPoints-1)/2+1+addGridPoints; gridIndex < numberGridPoints+2*addGridPoints; gridIndex++)
	{
		Discretization.at(gridIndex) = xGrid[gridIndex]-xGrid[gridIndex-1];
	}

	for (int i =0; i<(numberGridPoints+1)/2+2*addGridPoints;i++)
	{
		xGridG2[i] = xGrid[i]; 
	}
/////////////////////////////////////////// G1 spline (usually right branch)///////////////////////////////////
	for (int i =0; i<(numberGridPoints+1)/2+2*addGridPoints;i++)
	{
		xGridG1[i] = xGrid[i+(numberGridPoints-1)/2]; 
	}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////	
	std::cout << "gridType: " << gridSwitch << std::endl;

}
////////////////////////////////////////////////////////////////////////
void setNumericalParameters(double& Length, double& initialDiscretization, double& error_tol, string& gridSwitch, double& p1_o_p2)
{

	initialDiscretization = 0.1;
	gridSwitch = "tanh";
	error_tol = 5.*1e-4;
	Length = numberGridPoints*initialDiscretization;
	initialDisc = initialDiscretization; // initialDisc is used as globally declared parameter for normalShift and naming 
	pG1 = 1.0;
	pG2 = 1.0;
	deltaG1 = sqrt(PI*pG1)*exp(pG1)*erfc(sqrt(pG1));
	deltaG2 = sqrt(PI*pG2)*exp(pG2)*erfc(sqrt(pG2));
	p1_o_p2 = pG1/pG2;

	angle_two_phi = 0.5*PI;
	angle_delta = 0.;

	slopeG2 = -1./tan(angle_delta-0.5*angle_two_phi);
	slopeG1 = -1./tan(angle_delta+0.5*angle_two_phi);
	
	

	cout << " slopeG1: " << slopeG1  << " slopeG2: " << slopeG2 << endl;
	//exit(1); 
	
	//slopeG1 = -0.8; // right branch
	//slopeG2 =  0.8; // left branch
	tailFactor = 1.01;

}
/////////////////////////////////////////////////////////////////////////
int set_equationSystem_Meiron(const gsl_vector* x, void* params,gsl_vector* equationVector)
{
 
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	static size_t nr_access(0);
	//cout << " # accesses to setEqSystem " << nr_access << endl;
	++nr_access;
	
// for meiron problem, since slopes are not changed, not necessary yet to set prolongued splines for calculated slope in each step.
	double equation[numberGridPoints];	
for (int i=0;i<numberGridPoints;i++) {equation[i] = 100.;}


	double dydx(0.),yObs(0.),y_0,y_m1,y_p1,kappa(0.),/*kappaHMK(0.),*/d2ydx2(0.),sigma(1e5)/*,delta(1e5)*/;
	//double LocalShift_x, LocalShift_y;
	double xObs(0.);
	double exactIntegral(1e5),logIntegral(100.),error(1e5)/*,xComponent_of_n(10.)*/, /*yComponent_of_n(10.),*//*leftIntegral(10.),rightIntegral(10.),*/leftIntegral_locEq(10.),rightIntegral_locEq(10.),leftIntegral_conservationLaw(10.),rightIntegral_conservationLaw(10.),conservationLawIntegral(10.),locEqIntegral(10.);

	gsl_set_error_handler_off();
	
	gsl_spline* kappaRHP_splineG1 = gsl_spline_alloc(splineT,(numberGridPoints+1)/2);
	gsl_spline* kappaLHP_splineG2 = gsl_spline_alloc(splineT,(numberGridPoints+1)/2);
	gsl_interp_accel* kappaRHP_accelG1 = gsl_interp_accel_alloc();
	gsl_interp_accel* kappaLHP_accelG2 = gsl_interp_accel_alloc();

	gsl_spline* current_splineG1 = gsl_spline_alloc(splineT,(numberGridPoints+1)/2+2*addGridPoints);
	gsl_spline* current_splineG2 = gsl_spline_alloc(splineT,(numberGridPoints+1)/2+2*addGridPoints);
	gsl_interp_accel* current_accelG1 = gsl_interp_accel_alloc();
	gsl_interp_accel* current_accelG2 = gsl_interp_accel_alloc();

	sigma = gsl_vector_get(x,(numberGridPoints+1)/2-1);
	angle_delta = 0.;//gsl_vector_get(x,(numberGridPoints+1)/2);

	slopeG2 = -1./tan(angle_delta-0.5*angle_two_phi);
	slopeG1 = -1./tan(angle_delta+0.5*angle_two_phi);

	updateSplines(x, current_splineG1, current_splineG2, current_accelG1, current_accelG2,kappaRHP_splineG1, kappaLHP_splineG2, kappaRHP_accelG1, kappaLHP_accelG2);
	
	gsl_integration_workspace* w = gsl_integration_workspace_alloc(2000);
	//gsl_function integrateExactDiffKernel;
	
	struct integralParameters currentParams = {sigma,xObs,yObs,current_splineG1,current_accelG1,current_splineG2,current_accelG2,kappaRHP_splineG1, kappaLHP_splineG2, kappaRHP_accelG1, kappaLHP_accelG2};
	
	
	//////////////////////////////////////////////////////////////////////////////////////////////
	y_0 =  gsl_spline_eval(current_splineG2,0.,current_accelG2);
	y_m1 = gsl_spline_eval(current_splineG2,xGridG2[addGridPoints+(numberGridPoints+1)/2-2],current_accelG2);
	y_p1 = gsl_spline_eval(current_splineG1,xGridG1[addGridPoints+1],current_accelG1);
	//y_m1 = y_p1; 
// ////////////////////////////////////////////////// equationVectorSetting ////////////////////////	
	for (int i = 0; i< (numberGridPoints+1)/2-1;i++)
	{
		//testedShiftAtxObs = false;
		size_t n_evaluations;
		xObs = xGridG2[i+addGridPoints];	
		yObs = gsl_spline_eval(current_splineG2,xObs,current_accelG2);
		dydx = gsl_spline_eval_deriv(current_splineG2,xObs,current_accelG2);
		d2ydx2 = gsl_spline_eval_deriv2(current_splineG2,xObs,current_accelG2);
		//kappa = -d2ydx2/pow((1+dydx*dydx),1.5);
		kappa = gsl_spline_eval(kappaLHP_splineG2,xObs,kappaLHP_accelG2);
 /////////////////////////////////////////// arbitrary p integral ////////////////////////////////////////////////
		gsl_function integrateExactDiffKernel;
		integrateExactDiffKernel.function = &exactDiffusionIntegrand_conservationLaw;
		integrateExactDiffKernel.params =  &currentParams;
		currentParams.yObserv = yObs;
		currentParams.xObserv = xObs;
 //////////////////////////////////////////kappa from parabola at outer points ////////////////////////////
		if (i==0)
		{
			//kappa = -d2ydx2/pow((1.+(dydx)*(dydx)),1.5);
			kappa = -2.*C_N_L_end/pow((1.+(B_N_L_end-xObs*C_N_L_end)*(B_N_L_end-xObs*C_N_L_end)),1.5);
		}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		gsl_integration_qng(&integrateExactDiffKernel, tailFactor*xGridG2[addGridPoints], /*xGridG2[addGridPoints+(numberGridPoints+1)/2-2]*/-0.00,0,1e-7,&logIntegral/*exactIntegral*/,&error,&n_evaluations);
		leftIntegral_conservationLaw = logIntegral/*exactIntegral*/;

 		gsl_integration_qng(&integrateExactDiffKernel,/*xGridG1[addGridPoints+1]*/0.00, tailFactor*xGridG1[addGridPoints+(numberGridPoints+1)/2-1],0,1e-7,&logIntegral/*exactIntegral*/,&error,&n_evaluations);
		rightIntegral_conservationLaw = logIntegral/*exactIntegral*/;
		conservationLawIntegral = leftIntegral_conservationLaw + rightIntegral_conservationLaw;
//////////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////arbitrary p integral ////////////////////////////////////////////////////////////////
		integrateExactDiffKernel.function = &exactDiffusionIntegrand_locEq;
		integrateExactDiffKernel.params =  &currentParams;
		currentParams.yObserv = yObs;
		currentParams.xObserv = xObs;
 
/////////////////////////////////////////////////////////////////////////////////////////////////////////
// 		cout << " integr. limit check: " << xGrid[addGridPoints+(numberGridPoints-1)/2-1] << " " << xGrid[addGridPoints+(numberGridPoints-1)/2+1] << endl;
// 		exit(1);

		gsl_integration_qng(&integrateExactDiffKernel, tailFactor*xGridG2[addGridPoints], -0.00/*xGridG2[addGridPoints+(numberGridPoints+1)/2-2]*/,0,1e-7,&logIntegral/*exactIntegral*/,&error,&n_evaluations);
		leftIntegral_locEq = logIntegral/*exactIntegral*/;

		gsl_integration_qng(&integrateExactDiffKernel,0.00/*xGridG1[addGridPoints+1]*/, tailFactor*xGridG1[addGridPoints+(numberGridPoints+1)/2-1],0,1e-7,&logIntegral/*exactIntegral*/,&error,&n_evaluations);
		rightIntegral_locEq = logIntegral/*exactIntegral*/;

		locEqIntegral = leftIntegral_locEq + rightIntegral_locEq;
////////////////////////////////////////////////////////////////

		logIntegral/*exactIntegral*/ = conservationLawIntegral + locEqIntegral;

		equation[i] = 0.5*(deltaG2*(1./*+0.0*exp( -10.*abs(xObs/sigma) )*/) -kappa*sigma)-(pG1/(2.*PI))*logIntegral/*exactIntegral*/; // exact integral
		//equation[i] = -0.5*kappa*sigma - (1./(2.*PI))*logIntegral; // small p integral

		integrand.at(i) = pG1*(1./(2.*PI))*logIntegral/*exactIntegral*/; // integrand is checked for -Delta/2 - behaviour.
		curvature.at(i) = kappa;

 if(!(fabs(equation[i]) < 10000.))
{
	cout << " eq. error @ i=" << i+(numberGridPoints+1)/2 << endl;
	exit(1);
}


	}
// //exit(1);
// 
// 	
 	//equation[(numberGridPoints+1)/2-1] = y_m1-slopeG2*xGridG2[addGridPoints+(numberGridPoints+1)/2-2];
 	equation[(numberGridPoints-1)/2] = y_p1-slopeG1*xGridG1[addGridPoints+1];

//cout << " equation[" << (numberGridPoints+1)/2-1 << "]= "  <<  equation[(numberGridPoints+1)/2-1] << endl;
//cout << " equation[" << (numberGridPoints+1)/2<< "]= " << equation[(numberGridPoints+1)/2] << endl;
// 
	for (int i = 1; i< (numberGridPoints+1)/2;i++)
	{
		size_t n_evaluations;
		xObs = xGridG1[i+addGridPoints];	

		yObs = gsl_spline_eval(current_splineG1,xObs,current_accelG1);
		
		dydx = gsl_spline_eval_deriv(current_splineG1,xObs,current_accelG1);
		d2ydx2 = gsl_spline_eval_deriv2(current_splineG1,xObs,current_accelG1);
		//kappa = -d2ydx2/pow((1+dydx*dydx),1.5);
		kappa = gsl_spline_eval(kappaRHP_splineG1,xObs,kappaRHP_accelG1);
		 /////////////////////////////// arbitray p integral /////////////////////////////
		gsl_function integrateExactDiffKernel;
		integrateExactDiffKernel.function = &exactDiffusionIntegrand_conservationLaw;
		integrateExactDiffKernel.params =  &currentParams;
		currentParams.xObserv = xObs;
		currentParams.yObserv = yObs;
 
///////////////////////////take kappa not from saito spline for outer points /////////////
		if (i==(numberGridPoints+1)/2-1)
		{
			//kappa = -d2ydx2/pow((1.+(dydx)*(dydx)),1.5);
			kappa = -2.*C_N_R_end/pow((1.+(B_N_R_end-xObs*C_N_R_end)*(B_N_R_end-xObs*C_N_R_end)),1.5);
		}
//////////////////////////////////////////////////////////////////////////////////////////
		gsl_integration_qng(&integrateExactDiffKernel, tailFactor*xGridG2[addGridPoints], /*xGridG2[addGridPoints+(numberGridPoints+1)/2-2]*/-0.00,0,1e-7,&logIntegral/*exactIntegral*/,&error,&n_evaluations);
		leftIntegral_conservationLaw = logIntegral/*exactIntegral*/;

 		gsl_integration_qng(&integrateExactDiffKernel,0.00/*xGridG1[addGridPoints+1]*/, tailFactor*xGridG1[addGridPoints+(numberGridPoints+1)/2-1],0,1e-7,&logIntegral/*exactIntegral*/,&error,&n_evaluations);
		rightIntegral_conservationLaw = logIntegral/*exactIntegral*/;
		conservationLawIntegral = leftIntegral_conservationLaw + rightIntegral_conservationLaw;


////////////////////////////////////arbitrary p integral /////////////////////////////////////////////////////////////////////
		integrateExactDiffKernel.function = &exactDiffusionIntegrand_locEq;
		integrateExactDiffKernel.params =  &currentParams;
		currentParams.xObserv = xObs;
		currentParams.yObserv = yObs;
 
//////////////////////////////////////////////////////////////////////////////////////////
		gsl_integration_qng(&integrateExactDiffKernel, tailFactor*xGridG2[addGridPoints], /*xGridG2[addGridPoints+(numberGridPoints+1)/2-2]*/-0.00,0,1e-7,&logIntegral/*exactIntegral*/,&error,&n_evaluations);
		leftIntegral_locEq = logIntegral/*exactIntegral*/;

		gsl_integration_qng(&integrateExactDiffKernel,/*xGridG1[addGridPoints+1]*/0.00, tailFactor*xGridG1[addGridPoints+(numberGridPoints+1)/2-1],0,1e-7,&logIntegral/*exactIntegral*/,&error,&n_evaluations);

		rightIntegral_locEq = logIntegral/*exactIntegral*/;

		locEqIntegral = leftIntegral_locEq+rightIntegral_locEq;

		logIntegral/*exactIntegral*/ = conservationLawIntegral + locEqIntegral;


		integrand.at(i+(numberGridPoints-1)/2) = pG1*(1./(2.*PI))*logIntegral/*exactIntegral*/;
		curvature.at(i+(numberGridPoints-1)/2) = kappa;
		//equation[i+(numberGridPoints-1)/2+2] = (deltaG1-0.5*kappa*sigma)-(pG1/(2.*PI))*logIntegral/*exactIntegral*/; // 1sided exact 
		equation[i+(numberGridPoints-1)/2] = 0.5*(deltaG1*(1./*+0.0*exp( -10.*abs(xObs/sigma) )*/)-kappa*sigma)-(pG1/(2.*PI))*logIntegral;

 	}
 
	for(int i = 0; i< numberGridPoints; i++)
	{
		gsl_vector_set(equationVector,i,equation[i]);
	}
	gsl_integration_workspace_free(w);

	gsl_interp_accel_free(current_accelG1);
	gsl_interp_accel_free(current_accelG2);
	gsl_spline_free(current_splineG1);
	gsl_spline_free(current_splineG2);

	gsl_interp_accel_free(kappaRHP_accelG1);
	gsl_interp_accel_free(kappaLHP_accelG2);
	gsl_spline_free(kappaRHP_splineG1);
	gsl_spline_free(kappaLHP_splineG2);


	return GSL_SUCCESS;

}
////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////
int set_jacobian_equationSystem_Meiron(const gsl_vector* x, void* params, gsl_matrix* jacobian)
{
	double d_x(1e-8);
	gsl_vector* equationVectorAt_x = gsl_vector_alloc(numberGridPoints);
	gsl_vector* equationVectorAt_xdx = gsl_vector_alloc(numberGridPoints);
	gsl_vector* xVector = gsl_vector_alloc(numberGridPoints);
	gsl_vector* xdxVector = gsl_vector_alloc(numberGridPoints);
 	gsl_vector_memcpy(xVector,x);
	gsl_vector_memcpy(xdxVector,x);
	double jacobian_ij(1e5);
	
	set_equationSystem_Meiron(xVector, params, equationVectorAt_x);

	for (int j=0;j<numberGridPoints;j++)
	{
		gsl_vector_memcpy(xdxVector,xVector);
		gsl_vector_set(xdxVector,j,gsl_vector_get(xVector,j)+d_x);
		set_equationSystem_Meiron(xdxVector, params, equationVectorAt_xdx);
		for (int i=0;i<numberGridPoints;i++)
		{
			
			jacobian_ij = ( gsl_vector_get(equationVectorAt_xdx,i) - gsl_vector_get(equationVectorAt_x,i) )/d_x;
			gsl_matrix_set(jacobian, i, j, jacobian_ij);
		}
	}

	gsl_vector_free(equationVectorAt_x);
	gsl_vector_free(equationVectorAt_xdx);
	gsl_vector_free(xVector);
	gsl_vector_free(xdxVector);


	return GSL_SUCCESS;
}
/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////
int set_equationSystem_Meiron_fdf(const gsl_vector* x, void* params, gsl_vector* equationVector, gsl_matrix* jacobian)
{
	set_equationSystem_Meiron(x, params, equationVector);
////////////////////////////////////////////////////////////////////////////////////////////////
	set_jacobian_equationSystem_Meiron(x, params, jacobian);
	
	return GSL_SUCCESS;
}
/////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
void initializeSolutionVector(gsl_vector* x, double& p1_o_p2)
{
	double initialGuess[numberGridPoints+1], sigmaInitialGuess(.25),deltaInitialGuess(0.);
	list<double> solutionG1_G2,shortGrid;
	int counter(0);
/////////////////////////////////////////////////////////////////////////////////////////////////////////	
	for (int xGridIndex = addGridPoints; xGridIndex < numberGridPoints+1+addGridPoints; xGridIndex++)
	{
		//shortGrid.push_back(xGrid[xGridIndex]);	
		if (xGridIndex<(numberGridPoints+1)/2+addGridPoints) 
		{
			shortGrid.push_back(xGridG2[xGridIndex]);
		}
		else if (xGridIndex>=(numberGridPoints+1)/2+addGridPoints)
		{
			shortGrid.push_back(xGridG1[xGridIndex-(numberGridPoints+1)/2]);
		}
	}	 
	
// for (list<double>::iterator iter = shortGrid.begin(); iter != shortGrid.end(); ++iter)
// {
// 	cout << " in initializeSolutionVector 1: shortGrid= " << *iter << endl;	
// }
// exit(1);



	for (list<double>::iterator iter = shortGrid.begin(); iter != shortGrid.end(); ++iter)
	{
		if (*iter <= 0.)
		{
			solutionG1_G2.push_back((*iter)*slopeG2-p1_o_p2*0.5*(*iter)*(*iter));
		}
		else if (*iter >= 0.)
		{
			solutionG1_G2.push_back((*iter)*slopeG1-0.5*(*iter)*(*iter));
		}
		else 
		{
			cout << " error in initializing with Ivantsov solution " << endl;
			exit(1);
		}
	}

 	
	for ( list<double>::iterator iter = solutionG1_G2.begin(); iter != solutionG1_G2.end(); ++iter )
	{
		int xGridIndex(counter);
		if (xGridIndex < (numberGridPoints+1))
		{
			initialGuess[xGridIndex] = *iter; 
		}
		++counter;
		initialGuess[(numberGridPoints+1)/2-1] = deltaInitialGuess;
		initialGuess[(numberGridPoints+1)/2] = sigmaInitialGuess;
	}	 
	 
	//for (int i = 0; i<numberGridPoints+1; i++)
	//{
		 
	//	cout << i << " " << initialGuess[i] << endl;
	//}	
	for (int i = 0; i<=(numberGridPoints-1)/2-1; i++)
	{
		gsl_vector_set(x,i,initialGuess[i]);
	}	
	gsl_vector_set(x,(numberGridPoints-1)/2,initialGuess[(numberGridPoints-1)/2+1]);
	for (int i = (numberGridPoints-1)/2+1; i<=numberGridPoints-1; i++)
	{
		gsl_vector_set(x,i,initialGuess[i+1]);
	}	
	
	for (int i = 0; i<((numberGridPoints-1)/2); i++)
	{
		gsl_vector_set(x,i,initialGuess[i]);
		//cout << i << " L " << xGrid[i+addGridPoints] << " " << gsl_vector_get(x,i) << endl;
	}	
	for (int i = ((numberGridPoints-1)/2); i<numberGridPoints; i++)
	{
		gsl_vector_set(x,i,initialGuess[i+1]);
		//cout << i << " R " << xGrid[i+addGridPoints] <<  " " << gsl_vector_get(x,i) << endl;
	}	
	for (int i = 0; i<numberGridPoints; i++) 
	{
		//cout << i << " test ini x[]: x_i, x[i] = " << xGrid[i+addGridPoints] << " " << gsl_vector_get(x,i) << endl; 
	}
	//exit(1);

}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
void printState (gsl_multiroot_fdfsolver* solver)
{
	double sum_error(0.);
	for (int i=0; i<numberGridPoints+1;i++)	
	{
		cout << " f(i): " << gsl_vector_get(solver->f, i) <<  " @i= " << i << endl;
		if (i<(numberGridPoints+1)/2) 
		{
			cout << " x= " << xGridG2[addGridPoints+i] << endl; 
		}
		else if (i>=(numberGridPoints+1)/2) 
		{
			cout << " x= " << xGridG1[addGridPoints+i-(numberGridPoints+1)/2] << endl; 
		}

		sum_error += fabs(gsl_vector_get(solver->f, i));
	}
	cout << "iterations= " << iterations << " x_tip(sigma)= " << gsl_vector_get(solver->x, (numberGridPoints+1)/2-1) << " x_tip(delta)= " << gsl_vector_get(solver->x, (numberGridPoints+1)/2) << " sum_i |f_i(x)|= " << sum_error << endl;    


}
//////////////////////////////////////////////////////////////////////////////////////////////////////////
void setFileNames(string& pathName, string& YofXFileName,  string& shapeCheckFileName, string& fileName_gnuShape, int goodConvergeSwitch,string& fileName_IntegralsCheck,string& gridSwitch)
{
	
	string al = DoubleToString(alpha);
	string Gd = DoubleToString(gammaD);
	string PG1 = DoubleToString(pG1);
	string PG2 = DoubleToString(pG2);
	string n = IntToString(numberGridPoints);
	string tipS1 = DoubleToString(slopeG1); // Gamma 1 is right interface: x >0 , y<0
	string tipS2 = DoubleToString(slopeG2); // Gamma 2 is left interface: x<0,y<0
	string Ver = "asymm_arbitraryAngles_";
	string BC = IntToString(BC_switch);
	string Disc = DoubleToString(initialDisc);
	string YofX = "YofX_";
	string goodConverge = IntToString(goodConvergeSwitch); // goodConvergeSwitch = 1 means good Convergence
//////////////////////////////////////////////////////////
/////////// change path  /////////////////////
	pathName = "symmResults/";

	string configuration = Ver+"_n_"+n+"_disc_"+Disc+"_grid_"+gridSwitch+"_pG1_"+PG1+"_pG2_"+PG2+"_slopeG1_"+tipS1+"_slopeG2_"+tipS2+"_Gd_"+Gd+"_BC_"+BC+"_goodConverge_"+goodConverge+".dat";

/////////////////////////////////////////////////////////////////////////////////////////////////
	YofXFileName = YofX+configuration;
	shapeCheckFileName = "shapeCheck_"+Ver+"_n_"+n+"_disc_"+Disc+"_Gd_"+Gd+".dat";
	fileName_IntegralsCheck = "integralsCheck_"+Ver+"_n_"+n+"_disc_"+Disc+"_Gd_"+Gd+".dat";
	fileName_gnuShape = "checkShape.gnu";
}
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
string IntToString ( int number )
{
  std::ostringstream oss;
  
  oss<< number;
 
  return oss.str();
}
//////////////////////////////////////////////////////////
string DoubleToString ( double doubleNumber )
{
  std::ostringstream oss;
 
  oss<< doubleNumber;
  
  return oss.str();
}
////////////////////////////////////////////////////////////
void displayInitialGuess(gsl_vector* x)
{

}
//////////////////////////////////////////////////////////
void updateSplines(const gsl_vector* x, gsl_spline* localSplineG1, gsl_spline* localSplineG2, gsl_interp_accel* localAccelG1, gsl_interp_accel* localAccelG2, gsl_spline* kappaRHP_splineG1, gsl_spline* kappaLHP_splineG2, gsl_interp_accel* kappaRHP_accelG1, gsl_interp_accel* kappaLHP_accelG2)
{
static int nrAccess(0);
++nrAccess;

	double interpolatedDataG1[(numberGridPoints+1)/2+2*addGridPoints],interpolatedDataG2[(numberGridPoints+1)/2+2*addGridPoints],interpolatedKappaLHP[(numberGridPoints+1)/2],interpolatedKappaRHP[(numberGridPoints+1)/2];
	
	defineFittingParabola(x);


	for(int i =0; i< (numberGridPoints+1)/2+2*addGridPoints;i++)
	{	
		interpolatedDataG2[i] = prolonguedYGridG2[i];
		interpolatedDataG1[i] = prolonguedYGridG1[i];
	}


	gsl_spline_init(localSplineG2,xGridG2,interpolatedDataG2,(numberGridPoints+1)/2+2*addGridPoints);
	gsl_spline_init(localSplineG1,xGridG1,interpolatedDataG1,(numberGridPoints+1)/2+2*addGridPoints);


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	vector <double> s_vecRHP((numberGridPoints+1)/2,0.);
	vector <double> x_vecRHP((numberGridPoints+1)/2,0.);
	vector <double> y_vecRHP((numberGridPoints+1)/2,0.);
	vector <double> nx_vecRHP((numberGridPoints+1)/2,0.);
	vector <double> ny_vecRHP((numberGridPoints+1)/2,0.);
	vector <double> curv_vecRHP((numberGridPoints+1)/2,0.);
	vector <double> s_vecLHP((numberGridPoints+1)/2,0.);
	vector <double> x_vecLHP((numberGridPoints+1)/2,0.);
	vector <double> y_vecLHP((numberGridPoints+1)/2,0.);
	vector <double> nx_vecLHP((numberGridPoints+1)/2,0.);
	vector <double> ny_vecLHP((numberGridPoints+1)/2,0.);
	vector <double> curv_vecLHP((numberGridPoints+1)/2,0.); // before +1)/2

	for(int i =0; i< (numberGridPoints+1)/2;i++)
	{	
		x_vecRHP.at(i) = xGridG1[addGridPoints+i];
		x_vecLHP.at(i) = xGridG2[(numberGridPoints-1)/2+addGridPoints-i];
		y_vecRHP.at(i) = prolonguedYGridG1[addGridPoints+i];
		y_vecLHP.at(i) = prolonguedYGridG2[(numberGridPoints-1)/2+addGridPoints-i];
	}


	calc_curv_HMK_LHP(s_vecLHP, x_vecLHP, y_vecLHP, nx_vecLHP, ny_vecLHP, curv_vecLHP);
 	calc_curv_HMK_RHP(s_vecRHP, x_vecRHP, y_vecRHP, nx_vecRHP, ny_vecRHP, curv_vecRHP);
 	


	size_t xGridIndex(0);
	for (vector<double>::iterator iterRHP = curv_vecRHP.begin(), iterLHP = curv_vecLHP.end(); iterLHP != curv_vecLHP.begin(); ++iterRHP, --iterLHP) 
	{
 		interpolatedKappaLHP[xGridIndex] = *(iterLHP-1);
		interpolatedKappaRHP[xGridIndex] = *iterRHP;

		++xGridIndex;
	}
	interpolatedKappaLHP[(numberGridPoints-1)/2] = interpolatedKappaLHP[(numberGridPoints-1)/2-1];
	interpolatedKappaRHP[0] = interpolatedKappaLHP[1];	

	double shortXGridG1[(numberGridPoints+1)/2],shortXGridG2[(numberGridPoints+1)/2];
	for (int i=0; i< (numberGridPoints+1)/2;i++)
	{
		shortXGridG1[i] = xGridG1[addGridPoints+i];	
		shortXGridG2[i] = xGridG2[addGridPoints+i];
	}
	gsl_spline_init(kappaLHP_splineG2,shortXGridG2,interpolatedKappaLHP,(numberGridPoints+1)/2);
	gsl_spline_init(kappaRHP_splineG1,shortXGridG1,interpolatedKappaRHP,(numberGridPoints+1)/2);

//////////////////////////////////////////////////////////////////////////////
 
}
///////////////////////////////////////////////////////////////////////////
double exactDiffusionIntegrand_conservationLaw(double x, void* params)
{
	
	double xObs =((struct integralParameters*) params)-> xObserv;
	
	double KernelSource(10.),Kernel(10.),eta(10.),kappaLoc(10.)/*,Kernel_nGradGreensFunc(10.)*/;
	
	double yObs(1e5);
	double y(1e5),locDisc(10e-7);
	double dydx(0.),d2ydx2(0.);
	double y_x_m_locDisc(15.),y_x_p_locDisc(12.);

	if (xObs <= 0.)
	{
		//yObs = gsl_spline_eval(( ((struct integralParameters*) params)-> splineG2loc),xObs,( ((struct integralParameters*) params)-> accelG2loc));
		yObs = ((struct integralParameters*) params)-> yObserv;	

		if (x <= xGrid[addGridPoints])
		{
			y = A_N_L_end+B_N_L_end*x+C_N_L_end*x*x;
			y_x_m_locDisc = A_N_L_end+B_N_L_end*(x-locDisc)+C_N_L_end*(x-locDisc)*(x-locDisc);
			y_x_p_locDisc = A_N_L_end+B_N_L_end*(x+locDisc)+C_N_L_end*(x+locDisc)*(x+locDisc);

			dydx = -(y_x_m_locDisc-y)/locDisc;

			d2ydx2 = (y_x_m_locDisc-2.*y+y_x_p_locDisc)/(locDisc*locDisc);
			kappaLoc = - ( d2ydx2/pow((1.+dydx*dydx),1.5) );  // straightforward kappa, bad close to zero
			//cout << " 8.0.1 @ " << xObs << " " << x << " kappa= " << kappaLoc << endl;
			//exit(1);
		}
		else if ( (x>xGrid[addGridPoints]) && (x <= 0.))
		{
			y = gsl_spline_eval(( ((struct integralParameters*) params)-> splineG2loc),x,( ((struct integralParameters*) params)-> accelG2loc));

			dydx = -( gsl_spline_eval(( ((struct integralParameters*) params)-> splineG2loc),x-locDisc,( ((struct integralParameters*) params)-> accelG2loc)) - gsl_spline_eval(( ((struct integralParameters*) params)-> splineG2loc),x,( ((struct integralParameters*) params)-> accelG2loc)) )/locDisc;  // straightforward kappa, bad close to 0
			
			kappaLoc = gsl_spline_eval(( ((struct integralParameters*) params)-> kappaLHP_splineG2),x,( ((struct integralParameters*) params)-> kappaLHP_accelG2));
// 			cout << " 8.0.1a @ " << xObs << " " << x << " kappa= " << kappaLoc << endl;
// 			exit(1);
		}
		else if ( (x >= 0.) && (x<xGrid[addGridPoints+numberGridPoints-1]) )
		{
			//if (xObs == 0.) {cout << "xObs=0, y= " << yObs << endl;}
			y = gsl_spline_eval(( ((struct integralParameters*) params)-> splineG1loc),x,( ((struct integralParameters*) params)-> accelG1loc));

			dydx = ( gsl_spline_eval(( ((struct integralParameters*) params)-> splineG1loc),x,( ((struct integralParameters*) params)-> accelG1loc)) - gsl_spline_eval(( ((struct integralParameters*) params)-> splineG1loc),x-locDisc,( ((struct integralParameters*) params)-> accelG1loc)) )/locDisc;
			//cout << " 8.0.3 @ " << xObs << " " << x << endl;

			kappaLoc = gsl_spline_eval(( ((struct integralParameters*) params)-> kappaRHP_splineG1),x,( ((struct integralParameters*) params)-> kappaRHP_accelG1));
			//cout << " 8.0@" << xObs << " " << x << endl;
			
		}
		else if ( (x >= xGrid[addGridPoints+numberGridPoints-1]) )
		{
			y = A_N_R_end+B_N_R_end*x+C_N_R_end*x*x;
			y_x_m_locDisc = A_N_R_end+B_N_R_end*(x-locDisc)+C_N_R_end*(x-locDisc)*(x-locDisc);
			y_x_p_locDisc = A_N_R_end+B_N_R_end*(x+locDisc)+C_N_R_end*(x+locDisc)*(x+locDisc);

			dydx = (y_x_p_locDisc-y)/locDisc;

			d2ydx2 = (y_x_p_locDisc-2.*y+y_x_m_locDisc)/(locDisc*locDisc);
			
			kappaLoc = - ( d2ydx2/pow((1.+dydx*dydx),1.5) );  // straightforward kappa, bad close to zero
			//cout << " 8.0.4 @ " << xObs << " " << x << endl;	
			//if (fabs(xObs-11.0385)<=0.01) {cout << "xObs=11.0385, y= " << yObs << endl;}
		}
		else 
		{
			cout << " error in integration in testIvantsov relation b " << endl;
			exit(1);
		}
///////////////////////////////////////////////////// her xObs <= 0

/////////////////////////////////////////////////////////
	}

	else if (xObs >= 0.)
	{
		yObs = ((struct integralParameters*) params)-> yObserv;

		if (x <= xGrid[addGridPoints])
		{
			y = A_N_L_end+B_N_L_end*x+C_N_L_end*x*x;
			y_x_m_locDisc = A_N_L_end+B_N_L_end*(x-locDisc)+C_N_L_end*(x-locDisc)*(x-locDisc);
			y_x_p_locDisc = A_N_L_end+B_N_L_end*(x+locDisc)+C_N_L_end*(x+locDisc)*(x+locDisc);

			dydx = -(y_x_m_locDisc-y)/locDisc;

			d2ydx2 = (y_x_m_locDisc-2.*y+y_x_p_locDisc)/(locDisc*locDisc);
			
			kappaLoc = - ( d2ydx2/pow((1.+dydx*dydx),1.5) );  // straightforward kappa, bad close to zero
		//	cout << " 8.1 1@ " << xObs << " " << x << endl;
		}
		else if ( (x>xGrid[addGridPoints]) && (x <= 0.))
		{
			y = gsl_spline_eval(( ((struct integralParameters*) params)-> splineG2loc),x,( ((struct integralParameters*) params)-> accelG2loc));

			dydx = -( gsl_spline_eval(( ((struct integralParameters*) params)-> splineG2loc),x-locDisc,( ((struct integralParameters*) params)-> accelG2loc)) - gsl_spline_eval(( ((struct integralParameters*) params)-> splineG2loc),x,( ((struct integralParameters*) params)-> accelG2loc)) )/locDisc;  // straightforward kappa, bad close to 0

			kappaLoc = gsl_spline_eval(( ((struct integralParameters*) params)-> kappaLHP_splineG2),x,( ((struct integralParameters*) params)-> kappaLHP_accelG2));
		}
		else if ( (x >= 0.) && (x<xGrid[addGridPoints+numberGridPoints-1]) )
		{
			y = gsl_spline_eval(( ((struct integralParameters*) params)-> splineG1loc),x,( ((struct integralParameters*) params)-> accelG1loc));

			dydx = ( gsl_spline_eval(( ((struct integralParameters*) params)-> splineG1loc),x,( ((struct integralParameters*) params)-> accelG1loc)) - gsl_spline_eval(( ((struct integralParameters*) params)-> splineG1loc),x-locDisc,( ((struct integralParameters*) params)-> accelG1loc)) )/locDisc;
			
			kappaLoc = gsl_spline_eval(( ((struct integralParameters*) params)-> kappaRHP_splineG1),x,( ((struct integralParameters*) params)-> kappaRHP_accelG1));
		}
		else if ( (x >= xGrid[addGridPoints+numberGridPoints-1]) )
		{
			y = A_N_R_end+B_N_R_end*x+C_N_R_end*x*x;
			y_x_m_locDisc = A_N_R_end+B_N_R_end*(x-locDisc)+C_N_R_end*(x-locDisc)*(x-locDisc);
			y_x_p_locDisc = A_N_R_end+B_N_R_end*(x+locDisc)+C_N_R_end*(x+locDisc)*(x+locDisc);

			dydx = (y-y_x_m_locDisc)/locDisc;

			d2ydx2 = (y_x_p_locDisc-2.*y+y_x_m_locDisc)/(locDisc*locDisc);

			kappaLoc = - ( d2ydx2/pow((1.+dydx*dydx),1.5) );  // straightforward kappa, bad close to zero
		}
		else 
		{
			cout << " error in integration in testIvantsov relation c " << endl;
			exit(1);
		}
	}

	else 
	{
		cout << " error in integration in testIvantsov relation a " << endl;
		exit(1);
	}

	eta = sqrt( (x-xObs)*(x-xObs) + (y-yObs)*(y-yObs) );
 
	if ((pG1*eta) >= 0.0000001)
	{
		if ((pG1*eta) >= 10.)
		{
			//KernelSource = sqrt(PI/(2.*pG1*eta))*2.*exp( -pG1*( (yObs-y)+fabs(yObs-y)*(1.+0.5*pow((xObs-x),2.0)/pow((yObs-y),2.0)) )  ); problem: division by 0 if in some iteration step y(x) is symmetric possible. 
			KernelSource = sqrt(PI/(2.*pG1*eta))*2.*exp( -pG1*( (yObs-y)+ eta ));
		}
		else 
		{
			KernelSource = 2.*exp(-(pG1*(yObs-y)))*gsl_sf_bessel_K0(pG1*eta);
		}
		Kernel = KernelSource;
	}

	else 
	{
		Kernel = 0.;
	}

	if (!(fabs(Kernel) <= 10000.))
	{
		cout << " source kernel= " << Kernel << endl;
		cout << " exp(-(pG1*(yObs-y)))= " << exp(-(pG1*(yObs-y))) << " gsl_sf_bessel_K0(pG1*eta)= " << gsl_sf_bessel_K0(pG1*eta) << endl;
		cout << " xObs,x,yObs,y= " << xObs << " " << x << " " << yObs << " " << y << endl;
		cout << " sqrt(PI/(2.*pG1*eta))= " << sqrt(PI/(2.*pG1*eta)) << " exp( -pG1*( (yObs-y)+fabs(yObs-y)*(1.+0.5*pow((xObs-x),2.0)/pow((yObs-y),2.0)) )  )= " << exp( -pG1*( (yObs-y)+fabs(yObs-y)*(1.+0.5*pow((xObs-x),2.0)/pow((yObs-y),2.0)) )  ) <<endl;
		cout << " fabs(yObs-y)*(1.+0.5*pow((xObs-x),2.0)/pow((yObs-y),2.0))= "	<< fabs(yObs-y)*(1.+0.5*pow((xObs-x),2.0)/pow((yObs-y),2.0)) << endl;	
		cout << " eta= " << eta << endl; 
		exit(1);
	}

	return Kernel;
}
////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////
double exactDiffusionIntegrand_locEq(double x, void* params)
{
	
	double sigmaLoc = ((struct integralParameters*) params)-> sigma;
	double xObs =((struct integralParameters*) params)-> xObserv;
	
	double Kernel1sidedPart(10.),Kernel(10.),eta(10.),locEq(10.),kappaLoc(10.);
	
	double yObs(1e5);
	double y(1e5),locDisc(10e-4);
	double dydx(0.),d2ydx2(0.);
	double y_x_m_locDisc(15.),y_x_p_locDisc(12.);

	if (xObs <= 0.)
	{
		//yObs = gsl_spline_eval(( ((struct integralParameters*) params)-> splineG2loc),xObs,( ((struct integralParameters*) params)-> accelG2loc));
		yObs = ((struct integralParameters*) params)-> yObserv;	

		if (x <= xGrid[addGridPoints])
		{
			y = A_N_L_end+B_N_L_end*x+C_N_L_end*x*x;
			y_x_m_locDisc = A_N_L_end+B_N_L_end*(x-locDisc)+C_N_L_end*(x-locDisc)*(x-locDisc);
			y_x_p_locDisc = A_N_L_end+B_N_L_end*(x+locDisc)+C_N_L_end*(x+locDisc)*(x+locDisc);

			dydx = -(y_x_m_locDisc-y)/locDisc;

			d2ydx2 = (y_x_m_locDisc-2.*y+y_x_p_locDisc)/(locDisc*locDisc);
			kappaLoc = - ( d2ydx2/pow((1.+dydx*dydx),1.5) );  // straightforward kappa, bad close to zero
			//cout << " 8.0.1 @ " << xObs << " " << x << " kappa= " << kappaLoc << endl;
			//exit(1);
		}
		else if ( (x>xGrid[addGridPoints]) && (x <= 0.))
		{
			y = gsl_spline_eval(( ((struct integralParameters*) params)-> splineG2loc),x,( ((struct integralParameters*) params)-> accelG2loc));

			dydx = -( gsl_spline_eval(( ((struct integralParameters*) params)-> splineG2loc),x-locDisc,( ((struct integralParameters*) params)-> accelG2loc)) - gsl_spline_eval(( ((struct integralParameters*) params)-> splineG2loc),x,( ((struct integralParameters*) params)-> accelG2loc)) )/locDisc;  // straightforward kappa, bad close to 0
			
			kappaLoc = gsl_spline_eval(( ((struct integralParameters*) params)-> kappaLHP_splineG2),x,( ((struct integralParameters*) params)-> kappaLHP_accelG2));
// 			cout << " 8.0.1a @ " << xObs << " " << x << " kappa= " << kappaLoc << endl;
// 			exit(1);
		}
		else if ( (x >= 0.) && (x<xGrid[addGridPoints+numberGridPoints-1]) )
		{
			//if (xObs == 0.) {cout << "xObs=0, y= " << yObs << endl;}
			y = gsl_spline_eval(( ((struct integralParameters*) params)-> splineG1loc),x,( ((struct integralParameters*) params)-> accelG1loc));

			dydx = ( gsl_spline_eval(( ((struct integralParameters*) params)-> splineG1loc),x,( ((struct integralParameters*) params)-> accelG1loc)) - gsl_spline_eval(( ((struct integralParameters*) params)-> splineG1loc),x-locDisc,( ((struct integralParameters*) params)-> accelG1loc)) )/locDisc;
			//cout << " 8.0.3 @ " << xObs << " " << x << endl;

			kappaLoc = gsl_spline_eval(( ((struct integralParameters*) params)-> kappaRHP_splineG1),x,( ((struct integralParameters*) params)-> kappaRHP_accelG1));
			//cout << " 8.0@" << xObs << " " << x << endl;
			
		}
		else if ( (x >= xGrid[addGridPoints+numberGridPoints-1]) )
		{
			y = A_N_R_end+B_N_R_end*x+C_N_R_end*x*x;
			y_x_m_locDisc = A_N_R_end+B_N_R_end*(x-locDisc)+C_N_R_end*(x-locDisc)*(x-locDisc);
			y_x_p_locDisc = A_N_R_end+B_N_R_end*(x+locDisc)+C_N_R_end*(x+locDisc)*(x+locDisc);

			dydx = (y_x_p_locDisc-y)/locDisc;

			d2ydx2 = (y_x_p_locDisc-2.*y+y_x_m_locDisc)/(locDisc*locDisc);
			
			kappaLoc = - ( d2ydx2/pow((1.+dydx*dydx),1.5) );  // straightforward kappa, bad close to zero
			//cout << " 8.0.4 @ " << xObs << " " << x << endl;	
			//if (fabs(xObs-11.0385)<=0.01) {cout << "xObs=11.0385, y= " << yObs << endl;}
		}
		else 
		{
			cout << " error in integration in testIvantsov relation b " << endl;
			exit(1);
		}

/////////////////////////////////////////////////////////
	}

	else if (xObs >= 0.)
	{
		yObs = ((struct integralParameters*) params)-> yObserv;

		if (x <= xGrid[addGridPoints])
		{
			y = A_N_L_end+B_N_L_end*x+C_N_L_end*x*x;
			y_x_m_locDisc = A_N_L_end+B_N_L_end*(x-locDisc)+C_N_L_end*(x-locDisc)*(x-locDisc);
			y_x_p_locDisc = A_N_L_end+B_N_L_end*(x+locDisc)+C_N_L_end*(x+locDisc)*(x+locDisc);

			dydx = -(y_x_m_locDisc-y)/locDisc;

			d2ydx2 = (y_x_m_locDisc-2.*y+y_x_p_locDisc)/(locDisc*locDisc);
			
			kappaLoc = - ( d2ydx2/pow((1.+dydx*dydx),1.5) );  // straightforward kappa, bad close to zero
		//	cout << " 8.1 1@ " << xObs << " " << x << endl;
		}
		else if ( (x>xGrid[addGridPoints]) && (x <= 0.))
		{
			y = gsl_spline_eval(( ((struct integralParameters*) params)-> splineG2loc),x,( ((struct integralParameters*) params)-> accelG2loc));

			dydx = -( gsl_spline_eval(( ((struct integralParameters*) params)-> splineG2loc),x-locDisc,( ((struct integralParameters*) params)-> accelG2loc)) - gsl_spline_eval(( ((struct integralParameters*) params)-> splineG2loc),x,( ((struct integralParameters*) params)-> accelG2loc)) )/locDisc;  // straightforward kappa, bad close to 0

			kappaLoc = gsl_spline_eval(( ((struct integralParameters*) params)-> kappaLHP_splineG2),x,( ((struct integralParameters*) params)-> kappaLHP_accelG2));
		}
		else if ( (x >= 0.) && (x<xGrid[addGridPoints+numberGridPoints-1]) )
		{
			y = gsl_spline_eval(( ((struct integralParameters*) params)-> splineG1loc),x,( ((struct integralParameters*) params)-> accelG1loc));

			dydx = ( gsl_spline_eval(( ((struct integralParameters*) params)-> splineG1loc),x,( ((struct integralParameters*) params)-> accelG1loc)) - gsl_spline_eval(( ((struct integralParameters*) params)-> splineG1loc),x-locDisc,( ((struct integralParameters*) params)-> accelG1loc)) )/locDisc;
			
			kappaLoc = gsl_spline_eval(( ((struct integralParameters*) params)-> kappaRHP_splineG1),x,( ((struct integralParameters*) params)-> kappaRHP_accelG1));
		}
		else if ( (x >= xGrid[addGridPoints+numberGridPoints-1]) )
		{
			y = A_N_R_end+B_N_R_end*x+C_N_R_end*x*x;
			y_x_m_locDisc = A_N_R_end+B_N_R_end*(x-locDisc)+C_N_R_end*(x-locDisc)*(x-locDisc);
			y_x_p_locDisc = A_N_R_end+B_N_R_end*(x+locDisc)+C_N_R_end*(x+locDisc)*(x+locDisc);

			dydx = (y-y_x_m_locDisc)/locDisc;

			d2ydx2 = (y_x_p_locDisc-2.*y+y_x_m_locDisc)/(locDisc*locDisc);

			kappaLoc = - ( d2ydx2/pow((1.+dydx*dydx),1.5) );  // straightforward kappa, bad close to zero
		}
		else 
		{
			cout << " error in integration in testIvantsov relation c " << endl;
			exit(1);
		}
	}

	else 
	{
		cout << " error in integration in testIvantsov relation a " << endl;
		exit(1);
	}

	eta = sqrt( (x-xObs)*(x-xObs) + (y-yObs)*(y-yObs) );
 	if (x < 0.) 
 	{ 
		locEq = deltaG2-kappaLoc*sigmaLoc;
 	}
 	else if (x > 0.)
 	{
 		locEq = deltaG1-kappaLoc*sigmaLoc;
 	}

	//if (eta != 0.)
	if (((pG1*eta) >= 0.0000001) && (fabs(x)>0.00001))
	{
		if ((pG1*eta) >= 10.)
		{
			//Kernel1sidedPart = locEq*sqrt(PI/(2.*pG1*eta))*exp( -pG1*( (yObs-y)+fabs(yObs-y)*(1.+0.5*pow((xObs-x),2.0)/pow((yObs-y),2.0)) )  )*( ((-dydx*(xObs-x)+(yObs-y))/eta) - 1. );   problem: divison by zero possible if in iteration y(x) symmetric 
			Kernel1sidedPart = locEq*sqrt(PI/(2.*pG1*eta))*exp( -pG1*( (yObs-y)+eta ) )*( ((-dydx*(xObs-x)+(yObs-y))/eta)*(1.+(3./8.)*(1./(pG1*eta))) - 1. );
		}
		else
		{
			Kernel1sidedPart = locEq*exp(-(pG1*(yObs-y)))*( ((-dydx*(xObs-x)+(yObs-y))/eta)*gsl_sf_bessel_K1(pG1*eta) - gsl_sf_bessel_K0(pG1*eta) );
		}
		Kernel = Kernel1sidedPart;
	}

	else 
	{
		Kernel = 0.;
	}

	if (!(fabs(Kernel) <= 10000.)) {cout << " locEq kernel: " << Kernel << endl; exit(1);  }

	return Kernel;
}
 

/////////////////////////////////////////////////////////////////////////////////////////////////////
void defineFittingParabola(const gsl_vector* x)
{
double discl1_L,discl2_L,discm1_L,discm2_L,discu1_L,discu2_L,y_k_L,y_k_plus_L,y_k_minus_L,x_k_minus_L,x_k_L,x_k_plus_L;
double discl1_R,discl2_R,discm1_R,discm2_R,discu1_R,discu2_R,y_k_R,y_k_plus_R,y_k_minus_R,x_k_minus_R,x_k_R,x_k_plus_R;
////////////////////////////// for right branch x>0 //////////////////////
	discl1_R = -Discretization.at(numberGridPoints+addGridPoints-3);
	discl2_R = -Discretization.at(numberGridPoints+addGridPoints-3)-Discretization.at(numberGridPoints+addGridPoints-2);
	discm1_R = Discretization.at(numberGridPoints+addGridPoints-3);
	discm2_R = -Discretization.at(numberGridPoints+addGridPoints-2);
	discu1_R = Discretization.at(numberGridPoints+addGridPoints-2)+Discretization.at(numberGridPoints+addGridPoints-3);
	discu2_R = Discretization.at(numberGridPoints+addGridPoints-2);
	y_k_minus_R = gsl_vector_get(x,numberGridPoints+1-4);
	y_k_R =gsl_vector_get(x,numberGridPoints+1-3);
	y_k_plus_R = gsl_vector_get(x,numberGridPoints+1-2);
	x_k_minus_R = xGrid[numberGridPoints+addGridPoints-4];
	x_k_R = xGrid[numberGridPoints+addGridPoints-3];
	x_k_plus_R = xGrid[numberGridPoints+addGridPoints-2];

	A_N_R.at((numberGridPoints+1)/2-1) = y_k_minus_R*x_k_R*x_k_plus_R/(discl1_R*discl2_R) + y_k_R*x_k_minus_R*x_k_plus_R/(discm1_R*discm2_R) + y_k_plus_R*x_k_minus_R*x_k_R/(discu1_R*discu2_R);
	B_N_R.at((numberGridPoints+1)/2-1) = - y_k_minus_R*(x_k_R+x_k_plus_R) /(discl1_R*discl2_R) -  y_k_R*(x_k_minus_R+x_k_plus_R)/(discm1_R*discm2_R) -  y_k_plus_R*(x_k_minus_R+x_k_R)/(discu1_R*discu2_R);
	C_N_R.at((numberGridPoints+1)/2-1) =  y_k_minus_R/(discl1_R*discl2_R) + y_k_R/(discm1_R*discm2_R) + y_k_plus_R/(discu1_R*discu2_R);

	A_N_R_end = A_N_R.at((numberGridPoints+1)/2-1);
	B_N_R_end = B_N_R.at((numberGridPoints+1)/2-1);
	C_N_R_end = C_N_R.at((numberGridPoints+1)/2-1);

	///////////////////////////  for left branch x<0  /////////////////////////
	discl1_L = -Discretization.at(2+addGridPoints);
	discl2_L = -Discretization.at(2+addGridPoints)-Discretization.at(1+addGridPoints);
	discm1_L = Discretization.at(2+addGridPoints);
	discm2_L = -Discretization.at(1+addGridPoints);
	discu1_L = Discretization.at(1+addGridPoints)+Discretization.at(2+addGridPoints);
	discu2_L = Discretization.at(1+addGridPoints);
	y_k_minus_L = gsl_vector_get(x,3);
	y_k_L =gsl_vector_get(x,2);
	y_k_plus_L = gsl_vector_get(x,1);
	x_k_minus_L = xGrid[3+addGridPoints];
	x_k_L = xGrid[2+addGridPoints];
	x_k_plus_L = xGrid[1+addGridPoints];
	A_N_L.at(0+addGridPoints) = y_k_minus_L*x_k_L*x_k_plus_L/(discl1_L*discl2_L) + y_k_L*x_k_minus_L*x_k_plus_L/(discm1_L*discm2_L) + y_k_plus_L*x_k_minus_L*x_k_L/(discu1_L*discu2_L);
	B_N_L.at(0+addGridPoints) = - y_k_minus_L*(x_k_L+x_k_plus_L) /(discl1_L*discl2_L) -  y_k_L*(x_k_minus_L+x_k_plus_L)/(discm1_L*discm2_L) -  y_k_plus_L*(x_k_minus_L+x_k_L)/(discu1_L*discu2_L);
	C_N_L.at(0+addGridPoints) =  y_k_minus_L/(discl1_L*discl2_L) + y_k_L/(discm1_L*discm2_L) + y_k_plus_L/(discu1_L*discu2_L);

	A_N_L_end = A_N_L.at(0+addGridPoints);
	B_N_L_end = B_N_L.at(0+addGridPoints);
	C_N_L_end = C_N_L.at(0+addGridPoints);

	
	for (int sndGridIndex = 0; sndGridIndex < addGridPoints; ++sndGridIndex)
	{
		prolonguedYGridG2[sndGridIndex] = A_N_L_end+B_N_L_end*xGrid[sndGridIndex]+C_N_L_end*xGrid[sndGridIndex]*xGrid[sndGridIndex];
		prolonguedYGridG2[sndGridIndex+(numberGridPoints+1)/2+addGridPoints] =	slopeG2*xGrid[sndGridIndex+(numberGridPoints+1)/2+addGridPoints];
		// linker branch
	}

	for (int sndGridIndex = 0; sndGridIndex < addGridPoints; ++sndGridIndex)
	{
		prolonguedYGridG1[sndGridIndex] = slopeG1*xGrid[sndGridIndex+(numberGridPoints-1)/2];
		prolonguedYGridG1[sndGridIndex+(numberGridPoints+1)/2+addGridPoints] =	A_N_R_end+B_N_R_end*xGrid[sndGridIndex+numberGridPoints+addGridPoints-1]+C_N_R_end*xGrid[sndGridIndex+numberGridPoints+addGridPoints-1]*xGrid[sndGridIndex+numberGridPoints+addGridPoints-1];
		// rechter branch
	}

	for (int sndGridIndex = 0; sndGridIndex < (numberGridPoints+1)/2; ++sndGridIndex)
	{
		prolonguedYGridG1[sndGridIndex+addGridPoints] = gsl_vector_get(x,sndGridIndex+(numberGridPoints-1)/2);
		//prolonguedYGridG2[sndGridIndex+addGridPoints] =	gsl_vector_get(x,sndGridIndex);
	}
	for (int sndGridIndex = 0; sndGridIndex < (numberGridPoints-1)/2; ++sndGridIndex)
	{
		prolonguedYGridG2[sndGridIndex+addGridPoints] =	gsl_vector_get(x,sndGridIndex);
	}
	
	prolonguedYGridG1[addGridPoints] = 0.;
	prolonguedYGridG2[addGridPoints+(numberGridPoints-1)/2] = 0.;

}
/////////////////////////////////////////////////////////////////////////////////////////////////////
void displayFittingParabolaAndDisc()
{
	cout << " show fitting parabolas " << endl;
	cout << " parameters: " << endl;
	cout << " A_N_L_end,B_N_L_end,C_N_L_end: " << A_N_L_end << " " << B_N_L_end << " " << C_N_L_end << endl;
	cout << " A_N_R_end,B_N_R_end,C_N_R_end: " << A_N_R_end << " " << B_N_R_end << " " << C_N_R_end << endl;
	for (int i = 0; i< (numberGridPoints+1)/2+2*addGridPoints;i++)
	{
		cout << "(x,y_prol)G2 " << xGridG2[i]  << "," << prolonguedYGridG2[i] << endl;
	}
	cout << " ------------------------------------------------------- " << endl;
	for (int i = 0; i< (numberGridPoints+1)/2+2*addGridPoints;i++)
	{
		cout << "(x,y_prol)G1 " << xGridG1[i]  << "," << prolonguedYGridG1[i] << endl;
	}
	cout << " ------------------------------------------------------- " << endl;

}
////////////////////////////////////////////////////////////////
void printChangingParams(gsl_vector* x_solution)
{
	cout << " coefficients of parabolic prolongation " << endl;
	cout << " left branch: ANLend,BNLend,CNLend: " << A_N_L_end << " " << B_N_L_end << " " << C_N_L_end << endl;  
	cout << " right branch: ANRend,BNRend,CNRend: " << A_N_R_end << " " << B_N_R_end << " " << C_N_R_end << endl;
	cout << " coefficients of calculated shape last point of each branch/(-xObs*xObs/2) approx 1? " << endl;
	cout << " left side: y(-xMax)/y_Iv(-xMax)" <<  gsl_vector_get(x_solution,0)/(-0.5*(xGridG2[addGridPoints]*xGridG2[addGridPoints])) << endl;
	cout << " right side: y(xMax)/y_Iv(xMax)" <<  gsl_vector_get(x_solution,numberGridPoints-1)/(-0.5*(xGridG1[(numberGridPoints+1)/2+addGridPoints-1]*xGridG1[(numberGridPoints+1)/2+addGridPoints-1])) << endl; 
	//cout << " xSolution(xMax)= " << gsl_vector_get(x_solution,numberGridPoints) << " xMax= " << xGridG1[(numberGridPoints+1)/2+addGridPoints-1] << endl;
}
////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////
void setNormalVectorComponents (double & xComponent_of_n, double & yComponent_of_n,struct integralParameters* intParams, double & xp)
{
	double corr_dydx(0.);
	if (xp < 0.)
	{
		corr_dydx = gsl_spline_eval_deriv(intParams-> splineG2loc,xp,intParams-> accelG2loc);
	}
	else if (xp > 0.)
	{
		corr_dydx = gsl_spline_eval_deriv(intParams-> splineG1loc,xp,intParams-> accelG1loc);
	}
	double theta2 = atan(corr_dydx);   
	xComponent_of_n = -sin(theta2);
	yComponent_of_n = cos(theta2);

	if (!(fabs(xComponent_of_n)<10000.))
	{
		cout << " xp= " << xp << " yObs  " << gsl_spline_eval(intParams-> splineG2loc,xp,intParams-> accelG2loc) << endl;
		cout << " in set nX, nY: " << xComponent_of_n << " &nY= " << yComponent_of_n << endl;
		cout << " theta2 " << theta2 << endl;
		exit(1);
	}
}
////////////////////////////////////////////////////////////////////////////
void calc_sep_RHP(vector<double>& s,vector<double>& x,vector<double>&
y,vector<double>& xm,vector<double>& ym,vector<double>& sp)
{
	double dx, dy, ds,x1,y1;
	int n=x.size();
	x1=-x[1];
	y1=y[1];
	ds=s[0];
	dx=x[0]-x1;
	dy=y[0]-y1;
	sp[0]=sqrt(dx*dx+dy*dy);
	xm[0]=dx/sp[0];
	ym[0]=dy/sp[0];


  for(int m=1;m<n;++m)
     {
      ds=s[m]-s[m-1];
      dx=x[m]-x[m-1];
      dy=y[m]-y[m-1];
      sp[m]=sqrt(dx*dx+dy*dy);
      xm[m]=dx/sp[m];
      ym[m]=dy/sp[m];
  
     }
}
////////////////////////////////////////////////////////////////////////////
void calc_sep_LHP(vector<double>& s,vector<double>& x,vector<double>&
y,vector<double>& xm,vector<double>& ym,vector<double>& sp)
{
	double dx, dy, ds,x1,y1;
	int n=x.size();
	x1=-x[1];
	y1=y[1];
	ds=s[0];
	dx=x[0]-x1;
	dy=y[0]-y1;
	sp[0]=sqrt(dx*dx+dy*dy);
	xm[0]=dx/sp[0];
	ym[0]=dy/sp[0];


  for(int m=1;m<n;++m)
     {
      ds=s[m]-s[m-1];
      dx=x[m]-x[m-1];
      dy=y[m]-y[m-1];
      sp[m]=sqrt(dx*dx+dy*dy);
      xm[m]=dx/sp[m];
      ym[m]=dy/sp[m];
  
     }

}

////////////////////////////////////////////////////////////////////////////////////
void calc_curv_HMK_RHP(vector<double>& s,vector<double>& x,vector<double>&
y, vector<double>& nx, vector<double>& ny, vector<double>& curv)
{
	
  int n=x.size();
  vector<double> xm(n),ym(n),sp(n);
  
  double bunsen;  
  double bunshi;
  double bbx, bby;
  double bunbo;
	
  calc_sep_RHP(s,x,y,xm,ym,sp);
	
  for(int m=0;m<n-1;m++)
     {
      bunsen=xm[m]*xm[m+1]+ym[m]*ym[m+1];
      bunshi=xm[m]*ym[m+1]-ym[m]*xm[m+1];
      bbx=-sp[m]*ym[m+1]-sp[m+1]*ym[m];
      bby=sp[m]*xm[m+1]+ sp[m+1]*xm[m];

      bunbo = sqrt(bbx*bbx + bby*bby);
 
      curv[m]=-2.*bunshi/bunbo;
 
     if ((bunbo==0))
     {
	cout << " RHP/m=" << m << " bbx,bby= " << bbx << "," << bby << endl; 
	cout << " bunsen, bunshi= " << bunsen << " " << bunshi << endl;
	cout << " m=0: sp,ym,xm " << sp[m] << " " << ym[m] << " " << xm[m] << endl;
	cout << " m+1: sp,ym,xm " << sp[m+1] << " " << ym[m+1] << " " << xm[m+1] << endl; 
     	cout << "nan curv. Warning!\n";
 
        exit(1);
     }
     }
 
 
}
////////////////////////////////////////////////////////////////////////////////////
void calc_curv_HMK_LHP(vector<double>& s,vector<double>& x,vector<double>&
y, vector<double>& nx, vector<double>& ny, vector<double>& curv)
{
	
  int n=x.size();
  vector<double> xm(n),ym(n),sp(n);
  
  double bunsen;  
  double bunshi;
  double bbx, bby;
  double bunbo;
	
  calc_sep_LHP(s,x,y,xm,ym,sp);
	
  for(int m=0;m<n-1;m++)
     {
      bunsen=xm[m]*xm[m+1]+ym[m]*ym[m+1];
      bunshi=xm[m]*ym[m+1]-ym[m]*xm[m+1];
      bbx=-sp[m]*ym[m+1]-sp[m+1]*ym[m];
      bby=sp[m]*xm[m+1]+ sp[m+1]*xm[m];

      bunbo = sqrt(bbx*bbx + bby*bby);
 
      curv[m]=2.*bunshi/bunbo;
 
       if (bunbo==0)
        {
	cout << " LHP/m=" << m << " bbx,bby= " << bbx << "," << bby << endl; 
	cout << " bunsen, bunshi= " << bunsen << " " << bunshi << endl;
	cout << " m=0: sp,ym,xm " << sp[m] << " " << ym[m] << " " << xm[m] << endl;
	cout << " m+1: sp,ym,xm " << sp[m+1] << " " << ym[m+1] << " " << xm[m+1] << endl; 
        cout << "nan curv. Warning!\n";
 
        exit(1);
        }
     }
 
 
}
/////////////////////////////////////////////////////////////////////////////
float sgn(float val)
{
double outValue(10.);
    if (val<0)
    {
        outValue = -1.;
    }
    else if (val>0)
    {
        outValue= 1.;
    }
    else if (val==0)
    {
        outValue= 1.;
    }

	return outValue;
} 







