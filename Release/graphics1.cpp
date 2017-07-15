/*	Project:	IR-polarisation interaction modeling
	Author:		Josh Bowden
				Queensland University of Technology
	Date:		09/11/2006
	
	Program Description:
	
	Animates a distribution of 'theoretical' collagen fibers and calculates the % interaction
	possible with plane polarised IR light beam using a cos^2(theta) interaction function. Theta
	is the angle between the plane of IR-polarised light and the 'dipole'. Dipole is assumed to be
	between the two atomic centres.
*/


/* headers to libraries to be looked through */
#include<stdio.h>
#include<math.h>
#include <time.h>
#include<windows.h>
#include<stdlib.h> 
//#include<GL/gl.h>
//#include<GL/glu.h>
#include<GL/glut.h>


/* some constants that may or may not be needed */
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#define BACKBONELENGTH 30 
#define NUMBEROFFIBRES 18
#define AMIDEISPACING  5 



/* declare global window handler return values from glutCreateWindow() */
GLint triangleWin ;   

/* structures used */
struct vect2D {
	float x  ;
	float y  ;
	float z  ;
} ;


struct vect2Dflag {
	float x  ;
	float y  ;
	float z  ;
	int numConnect ;
} ;
// this allows 65536 vertecies to be indexed
// il, i2 and i3 are the index into the array vertex data for a triangle (3 points)
struct vertexIndex{
	short int i1 ;
	short int i2 ;
	short int i3 ;
	short int padding ;
};
struct lineEqutn{
	struct vect2D p1 ;
	struct vect2D p2 ;
} ;

struct spherCoords {
	float dist	;
	float theta		;  // in radians
	float phi		;  // in radians
} ;

struct polarCoords {
	float dist	;
	float alpha		;  // angle between x axis and line to point in radians
	float beta		;  // angle between y axis and line to point in radians
	float gamma		;  // angle between z axis and line to pointin radians
} ;

struct collagenIRInteraction{
	double carbonyl ;
	double amideII  ;
} ;


struct spherCoords	spherEyeCoords ; // These are pre-calculated from glEyePos coordinates with getSphericalCoords()
struct spherCoords	getSphericalCoords(struct vect2D p2) ;
struct vect2D		getXYZFromSpherCoords(struct spherCoords in) ;

/* prototypes for GLUT callback functions */
void RenderWindow(void) ;						/* for rendering callback function */
void ChangeViewport(GLsizei w, GLsizei h) ;		/* for resize callback function */
void WindowVisible(int vis) ;					/* glutVisibilityFunc() */
void setEyePos() ;								/* glutIdleFunc() */
void keyFunction(unsigned char key, int x, int y) ;
void mouseButtonFunction(int button, int state, int xPos, int yPos) ;
void mouseMoveFunction(int xPos, int yPos) ;

/* setup OpenGL initial state */
void initialiseGLState(void) ;

void DrawXYZScales(float length) ;


/* auxillary functions used by above*/
struct lineEqutn calcParametricLine(float parametricT, struct lineEqutn twoVect2D) ;
int				 inArray(struct vect2D *array, struct vect2D p1, int arrySize) ;
void			 sortByClosest(struct vect2D *linePoints, int arraySize) ;


double			 dToR(double deg) ;
void			 changeColor(int colorNum) ;
void			 sortByClosest(struct vect2D *linePoints, int arraySize) ;


/* global data storage */
struct vect2D glEyePos = { 60.0, 60.0, 60.0 };  // array storing current eye position for gluLookAt()
struct vect2D glEyeFocus = { 0.0, 0.0, 0.0 };
struct vect2D glUpDir = { 0.0, 1.0, 0.0 };

float twist = 0.0 ;
float elevation = 45.0 ;
float azimuth = 45.0 ;

char axisRotate = 'z' ;
float modelview[16] ;
// angle1 is the angle the fiber long axis makes with the z axis
// Is a global varible because keyboard used to change value
// Use 'a' or 's' to add or subtract 15 degrees
float angle1_as = 54.7 ;
float angle1_jk = 0 ;
float EVectAngle = 0; 
static long newSeed = 554433341 ;

int widthScreen, heightScreen ;

bool relToPerfectRadDist = true ;
int   CosOrSin = 0 ;


float	gaussStdev = 0 ;
int		fullOrHalf = 3 ;  // full radial distribution or half distribution (changed by 't' key)
float   minLength = 0.5 ;  // the minimum length of fibres along the split line

float x1[4] = {0,0,1,0} ;
float x2[4] = {0.5,0,1,0} ;

static int oldPos[] = {0,0} ;

enum  buttonState {buttonUp, buttonDown}  ;

struct buttons {
	int leftButton   ;
	int rightButton  ;
} leftRightButton ;

 


/*****************************************************************
* Function name	: main				     			  
* Description	: main entry point to program
*   
* Return type	: int									     
*										      
******************************************************************/
int main(int argc, char* argv[])
{

	for (int i = 1; i == 16 ; i++)
		modelview[i] = 0 ;

	spherEyeCoords = getSphericalCoords(glEyePos) ;

	glutInit(&argc, argv) ;
	glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGBA  | GLUT_DEPTH) ;

	printf("\nBrief instructions:") ;
	printf("\nkeys a and s change fibre angle WRT x (red) axis") ;
	printf("\nkeys e and w change IR beam plane around z (blue) axis") ;
	printf("\nkeys d and f changes standard deviation of fibres gaussian distribution") ;
	printf("\nkeys g and h re-randomises fiber distribution") ;
	printf("\nx,y,z and r changes viewing position to along designated axis") ;
	printf("\nt key changes from full radial distribution to half radial disrtibution") ;
	printf("\nc Rotates IR plane at each increment of fibre (0 to 90) ") ;
	printf("\n/ and ; keys change length of fibres in Sin/Cos distribution ") ;
	printf("\nv key toggles output from relative to perfect distribution to not relative ") ;
	printf("\nOutput:") ;
	printf("\n<IR angle> <Fibre angle WRT x axis> <% Amide I interaction>  <% Amide II interaction>  <% Amide I/II ratio> ") ;

	triangleWin = glutCreateWindow("collagen bond distribution") ;
	glutDisplayFunc(RenderWindow);
	glutKeyboardFunc(keyFunction);
	glutMouseFunc(mouseButtonFunction);
	glutMotionFunc(mouseMoveFunction);
	glutReshapeFunc(ChangeViewport) ;
	glutVisibilityFunc(WindowVisible);
	glutPositionWindow(400,50) ;  // set initial positon of viewport window (from top left of monitor)
	glutReshapeWindow(800,800);

	initialiseGLState() ;
	

	glutMainLoop();

	return 0;
}

/*****************************************************************
* Function name	: WindowVisible				     			  
* Description	: Calculates a new XYZ eye position for gluLookAt() using 
				  polar functions    
* Return type	: void - sets global variable glEyePos 										     
* Argument      : void										      
******************************************************************/
void WindowVisible(int vis)
{
  if (vis == GLUT_VISIBLE)
  {
//    glutIdleFunc(setEyePos);
//	glutIdleFunc(NULL) ;
	glutPostRedisplay();
  }
  else
    glutIdleFunc(NULL);
}


/*****************************************************************
* Function name	: setEyePos				     			  
* Description	: Calculates a new XYZ eye position for gluLookAt() using 
				  polar functions    
* Return type	: void - sets global variable glEyePos 										     
* Argument      : void										      
******************************************************************/
void setEyePos()
{
/*	static float angle = 0.0f ;
	int static i = 0 ;	


	glMatrixMode(GL_MODELVIEW);


	Sleep(2000) ;


	glutPostRedisplay(); */
}




/******************************************************************
* Function name	: transformPoint				     					  
* Description	: Multiplies a 3D point by a transormation matrix
*				  				          	
* Return type	: void									      	
* Arguments     : float * mat    - pointer to 16 float values
*				: float * point  - pointer to 3 float values
* Return value	: float * ret
*					  
******************************************************************/

void transformPoint(float * mat, float * p1, float * ret)
{

	ret[0] = mat[0]    * p1[0] +  mat[4]  * p1[1] +  mat[8]  * p1[2] +  mat[12]  * p1[3] ;
	ret[1] = mat[1]    * p1[0] +  mat[5]  * p1[1] +  mat[9]  * p1[2] +  mat[13]  * p1[3];
	ret[2] = mat[2]    * p1[0] +  mat[6]  * p1[1] +  mat[10] * p1[2] +  mat[14]  * p1[3];
	ret[3] = mat[3]    * p1[0] +  mat[7]  * p1[1] +  mat[11] * p1[2] +  mat[15]  * p1[3]; 

/*	ret[0] = mat[0]    * p1[0] +  mat[1]  * p1[1] +  mat[2]  * p1[2] +  mat[3]  * p1[3] ;
	ret[1] = mat[4]    * p1[0] +  mat[5]  * p1[1] +  mat[6]  * p1[2] +  mat[7]  * p1[3];
	ret[2] = mat[8]    * p1[0] +  mat[9]  * p1[1] +  mat[10] * p1[2] +  mat[11]  * p1[3];
	ret[3] = mat[12]   * p1[0] +  mat[13] * p1[1] +  mat[14] * p1[2] +  mat[15]  * p1[3]; */

}

/******************************************************************
* Function name	: printModelViewMat				     					  
* Description	: Multiplies a 3D point by a transormation matrix
*				  				          	
* Return type	: void									      	
* Arguments     : float * mat    - pointer to 16 float values
*				: float * point  - pointer to 3 float values
* Return value	: float * ret
*					  
******************************************************************/

void printModelViewMat(float * mat)
{

	printf( "\n\n %5.3f %5.3f %5.3f %5.3f", mat[0],  mat[4],  mat[8],   mat[12] ) ;
	printf( "\n %5.3f %5.3f %5.3f %5.3f",   mat[1],  mat[5],  mat[9],   mat[13] ) ;
	printf( "\n %5.3f %5.3f %5.3f %5.3f",   mat[2],  mat[6],  mat[10],  mat[14] ) ;
	printf( "\n %5.3f %5.3f %5.3f %5.3f",   mat[3],  mat[7],  mat[11],   mat[15] ) ;


}


/******************************************************************
* Function name	: getCosine				     					  
* Description	: calculates cosine of angle between two lines - may be the dot product
*				  				          										      	
* Arguments     :
*				: float * p1 & p2   	2 points defining lines starting from the origin
* Return value	: double
*					  
******************************************************************/
double  getCosine(float * p1, float * p2 )
{
	double ret ;
	double a, b, c, a1, a2 ;

		
	a = p1[0] * p2[0] ;
	b = p1[1] * p2[1] ;
	c = p1[2] * p2[2] ;

	a1 = sqrt(  p1[0]*p1[0] + p1[1]*p1[1] + p1[2]*p1[2]   ) ;
	a2 = sqrt(  p2[0]*p2[0] + p2[1]*p2[1] + p2[2]*p2[2]   ) ;

	ret = (a + b + c) / (a1 * a2) ;

	
	return ( ret ) ;
	
}





/******************************************************************
* Function name	: drawProjection				     					  
* Description	: projects (line between 2 points) onto 2D plane and draws it using OpenGL
*				  				          										      	
* Arguments     : float * mat    		Pointer to 16 float values - the current transformation matrix
*										returned from glGetFloatv(GL_MODELVIEW_MATRIX, float *) ;
*				: float * p1 & p2   	2 points defining the line representing the dipole in original coordiates
*				: float * planeEqutn 	The equation of the plane to be projected on. 
* Return value	: double
*					  
******************************************************************/
double  GetDipoleCosine(float * mat, float * p1, float * p2, float * polariserAngle )
{
	double cosTheta = 0 ;
	float rp1[4], rp2[4] ;
	
	transformPoint(mat,p1,rp1) ;
	transformPoint(mat,p2,rp2) ;


	rp1[0] = rp1[0] - rp2[0] ;
	rp1[1] = rp1[1] - rp2[1] ;
	rp1[2] = rp1[2] - rp2[2] ;
	rp1[3] = rp1[3] - rp2[3] ;



	cosTheta = getCosine( rp1 , polariserAngle) ;

//	printf("\n length  = %7.5f, original length  = %7.5f", length, origLength ) ;
	
	return ( cosTheta ) ;
	
}



/******************************************************************
* Function name	: drawCollagenBonds				     					  
* Description	: draws a stick of length 'direction.dist' with radiating bonds every 'spacing' (representing collagen C=O bonds) 
*                 with the backbone at direction, direction.rho, direction.theta (this is collagen backbone)
*				  				          	
* Return type	: double	 :  The net 'interaction' 									      	
* Arguments     : 
				: int numBonds		designates if 4 or 8 bands are radiating from single point
*					  
******************************************************************/
collagenIRInteraction  drawCollagenBonds(float spacing, float bondLength, struct spherCoords direction, int numBonds, float  IREvectAngle ) 
{
	
	float posz ;
	posz = direction.dist ;
    float cos1, cos2 ;
	float p1[4], p2[4], irEvector[3] ;
	double  totalThetaSqrd = 0;
	double  tempCosTheta = 0 ;
	collagenIRInteraction retVal ;
	retVal.amideII  = 0;
	retVal.carbonyl = 0 ;

	cos1 = cos( M_PI / 4 ) ;
	cos1 = cos( M_PI / 8 ) ;

	// direction of IR electric vextor in x,y plane (normal is in z)
	// This is vector to which angle of bond is determined
	irEvector[0] = cos(dToR(IREvectAngle))  ; irEvector[1] = sin(dToR(IREvectAngle)) ; irEvector[2] = 0 ;  

	

	p1[3] = 0 ;
	p2[3] = 0 ;
	for (int i = 0 ; posz >= 0  ; i++)
	{

	glBegin( GL_LINES ) ; 
	  //  glColor3f(1.0f, 0.0f, 0.0f) ;		
	    // radiating C=O bonds	
	    p1[0]=0;p1[1]=0;p1[2]=posz; // this is constant for each i

		p2[0]=bondLength;p2[1]=0;p2[2]=posz; // ps[0] and p2[1] change
    	glVertex3fv(p1) ; glVertex3fv(p2) ;
		tempCosTheta= GetDipoleCosine(modelview,p1,p2,irEvector) ;
		retVal.carbonyl += tempCosTheta * tempCosTheta ;

		if ( numBonds = 8) {
		p2[0]=bondLength*cos1;p2[1]=bondLength*cos1;p2[2]=posz; // ps[0] and p2[1] change
    	glVertex3fv(p1) ; glVertex3fv(p2) ;
		tempCosTheta= GetDipoleCosine(modelview,p1,p2,irEvector) ;
		retVal.carbonyl += tempCosTheta * tempCosTheta ;}

		p2[0]=0;p2[1]=bondLength;p2[2]=posz; // ps[0] and p2[1] change
    	glVertex3fv(p1) ; glVertex3fv(p2) ;
		tempCosTheta= GetDipoleCosine(modelview,p1,p2,irEvector) ;
		retVal.carbonyl += tempCosTheta * tempCosTheta ;

		if ( numBonds = 8) {
		p2[0]=-bondLength*cos1;p2[1]=bondLength*cos1;p2[2]=posz; // ps[0] and p2[1] change
    	glVertex3fv(p1) ; glVertex3fv(p2) ;
		tempCosTheta= GetDipoleCosine(modelview,p1,p2,irEvector) ;
		retVal.carbonyl += tempCosTheta * tempCosTheta ;}

		p2[0]=-bondLength;p2[1]=0;p2[2]=posz; // ps[0] and p2[1] change
    	glVertex3fv(p1) ; glVertex3fv(p2) ;
		tempCosTheta= GetDipoleCosine(modelview,p1,p2,irEvector) ;
		retVal.carbonyl += tempCosTheta * tempCosTheta ;

		if ( numBonds = 8) {
		p2[0]=-bondLength*cos1;p2[1]=-bondLength*cos1;p2[2]=posz; // ps[0] and p2[1] change
    	glVertex3fv(p1) ; glVertex3fv(p2) ;
		tempCosTheta= GetDipoleCosine(modelview,p1,p2,irEvector) ;
		retVal.carbonyl += tempCosTheta * tempCosTheta ;}

		p2[0]=0;p2[1]=-bondLength;p2[2]=posz; // ps[0] and p2[1] change
    	glVertex3fv(p1) ; glVertex3fv(p2) ;
		tempCosTheta= GetDipoleCosine(modelview,p1,p2,irEvector) ;
		retVal.carbonyl += tempCosTheta * tempCosTheta ;

		if ( numBonds = 8) {
		p2[0]=bondLength*cos1;p2[1]=-bondLength*cos1;p2[2]=posz; // ps[0] and p2[1] change
    	glVertex3fv(p1) ; glVertex3fv(p2) ;
		tempCosTheta= GetDipoleCosine(modelview,p1,p2,irEvector) ;
		retVal.carbonyl += tempCosTheta * tempCosTheta ; }
		// end of radiating C=O bonds		
	glEnd() ;

	posz = posz - spacing  ;
	
	}

	// backbone line - (Amide II)
	float distTemp = direction.dist ;
	if (distTemp != 0.0) 
	{
	   glBegin( GL_LINES ) ; 
	    p1[0]=0;p1[1]=0;p1[0]=0; 
		p2[0]=0;p2[1]=0;p2[2]=distTemp;
		glVertex3fv(p1) ; glVertex3fv(p2) ;
		tempCosTheta= GetDipoleCosine(modelview,p1,p2,irEvector) ;
		if (relToPerfectRadDist)
		{
			retVal.amideII = tempCosTheta * tempCosTheta * (direction.dist / BACKBONELENGTH) ;
		}
		else
			retVal.amideII = tempCosTheta * tempCosTheta   ;
			
	   glEnd() ;
    }
	else
	{
		retVal.amideII = 0.0 ;
	}
	

//	printf("\n single chain length  = %7.5f", totalLength ) ;

	return ( retVal ) ;
}



/*****************************************************************
* Function name	:	ran0to1		// From 'Numerical Recipies in C'		     			  
* Description	:   returns a random number between 0..1
					Minimal random number generator of Park and Miller. Returns a uniform random deviate
					between 0.0 and 1.0. Set or reset idum to any integer value (except the unlikely value MASK)
					to initialize the sequence; idum must not be altered between calls for successive deviates in
					a sequence.					
* Return type	:	float  - the random value										     
* Argument      :	void										      
******************************************************************/
#define IA 16807
#define IM 2147483647
#define AM (1.0/IM)
#define IQ 127773
#define IR 2836
#define MASK 123459876

float ran0to1(long *idum)
// idum should be defined as 'static long' in calling function i.e. static long seedVar = 5543341 ;

{
	long k;
	float ans;
	
	*idum ^= MASK;					// XORing with MASK allows use of zero and other simple bit patterns for idum.
	
	k=(*idum)/IQ; 
	
	*idum=IA*(*idum-k*IQ)-IR*k;		// Compute idum=(IA*idum) % IM without overflows by Schrage’s method
	
	if (*idum < 0) 
		*idum += IM;		// .

	ans=AM*(*idum);					// Convert idum to a floating result.
	
	*idum ^= MASK;					// Unmask before return.
	
	return ans;
}





/*****************************************************************
* Function name	:	ranGauss				     			  
* Description	:   Returns a normally distributed deviate with zero mean and unit variance, using ran1(idum)
*					as the source of uniform deviates. i.e. a gaussian curve with 0 mean and stddev of 1
*					
* Return type	:	float  - the value										     
* Argument      :	void										      
******************************************************************/
float ranGauss(long *idum)
// idum should be defined as 'static long' in calling function i.e. static long seedVar = 5543341 ;
{
//	float ran1(long *idum);
	static int iset=0;
	static float gset;
	float fac,rsq,v1,v2;

	// Reinitialize.
	if (*idum < 0) iset=0; 

	// We don’t have an extra deviate handy, so pick two uniform numbers in the square extending
	// from -1 to +1 in each direction, see if they are in the unit circle, and if they are not, try again.
	if (iset == 0) 
	{ 
	
		do {
		v1=2.0*ran0to1(idum)-1.0; 
		v2=2.0*ran0to1(idum)-1.0; 
		rsq=v1*v1+v2*v2; 
		
		} while (rsq >= 1.0 || rsq == 0.0); 

		fac=sqrt(-2.0*log(rsq)/rsq);
	
		// Now make the Box-Muller transformation to get two normal deviates. Return one and save the other for next time.

		gset=v1*fac;

		//  Set flag
		iset=1 ;

		return v2*fac;
	} 
	else   // We have an extra deviate handy, so unset the flag, and return it.
	{
		iset=0;  
		return gset; 
	}




/*	// this works
	randNum1 =  ( (float) rand() / (RAND_MAX+1) )   ;
	randNum2 =  ( (float) rand() / (RAND_MAX+1) )   ;
	randNum1 = sqrt( (-1 * 2 * log ( randNum1 )) / randNum1  ) * (cos (2*M_PI*randNum2) );

	return ( randNum1 ) ; */

}





/*****************************************************************
* Function name	:	SymetricFibers				     			  
* Description	:   Draws radialy symetric array of fibers each with
*					  radially symmetric bonds at 90 deg to fibre. 
*					Calculates an 'interaction factor' that is between 0-1 indicating how much a polarised IR beam 
*					interacts with 
* Return type	:	void - sets global variable glEyePos 										     
* Argument      :	void										      
******************************************************************/

collagenIRInteraction SymetricFibers(float angleWRTZ1, float angleWRTZ2, long seedIn) 
{
	struct spherCoords colBackbone ;
	colBackbone.dist =  BACKBONELENGTH ;
	float initialDist = colBackbone.dist ;
	float spacing =  AMIDEISPACING  ; // was 8
	int numBonds = 8 ;  // number radiating around the length
	int numFibers= NUMBEROFFIBRES ; // 18 ;
	static long seed = 554433341 ;
	float randomGauss1, randomGauss2 ;
	float eVectPos = 0, eVectLength = 0 ;
	float totalLength = 0.0 ;

	float minLength = 0.5 ;

	seed = seedIn ;

	double sectionCarbonylThetaSqrd = 0 ;
	double sectionAmideIIThetaSqrd = 0 ;
	collagenIRInteraction retColInteractn2 ;
	retColInteractn2.amideII  = 0;
	retColInteractn2.carbonyl = 0 ;

 
    float angle2 = 0 ;  // this is angle around the z axis (as compared with: angleWRTZ, the angle with the z axis)
	float stepAngle = 2* angle2 / (numFibers -1) ;
	
	
	DrawXYZScales(40.0) ;


	// draw cos curve over plane of IR beam
	glBegin( GL_LINES ) ;
	for (int i = 0 ; i < 40 ; i++ )
	{
	// draw the E vector
    	glVertex3f(cos(dToR(EVectAngle))*eVectLength,sin(dToR(EVectAngle))*eVectLength,eVectPos) ; glVertex3f(0,0,eVectPos) ;
		eVectPos +=  1.0 ;
		eVectLength = cos(eVectPos / 4) * 10 ;
	}
	glEnd() ;

	int numFiberPositions = 0 ;  // this is the number of 'nodes' where carbonyl bonds radiate from
	for ( i = 0 ; i < numFibers ; i++ )
	{
		glColor3f(1.0f, 0.0f, 0.0f) ;	

		glPushMatrix( ) ;
			randomGauss1 = ranGauss(&seed) * gaussStdev ;
			glRotatef(90,0.0,1.0,0.0) ;
			glRotatef(angle2,0.0,0.0,1.0) ;  // this is number that goes around z axis
			glRotatef(angleWRTZ1+randomGauss1,0.0,1.0,0.0) ;  // this is angle wrt z axis (splay angle)  
		  	
		  // repeat identical transform but without the gluLookAt
		  glPushMatrix( ) ; 
		    glLoadIdentity();  // start from fresh matrix
			glRotatef(90,0.0,1.0,0.0) ;
			glRotatef(angle2,0.0,0.0,1.0) ;  // this is number that goes around z axis
			glRotatef(angleWRTZ1+randomGauss1,0.0,1.0,0.0) ;  // this is angle wrt z axis (splay angle)   
			glGetFloatv(GL_MODELVIEW_MATRIX, modelview) ;
		  glPopMatrix() ;

		  
		    retColInteractn2  = drawCollagenBonds( spacing , 5, colBackbone, numBonds, EVectAngle ) ;
		    sectionCarbonylThetaSqrd += retColInteractn2.carbonyl ;
		    sectionAmideIIThetaSqrd  += (retColInteractn2.amideII  * (colBackbone.dist / initialDist)) ;
	    glPopMatrix() ;
	   	 

		angle2 = angle2 + (360 / numFibers) ;
		if ( angle2 > 360) angle2 = 0 ;

		totalLength = totalLength + colBackbone.dist ;

	    numFiberPositions += (int) (colBackbone.dist / spacing ) ;
		numFiberPositions++ ;

	}


	int totalBonds ;	
	// Amide I calculation	
	if (relToPerfectRadDist)
	{
		totalBonds =  ((int) (initialDist / spacing) * numBonds ) * numFibers ; 
		sectionCarbonylThetaSqrd  = sectionCarbonylThetaSqrd / totalBonds ;
	}
	else
		sectionCarbonylThetaSqrd = (sectionCarbonylThetaSqrd / ( numBonds * numFiberPositions )) ;  // used to be * numFibres in denominator

	
	// Amide II calculation	
	totalLength = totalLength / initialDist ;
	sectionAmideIIThetaSqrd  = sectionAmideIIThetaSqrd / totalLength ; //  sectionAmideIIThetaSqrd / numFibers  ; //* initialDist) / totalLength) ;
	
	float ratio ;
	ratio = sectionCarbonylThetaSqrd / sectionAmideIIThetaSqrd ;


	if (leftRightButton.leftButton == buttonUp  && leftRightButton.rightButton == buttonUp && angleWRTZ2 != 0.0 ) 
	{printf("%7.5f %7.5f %7.5f\n", sectionCarbonylThetaSqrd , sectionAmideIIThetaSqrd , ratio) ;}

	retColInteractn2.carbonyl = sectionCarbonylThetaSqrd ;
	retColInteractn2.amideII  = sectionAmideIIThetaSqrd  ;

	return (retColInteractn2) ;


}  // end of SymetricSplitLineFibers



/*****************************************************************
* Function name	:	HalfSymetricFibers				     			  
* Description	:   Draws radialy symetric array of fibers each with
*					  radially symmetric bonds at 90 deg to fibre.   
*					Calculates an 'interaction factor' that is between 0-1 indicating how much a polarised IR beam 
*					interacts with 
* Return type	:	void - sets global variable glEyePos 										     
* Argument      :	void										      
******************************************************************/
collagenIRInteraction HalfSymetricFibers(float angleWRTZ1, float angleWRTZ2, long seedIn) 
{
	struct spherCoords colBackbone ;
	colBackbone.dist =  BACKBONELENGTH ;
	float spacing =  AMIDEISPACING  ;
	int numBonds = 8 ;
	int numFibers= NUMBEROFFIBRES ;
	static long seed = 554433341 ;
	float randomGauss1, randomGauss2 ;
	float eVectPos = 0, eVectLength = 0 ;

	seed = seedIn ;

	double sectionCarbonylThetaSqrd = 0 ;
	double sectionAmideIIThetaSqrd = 0 ;
	collagenIRInteraction retColInteractn ;
	retColInteractn.amideII  = 0;
	retColInteractn.carbonyl = 0 ;

 
    float angle2 = 0 ;  // this is angle around the z axis (as compared with: angleWRTZ, the angle with the z axis)
	float stepAngle = 180 / (numFibers -1) ;
	
	
	DrawXYZScales(40.0) ;


	// draw cos curve over plane of IR beam
	glBegin( GL_LINES ) ;
	for (int i = 0 ; i < 40 ; i++ )
	{
	// draw the E vector
    	glVertex3f(cos(dToR(EVectAngle))*eVectLength,sin(dToR(EVectAngle))*eVectLength,eVectPos) ; glVertex3f(0,0,eVectPos) ;
		eVectPos +=  1.0 ;
		eVectLength = cos(eVectPos / 4) * 10 ;
	}
	glEnd() ;

//	if (leftRightButton.leftButton == buttonUp  && leftRightButton.rightButton == buttonUp && angleWRTZ2 != 0.0 ) 
//	printf("\nE vect: %7.5f  %7.5f  : ", EVectAngle, angleWRTZ1 ) ;
	for ( i = 0 ; i < numFibers ; i++ )
	{
		glColor3f(1.0f, 0.0f, 0.0f) ;

		glPushMatrix( ) ;
			randomGauss1 = ranGauss(&seed) * gaussStdev ;  ;
		//	if (ran0to1(&seed) > 0.5) randomGauss1 = randomGauss1 * -1 ;

			glRotatef(90,0.0,1.0,0.0) ;
			glRotatef(angle2,0.0,0.0,1.0) ;  // this is number that goes around z axis
			// angle1 is angle fiber long axis makes with z axis
			glRotatef(angleWRTZ1+randomGauss1,0.0,1.0,0.0) ;  // this is angle wrt z axis (splay angle)  
		  	
		  // repeat identical transform but without the gluLookAt
		  glPushMatrix( ) ; 
		    glLoadIdentity();  // start from fresh matrix
			glRotatef(90,0.0,1.0,0.0) ;
			glRotatef(angle2,0.0,0.0,1.0) ;  // this is number that goes around z axis
			// angle1 is angle fiber long axis makes with z axis
			glRotatef(angleWRTZ1+randomGauss1,0.0,1.0,0.0) ;  // this is angle wrt z axis (splay angle)   
		
			glGetFloatv(GL_MODELVIEW_MATRIX, modelview) ;
		//	printModelViewMat(modelview) ;
		  glPopMatrix() ;

		  
		    retColInteractn  = drawCollagenBonds( spacing , 5, colBackbone, numBonds, EVectAngle ) ;
		    sectionCarbonylThetaSqrd += retColInteractn.carbonyl ;
		    sectionAmideIIThetaSqrd  += retColInteractn.amideII  ;
	    glPopMatrix() ;
	   	 

		angle2 = angle2 + stepAngle ;
		if ( angle2 > 180) angle2 = 0 ;
		
	}

    int numFiberPositions = (int) (colBackbone.dist / spacing ) ;
    numFiberPositions++ ;



	float ratio ;
	sectionCarbonylThetaSqrd = (sectionCarbonylThetaSqrd / ( numBonds * numFiberPositions * numFibers)) ;
	sectionAmideIIThetaSqrd  = sectionAmideIIThetaSqrd / numFibers ;
	ratio = sectionCarbonylThetaSqrd / sectionAmideIIThetaSqrd ;

	if (leftRightButton.leftButton == buttonUp  && leftRightButton.rightButton == buttonUp && angleWRTZ2 != 0.0 ) 
	{printf("%7.5f %7.5f %7.5f\n", sectionCarbonylThetaSqrd , sectionAmideIIThetaSqrd , ratio) ;}



	retColInteractn.carbonyl = sectionCarbonylThetaSqrd ;
	retColInteractn.amideII  = sectionAmideIIThetaSqrd  ;

	return (retColInteractn) ;


}

/* For radial distribution:
		glRotatef(90,0.0,1.0,0.0) ;
		glRotatef(angle2,0.0,0.0,1.0) ;  // this is number that goes around z axis
		// angle1 is angle fiber long axis makes with z axis
		glRotatef(angleWRTZ1,0.0,1.0,0.0) ;  // this is angle wrt z axis (splay angle)  

  		angle2 = angle2 + (360 / numFibers) ;
		if ( angle2 > 360) angle2 = 0 ;
*/

/*
			randomGauss = ranGauss(&seed) * angleWRTZ ;
			glRotatef(90,0.0,1.0,0.0) ;
			glRotatef(angle2 ,0.0,0.0,1.0) ;  // this is number that goes around z axis
			// angle1 is angle fiber long axis makes with z axis
			glRotatef(randomGauss,0.0,1.0,0.0) ;  // this is angle wrt z axis (splay angle)  
			angle2 = angle2 + (7200/ numFibers) ;
		if ( angle2 > 360) angle2 = 0 ;
*/

/*
			// that produces completley random fibre array with range = 0 and amide I/II of 1.5
			randomGauss1 = ran0to1(&seed) * angleWRTZ1 ;
			randomGauss2 = ran0to1(&seed) * angleWRTZ2 ;

			glRotatef(90,0.0,0.0,1.0) ;
			glRotatef(randomGauss1 ,0.0,0.0,1.0) ;  // this is number that goes around z axis
			// angle1 is angle fiber long axis makes with z axis
			glRotatef(randomGauss2,0.0,1.0,0.0) ;  // this is angle wrt z axis (splay angle) 
*/


/*****************************************************************
* Function name	:	NonSymetricSplitLineFibers				     			  
* Description	:   Draws mirror plane symetric array of fibers each with
*					radially symmetric bonds at 90 deg to fibre.   
*					Calculates an 'interaction factor' that is between 0-1 indicating how much a polarised IR beam 
*					interacts with 
* Return type	:	void - sets global variable glEyePos 										     
* Argument      :	void										      
******************************************************************/
collagenIRInteraction NonSymetricSplitLineFibers(float angleWRTZ1, float angleWRTZ2, long seedIn) 
{
	struct spherCoords colBackbone ;
	colBackbone.dist =  BACKBONELENGTH ;
	float initialDist = colBackbone.dist ;
	float spacing =  AMIDEISPACING  ; // was 8
	int numBonds = 8 ;  // number radiating around the length
	int numFibers= NUMBEROFFIBRES ; // 18 ;
	static long seed = 554433341 ;
	float randomGauss1, randomGauss2 ;
	float eVectPos = 0, eVectLength = 0 ;
	float totalLength = 0.0 ;

	

	seed = seedIn ;

	double sectionCarbonylThetaSqrd = 0 ;
	double sectionAmideIIThetaSqrd = 0 ;
	collagenIRInteraction retColInteractn2 ;
	retColInteractn2.amideII  = 0;
	retColInteractn2.carbonyl = 0 ;

 
//   float angle2 = 0 ;  // this is angle around the z axis (as compared with: angleWRTZ, the angle with the z axis)
//	float stepAngle = 2* angle2 / (numFibers -1) ;
	float angle2 = 0 ;  // this is angle around the z axis (as compared with: angleWRTZ, the angle with the z axis)
	float stepAngle = 360 / (numFibers -1) ;
	
	
	DrawXYZScales(40.0) ;


	// draw cos curve over plane of IR beam
	glBegin( GL_LINES ) ;
	for (int i = 0 ; i < 40 ; i++ )
	{
	// draw the E vector
    	glVertex3f(cos(dToR(EVectAngle))*eVectLength,sin(dToR(EVectAngle))*eVectLength,eVectPos) ; glVertex3f(0,0,eVectPos) ;
		eVectPos +=  1.0 ;
		eVectLength = cos(eVectPos / 4) * 10 ;
	}
	glEnd() ;

	int numFiberPositions = 0 ;  // this is the number of 'nodes' where carbonyl bonds radiate from
	for ( i = 0 ; i < numFibers ; i++ )
	{
		glColor3f(1.0f, 0.0f, 0.0f) ;
		
		// this shortens the length of fibers in one direction
		if (CosOrSin == 1)
			colBackbone.dist = cos((angle2 / 360.0) * M_PI) ;  // change this from sin to cos to change direction if split line by 90 degrees
		else
			colBackbone.dist = sin((angle2 / 360.0) * M_PI) ;  // change this from sin to cos to change direction if split line by 90 degrees	 

		colBackbone.dist = colBackbone.dist * colBackbone.dist ;
		colBackbone.dist = minLength +  (minLength * colBackbone.dist ) ;	// minLength is variable by ''' and '/' keys
		if (	colBackbone.dist > 1 ) 	colBackbone.dist = 1.0 ;
		colBackbone.dist = initialDist * colBackbone.dist ;
		

	/*	if (angle2 > 180) 
			colBackbone.dist = 25 ;  */

		glPushMatrix( ) ;
			randomGauss1 = ranGauss(&seed) * gaussStdev ;
			glRotatef(90,1.0,0.0,0.0) ;
			glRotatef(90,0.0,1.0,0.0) ;
			glRotatef(angle2,0.0,0.0,1.0) ;  // this is number that goes around z axis
			glRotatef(angleWRTZ1+randomGauss1,0.0,1.0,0.0) ;  // this is angle wrt z axis (splay angle)  
		  	
		  // repeat identical transform but without the gluLookAt
		  glPushMatrix( ) ; 
		    glLoadIdentity();  // start from fresh matrix
			glRotatef(90,1.0,0.0,0.0) ;
			glRotatef(90,0.0,1.0,0.0) ;
			glRotatef(angle2,0.0,0.0,1.0) ;  // this is number that goes around z axis
			glRotatef(angleWRTZ1+randomGauss1,0.0,1.0,0.0) ;  // this is angle wrt z axis (splay angle)   
			glGetFloatv(GL_MODELVIEW_MATRIX, modelview) ;
		  glPopMatrix() ;

		  
		    retColInteractn2  = drawCollagenBonds( spacing , 5, colBackbone, numBonds, EVectAngle ) ;
		    sectionCarbonylThetaSqrd += retColInteractn2.carbonyl ;
		    sectionAmideIIThetaSqrd  += (retColInteractn2.amideII  * (colBackbone.dist / initialDist)) ;
	    glPopMatrix() ;
	   	 

	//	angle2 = angle2 + (360 / numFibers) ;
	//	if ( angle2 > 360) angle2 = 0 ;

		totalLength = totalLength + colBackbone.dist ;

	    numFiberPositions += (int) (colBackbone.dist / spacing ) ;
		numFiberPositions++ ;

		angle2 = angle2 + stepAngle ;
		if ( angle2 > 360) angle2 = 0 ;

	}
	
		
	// Amide I calculation	
	if (relToPerfectRadDist)
	{
		int totalBonds ;
		totalBonds =  ((int) (initialDist / spacing) * numBonds ) * numFibers ; 
		sectionCarbonylThetaSqrd  = sectionCarbonylThetaSqrd / totalBonds ;
	}
	else
		sectionCarbonylThetaSqrd = (sectionCarbonylThetaSqrd / ( numBonds * numFiberPositions )) ;  // used to be * numFibres in denominator

	// Amide II calculation	
	totalLength = totalLength / initialDist ;
	sectionAmideIIThetaSqrd  = sectionAmideIIThetaSqrd / totalLength ; //  sectionAmideIIThetaSqrd / numFibers  ; //* initialDist) / totalLength) ;
	
	float ratio ;
	ratio = sectionCarbonylThetaSqrd / sectionAmideIIThetaSqrd ;

	if (leftRightButton.leftButton == buttonUp  && leftRightButton.rightButton == buttonUp && angleWRTZ2 != 0.0 ) 
	{printf("%7.5f %7.5f %7.5f\n", sectionCarbonylThetaSqrd , sectionAmideIIThetaSqrd , ratio) ;}

	retColInteractn2.carbonyl = sectionCarbonylThetaSqrd ;
	retColInteractn2.amideII  = sectionAmideIIThetaSqrd  ;

	return (retColInteractn2) ;


}  // end of SymetricSplitLineFibers

/*****************************************************************
* Function name	:	SymetricRandomSplitLineFibers				     			  
* Description	:   Draws mirror plane symetric array of fibers each with
*					radially symmetric bonds at 90 deg to fibre.   
*					Calculates an 'interaction factor' that is between 0-1 indicating how much a polarised IR beam 
*					interacts with 
* Return type	:	void - sets global variable glEyePos 										     
* Argument      :	void										      
******************************************************************/
collagenIRInteraction SymetricRandomSplitLineFibers(float angleWRTZ1, float angleWRTZ2, long seedIn) 
{
	struct spherCoords colBackbone ;
	colBackbone.dist =  BACKBONELENGTH ;
	float initialDist = colBackbone.dist ;
	float spacing =  AMIDEISPACING ; // was 8
	int numBonds = 8 ;  // number radiating around the length
	int numFibers= NUMBEROFFIBRES ; // 18 ;
	static long seed = 554433341 ;
	float randomGauss1, randomGauss2 ;
	float eVectPos = 0, eVectLength = 0 ;
	float totalLength = 0.0 ;

	

	seed = seedIn ;

	double sectionCarbonylThetaSqrd = 0 ;
	double sectionAmideIIThetaSqrd = 0 ;
	collagenIRInteraction retColInteractn2 ;
	retColInteractn2.amideII  = 0;
	retColInteractn2.carbonyl = 0 ;

 
    float angle2 = 0 ;  // this is angle around the z axis (as compared with: angleWRTZ, the angle with the z axis)
	float stepAngle = 2* angle2 / (numFibers -1) ;
	
	
	DrawXYZScales(40.0) ;
	int i ;

	// draw cos curve over plane of IR beam
	glBegin( GL_LINES ) ;
	for ( i = 0 ; i < 40 ; i++ )
	{
	// draw the E vector
    	glVertex3f(cos(dToR(EVectAngle))*eVectLength,sin(dToR(EVectAngle))*eVectLength,eVectPos) ; glVertex3f(0,0,eVectPos) ;
		eVectPos +=  1.0 ;
		eVectLength = cos(eVectPos / 4) * 10 ;
	}
	glEnd() ;

	int numFiberPositions = 0 ;  // this is the number of 'nodes' where carbonyl bonds radiate from
	for ( i = 0 ; i < numFibers ; i++ )
	{
		glColor3f(1.0f, 0.0f, 0.0f) ;
		
		// this shortens the length of fibers in one direction
		if (CosOrSin == 1)
			colBackbone.dist = cos((angle2 / 180.0) * M_PI) ;  // change this from sin to cos to change direction if split line by 90 degrees
		else
			colBackbone.dist = sin((angle2 / 180.0) * M_PI) ;  // change this from sin to cos to change direction if split line by 90 degrees

		colBackbone.dist = colBackbone.dist * colBackbone.dist ;
		colBackbone.dist = minLength +  (minLength * colBackbone.dist ) ;	// minLength is variable by ''' and '/' keys
		if (	colBackbone.dist > 1 ) 	colBackbone.dist = 1.0 ;
		colBackbone.dist = initialDist * colBackbone.dist ;



		glPushMatrix( ) ;
			randomGauss1 = ranGauss(&seed) * gaussStdev ;
			randomGauss2 = ranGauss(&seed) * gaussStdev ;
			glRotatef(angleWRTZ2,1.0,0.0,0.0) ;	
			glRotatef(90,0.0,1.0,0.0) ;
	//		glRotatef(45,1.0,0.0,0.0) ;
			glRotatef(angle2+randomGauss2,0.0,0.0,1.0) ;  // this is number that goes around z axis
			glRotatef(angleWRTZ1+randomGauss1,0.0,1.0,0.0) ;  // this is angle wrt z axis (splay angle)  
		  	
		  // repeat identical transform but without the gluLookAt
		  glPushMatrix( ) ; 
		    glLoadIdentity();  // start from fresh matrix
			glRotatef(angleWRTZ2,1.0,0.0,0.0) ;
			glRotatef(90,0.0,1.0,0.0) ;
		//	glRotatef(45,1.0,0.0,0.0) ;
			glRotatef(angle2+randomGauss2,0.0,0.0,1.0) ;  // this is number that goes around z axis
			glRotatef(angleWRTZ1+randomGauss1,0.0,1.0,0.0) ;  // this is angle wrt z axis (splay angle)  
			// get the rotation matrix to use later
			glGetFloatv(GL_MODELVIEW_MATRIX, modelview) ;
		  glPopMatrix() ;

		  
		    retColInteractn2  = drawCollagenBonds( spacing , 5, colBackbone, numBonds, EVectAngle ) ;
		    sectionCarbonylThetaSqrd += retColInteractn2.carbonyl ;
		    sectionAmideIIThetaSqrd  += (retColInteractn2.amideII  * (colBackbone.dist / initialDist)) ;
	    glPopMatrix() ;
	   	 


		angle2 = angle2 + (360 / numFibers) ;
		if ( angle2 > 360) angle2 = 0 ;

		totalLength = totalLength + colBackbone.dist ;

	    numFiberPositions += (int) (colBackbone.dist / spacing ) ;
		numFiberPositions++ ;

	}

    
	// Amide I calculation	
	if (relToPerfectRadDist)
	{
		int totalBonds ;
		totalBonds =  ((int) (initialDist / spacing) * numBonds ) * numFibers ; 
		sectionCarbonylThetaSqrd  = sectionCarbonylThetaSqrd / totalBonds ;
	}
	else
		sectionCarbonylThetaSqrd = (sectionCarbonylThetaSqrd / ( numBonds * numFiberPositions )) ;  // used to be * numFibres in denominator

	// Amide II calculation	
	totalLength = totalLength / initialDist ;
	sectionAmideIIThetaSqrd  = sectionAmideIIThetaSqrd / totalLength ; //  sectionAmideIIThetaSqrd / numFibers  ; //* initialDist) / totalLength) ;
	
	float ratio ;
	ratio = sectionCarbonylThetaSqrd / sectionAmideIIThetaSqrd ;


//	if (leftRightButton.leftButton == buttonUp  && leftRightButton.rightButton == buttonUp && angleWRTZ2 != 0.0 ) 
//	{printf("%7.5f %7.5f %7.5f\n", sectionCarbonylThetaSqrd , sectionAmideIIThetaSqrd , ratio) ;}

	retColInteractn2.carbonyl = sectionCarbonylThetaSqrd ;
	retColInteractn2.amideII  = sectionAmideIIThetaSqrd  ;

	return (retColInteractn2) ;


}  // end of SymetricSplitLineFibers



/*****************************************************************
* Function name	: RenderWindow				     			  
* Description	: callback function for glutDisplayFunc().
				  it is used for windows resizing events   
* Return type	: void - sets global variable glEyePos 										     
* Argument      : void										      
******************************************************************/
void RenderWindow(void) 
{

	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT) ; // this does the clearing to the above clear color
	glLoadIdentity();
//	gluLookAt(glEyePos.x, glEyePos.y, glEyePos.z, glEyeFocus.x, glEyeFocus.y, glEyeFocus.z, glUpDir.x, glUpDir.y, glUpDir.z);


	glTranslated(0.0,0.0,-120) ;
	glRotated(-twist,0.0,0.0,1.0) ;
	glRotated(-elevation,1.0,0.0,0.0) ;
	glRotated(azimuth,0.0,0.0,1.0) ;


	if ( fullOrHalf == 1 ) 
		SymetricFibers(angle1_as, 1.0, newSeed) ;
	else if ( fullOrHalf == 2 ) 
		HalfSymetricFibers(angle1_as, 1.0, newSeed) ;
	else if ( fullOrHalf == 3 ) 
		NonSymetricSplitLineFibers(angle1_as, 1.0, newSeed) ;
	else if ( fullOrHalf == 4 ) 
	    SymetricRandomSplitLineFibers(angle1_as, 1.0, newSeed) ; // DrawXYZScales(40.0) ; 


	glutSwapBuffers() ;		// copies back buffer to front
}







float getMinVal(int AmideIorII, float graphvals[2][36])   // when press 'c' key
{
  float currentMin = 1.0 ;
  
  for (int i = 0; i < 36; i++)
  {
		if (graphvals[AmideIorII-1][i] < currentMin) 
			currentMin = graphvals[AmideIorII-1][i] ;

  }

  return currentMin ;

}

float getMaxVal(int AmideIorII, float graphvals[2][36])   // when press 'c' key
{
  float currentMax = 0.0 ;
  
  for (int i = 0; i < 36; i++)
  {
		if (graphvals[AmideIorII-1][i] > currentMax)
			currentMax = graphvals[AmideIorII-1][i] ;

  }

  return currentMax ;

}



void CalculateMultiElectricVectAngle()   // when press 'c' key
{
	float graphvals[2][36] ;  // 0 is list of C=O band interactions ; 1 is amide II interactions
	float minOut, maxOut ;
	collagenIRInteraction retColInteractn ;
	float s1 = 0.0 ;
	float s2 = 2.5 ;

	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT) ; // this does the clearing to the above clear color


	printf("\nstddev of Gaussian dist = %7.5f", gaussStdev) ;
	if ( (CosOrSin == 0) && (fullOrHalf == 3)) printf("  Split line: Sin distribution\n") ;
	else if ( (CosOrSin == 1) && (fullOrHalf == 3)) printf("  Split line: Cos distribution\n") ;
	else if ( (fullOrHalf == 2)) printf("  Half symmetric distribution\n") ;
	else if ( (fullOrHalf == 1)) printf("  Symmetric distribution\n") ;
	printf("                Amide I                        Amide II \n" ) ;
	printf("<Fibre angle> <IR angle> <Amide I abs>  <Amide II abs> \n" ) ;
	
	angle1_as = 0.0 ;
	for (int i = 0; i <= 95; i += 5)
	{
	//	printf("Angle = %5.1f	",angle1_as) ;
		EVectAngle = 0.0 ;

		for (int i2 = 0; i2 <= 180; i2 += 5 )
		{
		  EVectAngle = i2  ;
			
		  s1 = i / s2 ;
		  if ( fullOrHalf == 1 ) 
		    retColInteractn =	SymetricFibers(angle1_as, 0, newSeed) ;
	  	  else if ( fullOrHalf == 2 ) 
		    retColInteractn =	HalfSymetricFibers(angle1_as, 0, newSeed) ;
		  else if ( fullOrHalf == 3 ) 
		    retColInteractn =	NonSymetricSplitLineFibers(angle1_as, 0.0, newSeed) ;
		  else if ( fullOrHalf == 4 ) 
			retColInteractn = SymetricRandomSplitLineFibers(angle1_as, 0.0, newSeed) ;

		  glutSwapBuffers() ;
	  	  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT) ; // this does the clearing to the above clear color

		  printf("%3.1f  %4.1f %7.5f %7.5f\n",angle1_as,  EVectAngle,  retColInteractn.carbonyl, retColInteractn.amideII ) ;

		  graphvals[0][i2 / 5] =  retColInteractn.carbonyl ;
		  graphvals[1][i2 / 5] =  retColInteractn.amideII ;
		}
		
		minOut = getMinVal(1, graphvals) ;
		maxOut = getMaxVal(1, graphvals) ;
		
//		printf("%2.1f     %7.5f %7.5f %7.5f", angle1_as,  minOut, maxOut ,  maxOut-minOut) ;
		minOut = getMinVal(2, graphvals) ;
		maxOut = getMaxVal(2, graphvals) ;
//		printf("       %7.5f %7.5f %7.5f\n",  minOut, maxOut ,  maxOut-minOut) ;
		
		if (angle1_as == 50) 
			angle1_as = 54.7 ;
		else if ((angle1_as > 50) && (angle1_as < 55))
			angle1_as = 55 ;
		else 
			angle1_as = angle1_as + 5 ;
		
		
	}

	printf("\n\n") ;

}



void keyFunction(unsigned char key, int x, int y)
{
    if(key == '\033')   
        exit(0);
	else if(key == 'x' || key == 'X') 
	{  
		azimuth = 90 ;
		elevation = 270.0 ;
		twist = 90 ;

	}

	else if(key == 'y' || key == 'Y') 
	{ 
		azimuth = 0 ;
		elevation = 270.0 ;	
		twist = 0 ;
	}

	else if(key == 'z' || key == 'Z')  
	{
		azimuth = 0 ;
	    elevation = 0.0 ;
		twist = 0 ;
	}
		
	else if(key == 'r' || key == 'R') 
	{ 
		azimuth = 315.0 ;
	    elevation = 45.0 ;
		twist = 0 ;
	}

	else if(key == 'g' || key == 'G')  
	{ 
		newSeed += 1 ;
	}
		
	else if(key == 'h' || key == 'H') 
	{ 
		newSeed -= 1 ;
	}

	else if(key == 'c' || key == 'C') 
	{ 
		CalculateMultiElectricVectAngle() ;
	}

	else if(key == 'a' || key == 'A') 
	{ 
		angle1_as = angle1_as + 5 ;
		if ( angle1_as > 360) angle1_as = 0 ; }

	else if(key == 's' || key == 'S') 
	{ angle1_as = angle1_as - 5 ;
		if ( angle1_as <= 0) angle1_as = 0 ; }

	else if(key == 'j' || key == 'J') 
	{ angle1_jk = angle1_jk + 5 ;
		if ( angle1_jk > 360) angle1_jk = 0 ; }

	else if(key == 'k' || key == 'K') 
	{ angle1_jk = angle1_jk - 5 ;
		if ( angle1_jk <= 0) angle1_jk = 0 ; }

	else if(key == 'e' || key == 'E') 
	{ EVectAngle = EVectAngle + 15 ;
		if ( EVectAngle >= 360) EVectAngle = 0 ; }

	else if(key == 'w' || key == 'W') 
	{ EVectAngle = EVectAngle - 15 ;
		if ( EVectAngle <= 0) EVectAngle = 360 ; }

	else if(key == 'd' || key == 'D') 
	{ gaussStdev = gaussStdev + 2 ;
		if ( gaussStdev > 400) gaussStdev = 0 ; }

	else if(key == 'f' || key == 'F') 
	{ gaussStdev = gaussStdev - 2 ;
		if ( gaussStdev < 0) gaussStdev = 0 ; }

	else if(key == 't' || key == 'T') 
	{	
		fullOrHalf += 1 ;
		if ( fullOrHalf > 4 ) fullOrHalf = 1 ;
	}
	else if(key == 'l' || key == 'L') 
	{	
		CosOrSin += 1 ;
		if ( CosOrSin > 1 ) CosOrSin = 0 ;
	}
	else if(key == ';' || key == ':') 
	{	
		minLength += 0.1 ;
		if ( minLength > 1 ) minLength = 1 ;
	}
	else if(key == '/' || key == '?') 
	{	
		minLength -= 0.1 ;
		if ( minLength < 0.0 ) minLength = 0 ;
	}
	else if(key == 'v' || key == 'V') 
	{	
		relToPerfectRadDist = !relToPerfectRadDist ;
		if (relToPerfectRadDist == true)
		printf("relToPerfectRadDist = true\n");
		else
		printf("relToPerfectRadDist = false\n");
	}


//	spherEyeCoords = getSphericalCoords(glEyePos) ;

	glutPostRedisplay();

}


void mouseButtonFunction(int button, int state, int xPos, int yPos)
{
	// for tracking mouse //

	if (state == GLUT_DOWN)
	{
		oldPos[0] = xPos;  
		oldPos[1] = yPos; 	
		if (button == GLUT_LEFT_BUTTON)
		{
			leftRightButton.leftButton = buttonDown ;
			//printf("\n GLUT_LEFT_BUTTON" ) ;
		}
		else if (button == GLUT_RIGHT_BUTTON)
		{
			leftRightButton.rightButton = buttonDown ;
			//printf("\n GLUT_RIGHT_BUTTON" ) ;
		}

	}
	else if (state == GLUT_UP)
	{
	
		if (button == GLUT_LEFT_BUTTON)
		{
			leftRightButton.leftButton = buttonUp ;
		//	printf("\n GLUT_LEFT_BUTTON" ) ;
		}
		else if (button == GLUT_RIGHT_BUTTON)
		{
			leftRightButton.rightButton = buttonUp ;
		//	printf("\n GLUT_RIGHT_BUTTON" ) ;
		}
	}
	
}





void mouseMoveFunction(int xPos, int yPos) 
{	
	if (leftRightButton.leftButton == buttonDown && leftRightButton.rightButton == buttonUp) 
	{
		azimuth   -= (float) (oldPos[0] - xPos) ;
		elevation += (float) (oldPos[1] - yPos) ;

		if (azimuth > 360) azimuth = 0 ;
		if (azimuth < 0)   azimuth = 360 ;


		if (elevation > 360) elevation = 0 ;
		if (elevation < 0)   elevation = 360 ;

		spherEyeCoords.theta  += ((float)(oldPos[0] - xPos) * 0.01) ;

		spherEyeCoords.phi += ((float)(oldPos[1] - yPos) * 0.01) ;

	
	} 
	else if (leftRightButton.leftButton == buttonUp && leftRightButton.rightButton == buttonDown) 
	{					
		twist += (float) (oldPos[0] - xPos) ;

		if (twist > 360) twist = 0 ;
		if (twist < 0)   twist = 360 ;	

	}
	else if (leftRightButton.leftButton == buttonDown && leftRightButton.rightButton == buttonDown) 
	{
		
		
	}
	oldPos[0] = xPos ;
	oldPos[1] = yPos ;
//	printf("\n phi = %5.3f, theta = %5.3f" , spherEyeCoords.phi, spherEyeCoords.theta) ;

//	glEyePos = getXYZFromSpherCoords(spherEyeCoords) ;

	glutPostRedisplay(); 
}


double dToR(double deg)
{
	double rad ;
	rad = (deg / 180) * M_PI  ;
	return rad ;
}




/*****************************************************************
* Function name	: ChangeViewport				     			  
* Description	: Callback function for glutReshapeFunc(). Handles 
				  resizing events    
* Return type	: void										     
* Argument      : void										      
******************************************************************/
void ChangeViewport(GLsizei w, GLsizei h)
{
	if (h==0) h=1 ;

	widthScreen = w ;
	heightScreen = h ;

	glMatrixMode(GL_PROJECTION);	// change matrix mode to effect projection matrix
	glLoadIdentity();				// load identity matrix for projection

	if (w <= h)
	//	glOrtho(-120.0, 120.0, -120.0, 120.0*h/w, 120.0, -120.0) ;  // left,right, bottom,top, near,far
		gluPerspective(70.0,w/h,0.2,1000) ;
	else
	//	glOrtho(-120.0, 120.0*w/h, -120.0, 120.0, 120.0, -120.0) ;
		gluPerspective(70.0,w/h,0.2,1000) ;

	glViewport(0,0,w, h) ;
	
	glMatrixMode(GL_MODELVIEW) ;	// revert back to model view matrix mode
	glLoadIdentity() ;				// load identity matrix

	
	glClearColor(1.0,1.0,1.0,1.0); // background color of the context
}





/******************************************************************
* Function name	: initialiseGLState				     			 
* Description	: Sets common OpenGL rendering states.			  
*				  Can be called to re-set parameters based on     
*				  global variables to enable/disable things       
* Return type	: void										      
* Argument      : void										      
******************************************************************/
void  initialiseGLState(void)
{
	// set initial background color
	glClearColor(1.0f,1.0f,1.0f,1.0f) ;
	// We draw to the back buffer, then we flip to the front
	glDrawBuffer(GL_BACK);
	// Use flat/gouraud shading
	glShadeModel(GL_FLAT);
	// Enable Z-Buffering
	glEnable(GL_DEPTH_TEST);
	// set line width
	glLineWidth(1.0f) ;
	// set point size
	glPointSize(2.0f) ;
	// turn on antialiased point drawing and hint at quality
	glEnable(GL_POINT_SMOOTH);
	glHint(GL_POINT_SMOOTH_HINT, GL_NICEST);
	// turn on antialiased line drawing and hint at quality
	glEnable(GL_LINE_SMOOTH); 
	glHint(GL_POINT_SMOOTH_HINT, GL_NICEST);
	

}



/******************************************************************
* Function name	: DrawXYZScales				     					  
* Description	: Draws XYZ scales of length = parameter 1
*				  with full line in +ve direction and dashed line in -ve direction
*			      red = x, green = y, blue = z
*				  Usefull for debuging					          	
* Return type	: void										      	
* Argument      : length - length of line for scales in object space
*					  
******************************************************************/
void DrawXYZScales(float length)
{
	glBegin(GL_LINES) ; 
		glColor3f(0.0f, 1.0f, 0.0f) ; // X
		glVertex3f(0,0,0) ; glVertex3f(length,0,0) ;
		glColor3f(0.0f, 1.0f, 0.0f) ; // Y
		glVertex3f(0,0,0) ; glVertex3f(0,length,0) ;	
		glColor3f(0.0f, 1.0f, 0.0f) ; // Z
		glVertex3f(0,0,0) ; glVertex3f(0,0,length) ;
	glEnd() ;

	glEnable(GL_LINE_STIPPLE) ;
	glLineStipple(1,255) ;

	glBegin(GL_LINES) ;
		glColor3f(0.0f, 1.0f, 0.0f) ; // X
		glVertex3f(0,0,0) ; glVertex3f(-length,0,0) ;
		glColor3f(0.0f, 1.0f, 0.0f) ; // Y
		glVertex3f(0,0,0) ;	glVertex3f(0,-length,0) ;
		glColor3f(0.0f, 1.0f, 0.0f) ; // Z
		glVertex3f(0,0,0) ;	glVertex3f(0, 0,-length) ;
	glEnd() ;

	glLineStipple(1,65535) ;
	glDisable(GL_LINE_STIPPLE) ;
}


// not used:
/******************************************************************
* Function name	: getXYZFromSpherCoords				     					  
* Description	: Converts an spherical (rho, theta, phi) coord to (x,y,z) cartesian equivalent
*				  				          	
* Return type	: struct vect2D in - cartesian (x,y,z) coordiantes of a point 										      	
* Argument      : struct spherCoords - spherical (rho, theta, phi) coordinate to convert 
*					  
******************************************************************/
struct vect2D getXYZFromSpherCoords(struct spherCoords in) 
{
	struct vect2D ret ;
	double sinPhi, sinTheta, cosTheta, cosPhi ;
	sinTheta	= sin(in.theta) ;
	sinPhi		= sin(in.phi) ;
	cosTheta	= cos(in.theta) ;
	cosPhi		= cos(in.phi) ;

	ret.x = in.dist * sinPhi * cosTheta ;
	ret.y = in.dist * sinPhi * sinTheta ;
	ret.z = in.dist * cosPhi ;
//	fprintf(stderr, "eyeDist = %f, theta = %f, phi = %f\n", ret.eyeDist, rToD((double)ret.theta), rToD((double)ret.phi) ) ;
	return (ret) ;
}


// not used:
/******************************************************************
* Function name	: getSphericalCoords				     					  
* Description	: Converts an (x,y,z) coordinate to spherical (rho, theta, phi) equivalent
*				  				          	
* Return type	: struct spherCoords - spherical (rho, theta, phi) coordinates 										      	
* Argument      : struct vect2D in - cartesian (x,y,z) coordiantes of a point
*					  
******************************************************************/
struct spherCoords getSphericalCoords(struct vect2D in) 
{
	struct spherCoords ret ;
	ret.dist = sqrt((double) ( (in.x*in.x) + (in.y*in.y) + (in.z*in.z) )) ;
	ret.theta = atan2(in.y, in.x) ;
	ret.phi = acos( in.z / ret.dist ) ;
//	fprintf(stderr, "eyeDist = %f, theta = %f, phi = %f\n", ret.eyeDist, rToD((double)ret.theta), rToD((double)ret.phi) ) ;
	return (ret) ;
}






/*****************************************************************
* Function name	:	SymetricFibersOld		 - replaced with SymetricSplitLineFibers code to test we end up with same result		     			  
* Description	:   Draws radialy symetric array of fibers each with
*					  radially symmetric bonds at 90 deg to fibre.   
*					Calculates an 'interaction factor' that is between 0-1 indicating how much a polarised IR beam 
*					interacts with 
* Return type	:	void - sets global variable glEyePos 										     
* Argument      :	void										      
******************************************************************/

collagenIRInteraction SymetricFibersOld(float angleWRTZ1, float angleWRTZ2, long seedIn) 
{
	struct spherCoords colBackbone ;
	colBackbone.dist =  BACKBONELENGTH ;
	float spacing =  1.25  ;
	int numBonds = 8 ;
	int numFibers= 18 ; // 18 ;
	static long seed = 554433341 ;
	float randomGauss1, randomGauss2 ;
	float eVectPos = 0, eVectLength = 0 ;

	seed = seedIn ;

	double sectionCarbonylThetaSqrd = 0 ;
	double sectionAmideIIThetaSqrd = 0 ;
	collagenIRInteraction retColInteractn ;
	retColInteractn.amideII  = 0;
	retColInteractn.carbonyl = 0 ;

 
    float angle2 = 0 ;  // this is angle around the z axis (as compared with: angleWRTZ, the angle with the z axis)
	float stepAngle = 2* angle2 / (numFibers -1) ;
	
	
	DrawXYZScales(40.0) ;


	// draw cos curve over plane of IR beam
	glBegin( GL_LINES ) ;
	for (int i = 0 ; i < 40 ; i++ )
	{
	// draw the E vector
    	glVertex3f(cos(dToR(EVectAngle))*eVectLength,sin(dToR(EVectAngle))*eVectLength,eVectPos) ; glVertex3f(0,0,eVectPos) ;
		eVectPos +=  1.0 ;
		eVectLength = cos(eVectPos / 4) * 10 ;
	}
	glEnd() ;

//	if (leftRightButton.leftButton == buttonUp  && leftRightButton.rightButton == buttonUp && angleWRTZ2 != 0.0 ) 
//	printf("\nE vect: %7.5f  %7.5f  : ", EVectAngle, angleWRTZ1 ) ;
	for ( i = 0 ; i < numFibers ; i++ )
	{

		glColor3f(1.0f, 0.0f, 0.0f) ;

		glPushMatrix( ) ;
			randomGauss1 = ranGauss(&seed) * gaussStdev ;
;

			glRotatef(90,0.0,1.0,0.0) ;
			glRotatef(angle2,0.0,0.0,1.0) ;  // this is number that goes around z axis
			// angle1 is angle fiber long axis makes with z axis
			glRotatef(angleWRTZ1+randomGauss1,0.0,1.0,0.0) ;  // this is angle wrt z axis (splay angle)  
		  	
		  // repeat identical transform but without the gluLookAt
		  glPushMatrix( ) ; 
		    glLoadIdentity();  // start from fresh matrix
			glRotatef(90,0.0,1.0,0.0) ;
			glRotatef(angle2,0.0,0.0,1.0) ;  // this is number that goes around z axis
			// angle1 is angle fiber long axis makes with z axis
			glRotatef(angleWRTZ1+randomGauss1,0.0,1.0,0.0) ;  // this is angle wrt z axis (splay angle)   
		
			glGetFloatv(GL_MODELVIEW_MATRIX, modelview) ;
		//	printModelViewMat(modelview) ;
		  glPopMatrix() ;

		  
		    retColInteractn  = drawCollagenBonds( spacing , 5, colBackbone, numBonds, EVectAngle ) ;
		    sectionCarbonylThetaSqrd += retColInteractn.carbonyl ;
		    sectionAmideIIThetaSqrd  += retColInteractn.amideII  ;
	    glPopMatrix() ;
	   	 

		angle2 = angle2 + (360 / numFibers) ;
		if ( angle2 > 360) angle2 = 0 ;


	}

    int numFiberPositions = (int) (colBackbone.dist / spacing ) ;
    numFiberPositions++ ;



	float ratio ;
	sectionCarbonylThetaSqrd = (sectionCarbonylThetaSqrd / ( numBonds * numFiberPositions * numFibers)) ;
	sectionAmideIIThetaSqrd  = sectionAmideIIThetaSqrd / numFibers ;
	ratio = sectionCarbonylThetaSqrd / sectionAmideIIThetaSqrd ;

	if (leftRightButton.leftButton == buttonUp  && leftRightButton.rightButton == buttonUp && angleWRTZ2 != 0.0 ) 
	{printf("%7.5f %7.5f %7.5f\n", sectionCarbonylThetaSqrd , sectionAmideIIThetaSqrd , ratio) ;}

	retColInteractn.carbonyl = sectionCarbonylThetaSqrd ;
	retColInteractn.amideII  = sectionAmideIIThetaSqrd  ;

	return (retColInteractn) ;

}



