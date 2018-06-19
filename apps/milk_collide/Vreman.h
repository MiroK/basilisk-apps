// Solving for filtered quantity
/* Eddy viscosity from velocityfield according to SGS-model proposed 
by Vreman in Physics of Fluids (13/10) from 2004 under the title 
“An eddy-viscosity subgrid-scale model for turbulent shear flow: Algebraic theory and applications”.

Important aspect of this closure formulation is that the SGS-fluxes vanish 
in (some) cases of a laminar flow (unlike Standard Smagorinsky). Function 
requires a double ‘Cs’ for (classic) Smagorinsky constant, (3D) flow field 
u and the molecular viscosity as input and writes the Eddyviscosity field 
‘Evis’ (i.e. the eddy viscosity for momentum). This function returns centered 
(i.e. scalar) values of eddy viscosity field rather than in face vector form. 
Therefore, in the “.c” script where you call this function you way want to 
average between two centered values on either side of the respective face.
*/

#include "common.h"
                                      
void eddyviscosity(double Cs, vector u, scalar Evis){	
  double d1v1, d2v1, d3v1, d1v2, d2v2, d3v2, d1v3, d2v3, d3v3; 
  double b11, b12, b13, b22, b23, b33;
  double abeta, bbeta;	
  foreach(){
    d1v1=(u.x[1,0,0]-u.x[-1,0,0])/2/Delta;
    d2v1=(u.x[0,1,0]-u.x[0,-1,0])/2/Delta;
    d3v1=(u.x[0,0,1]-u.x[0,0,-1])/2/Delta;
    d1v2=(u.y[1,0,0]-u.y[-1,0,0])/2/Delta;
    d2v2=(u.y[0,1,0]-u.y[0,-1,0])/2/Delta;
    d3v2=(u.y[0,0,1]-u.y[0,0,-1])/2/Delta;
    d1v3=(u.z[1,0,0]-u.z[-1,0,0])/2/Delta;
    d2v3=(u.z[0,1,0]-u.z[0,-1,0])/2/Delta;
    d3v3=(u.z[0,0,1]-u.z[0,0,-1])/2/Delta;
    b11 = Delta*Delta*(d1v1*d1v1+d2v1*d2v1+d3v1*d3v1);
    b12 = Delta*Delta*(d1v1*d1v2+d2v1*d2v2+d3v1*d3v2);
    b13 = Delta*Delta*(d1v1*d1v3+d2v1*d2v3+d3v1*d3v3);
    b22 = Delta*Delta*(d1v2*d1v2+d2v2*d2v2+d3v2*d3v2);
    b23 = Delta*Delta*(d1v2*d1v3+d2v2*d2v3+d3v2*d3v3);
    b33 = Delta*Delta*(d1v3*d1v3+d2v3*d2v3+d3v3*d3v3);
    abeta = sq(d1v1)+sq(d2v1)+sq(d3v1)+
      sq(d1v2)+sq(d2v2)+sq(d3v2)+
      sq(d1v3)+sq(d2v3)+sq(d3v3);
    bbeta =	b11*b22-sq(b12)+b11*b33-sq(b13)+b22*b33-sq(b23);

    double eddy = abeta>0 ? ( 2.5*Cs*Cs*sqrt(bbeta/(abeta))) : 0.;
    Evis[] = max(eddy,  0);
  }	
}
