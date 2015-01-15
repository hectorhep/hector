/*#include "H_RectangularAperture.h"
#include "H_EllipticAperture.h"
#include "H_RectEllipticAperture.h"
#include "H_CircularAperture.h"
#include "H_SectorDipole.h"
*/
#include "H_BeamLine.h"
//#include "H_BeamParticle.h"
#include <iostream>
using namespace std;

int main() {
 /* 
	H_Aperture h0; h0.printProperties();
	H_RectangularAperture h1(1,2,3,4); h1.printProperties();
	H_EllipticAperture h2(11,22,33,44); h2.printProperties();
	H_CircularAperture h3(111,222,333); cout << h3;
	H_RectEllipticAperture h4(2,3,4,5,6,7); cout << h4;

	H_Aperture * ap;
	ap = new H_RectEllipticAperture(h4);
	cout << "essai1\n";
	ap->printProperties();
	cout << "essai2\n";
	cout << *ap;
cout << endl;
	H_SectorDipole * s = new H_SectorDipole(1,2,3);
	s->setAperture(ap);
	cout << *s << endl;
	

	delete ap;
*/

  H_BeamLine *beamline1 = new H_BeamLine(1,420);
  beamline1->fill("data/LHCB1IR5_7TeV.tfs",1,"IP5");
  cout << *beamline1;
/*
  H_BeamParticle *p = new H_BeamParticle;
  p->smearAng(1,1);
  p->smearPos(2,2);
  p->computePath(beamline1);
  if(p->stopped(beamline1)) cout << "stopped\n";
  float a= p->getE();
  cout << a << endl;
  a++;
  cout << a << endl;

  delete p;
*/
  delete beamline1;

  return 0;
}
