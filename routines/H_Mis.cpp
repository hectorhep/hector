
// ROOT includes
#include "TCanvas.h"
#include "TLegend.h"
#include "TFile.h"
#include "TTree.h"
#include "THStack.h"
#include "TGraph.h"
#include "TGraph2D.h"
#include "TGraphErrors.h"
#include "TMultiGraph.h"
#include "TF1.h"
#include "TLorentzVector.h"

// local includes
#include "H_BeamLine.h"
#include "H_BeamParticle.h"
#include "H_RomanPot.h"
#include "H_RecRPObject.h"
#include <iostream>
using namespace std;



/// choice between shift and tilt
const bool shift = false; // if false, it is then tilt

/// Impact on E reconstruction of the misalignment of each quadrupole (separately), for a RP at 420.
void misalignment_quad_bias_old(float energy = 100., float displacement = 0.0005, bool corrected = true, int method = TM, bool testtag = false) {
	/// @param energy : gamma enegy
	/// @param displacement : amplitude of the misalignment (\f$\mu m\f$ for shifts, \f$\mu rad\f$ for tilts)
	/// @param method : reconstruction method (TM, ACM)
	/// @param testtag : checks if the particle is still in the RP acceptance
	extern bool relative_energy;
	relative_energy = true;
	if(!relative_energy) {
		cout << "You should be in relative energy" << endl;
		return;
	}

	float x1,y1,x2,y2;
	float x01,y01,x02,y02;

	H_BeamLine* beam1 = new H_BeamLine( 1,450);
	beam1->fill("data/LHCB1IR5_v6.500.tfs");
	H_RecRPObject recp1(420,428,*beam1);
	H_RomanPot * rp = new H_RomanPot("RP",420.,3500.);
	beam1->add(rp);

	const int nelements_1 = beam1->getNumberOfElements();

	int qi = 0;

	int nquadrupoles_1 = 0;
	for(int i = 0; i < nelements_1; i++) {
		if(beam1->getElement(i)->getType() == VQUADRUPOLE || beam1->getElement(i)->getType() == HQUADRUPOLE) {
			nquadrupoles_1++;
		}
	}
	cout<<"Testing "<<nquadrupoles_1<<" quadrupoles of the beam"<<endl;
	float s1[nquadrupoles_1];
	float e1p[nquadrupoles_1];
	float e1m[nquadrupoles_1];

	cout<<"---------------------------"<<endl;
	cout<<"Checking beam 1 quadrupoles"<<endl;
	cout<<"---------------------------"<<endl;
	for(int i = 0; i < nelements_1; i++) {
		if(beam1->getElement(i)->getType() == VQUADRUPOLE || beam1->getElement(i)->getType() == HQUADRUPOLE) {
			cout<<beam1->getElement(i)->getName();
			s1[qi] = beam1->getElement(i)->getS();
			H_BeamParticle p1; /* test proton */
			H_BeamParticle p0; /* reference proton, no energy loss */
			p1.setE(7000. - energy);
			if (shift) beam1->alignElement(beam1->getElement(i)->getName(), displacement,0);
			else 	beam1->tiltElement(beam1->getElement(i)->getName(), displacement,0);		
			p1.computePath(beam1,1);
			p0.computePath(beam1,1);
			if(p1.stopped(beam1) || !testtag) {
				if(p1.getStoppingElement()->getName()>="RP" || !testtag) {
					cout<<"\t d+ check";
					p1.propagate(420);
					p0.propagate(420);
					x1 = p1.getX();
					x01 = p0.getX();
					y1 = p1.getY();
					y01 = p0.getY();
					p1.propagate(428);
					p0.propagate(428);
					x2 = p1.getX();
					x02 = p0.getX();
					y2 = p1.getY();
					y02 = p0.getY();
					if(corrected) {
						recp1.setPositions(x1-x01,y1-y01,x2-x02,y2-y02);
					} else {
						recp1.setPositions(x1,y1,x2,y2);
					}
					e1p[qi] = recp1.getE(method)/energy;
				} else {
					e1p[qi] = 0;
					cout<<"\t d+ no hit";
				}
			}
			if (shift) beam1->alignElement(beam1->getElement(i)->getName(),-2*displacement,0);
			else beam1->tiltElement(beam1->getElement(i)->getName(),-2*displacement,0);
			p1.resetPath();
			p0.resetPath();
			p1.computePath(beam1,1);
			p0.computePath(beam1,1);
			if(p1.stopped(beam1) || !testtag) {
				if(p1.getStoppingElement()->getName()>="RP" || !testtag) {
					cout<<"\t d- check"<<endl;;
					p1.propagate(420);
					p0.propagate(420);
					x1 = p1.getX();
					x01 = p0.getX();
					y1 = p1.getY();
					y01 = p0.getY();
					p1.propagate(428);
					p0.propagate(428);
					x2 = p1.getX();
					x02 = p0.getX();
					y2 = p1.getY();
					y02 = p0.getY();
					if(corrected) {
						recp1.setPositions(x1-x01,y1-y01,x2-x02,y2-y02);
					} else {
						recp1.setPositions(x1,y1,x2,y2);
					}
                    e1m[qi] = recp1.getE(TM)/energy;
                } else {
					e1m[qi] = 0;
					cout<<"\t d- no hit"<<endl;
				}
            }
			if (shift) beam1->alignElement(beam1->getElement(i)->getName(), displacement,0);
			else beam1->tiltElement(beam1->getElement(i)->getName(), displacement,0);
			qi++;
		}
	}
	// drawing
	
	TLegend* leg=new TLegend(0.68,0.75,0.99,0.88);
	TCanvas* c1 = new TCanvas("c1","Misalignment of quadrupoles",1);
	c1->SetGrid();
	c1->cd();
	TGraph* g1p = new TGraph(nquadrupoles_1,s1,e1p);
	g1p->SetLineStyle(1);
	g1p->SetLineColor(1);
	g1p->SetMarkerColor(1);
	g1p->SetMarkerStyle(22);
	leg->AddEntry(g1p,"positive diplacement","pl");
	TGraph* g1m = new TGraph(nquadrupoles_1,s1,e1m);
	g1m->SetLineStyle(2);
	g1m->SetLineColor(2);
	g1m->SetMarkerColor(2);
	g1m->SetMarkerStyle(23);
	leg->AddEntry(g1m,"negative diplacement","pl");
	char mytitle[150];
	if (shift) sprintf(mytitle,"Reconstructed energy for quadrupole displacement of %.0f#mum (E_{#gamma} = %.0f GeV)",displacement*1000000.,energy);
	else sprintf(mytitle,"Reconstructed energy for quadrupole displacement of %.0f#murad (E_{#gamma} = %.0f GeV)",displacement*1000000.,energy);
	TMultiGraph* mg1 = new TMultiGraph("mg1",mytitle);
	mg1->Add(g1p,"lp");
	mg1->Add(g1m,"lp");
	mg1->Draw("a");
	mg1->GetYaxis()->SetTitle("E_{rec}/E_{#gamma}");
	mg1->GetYaxis()->SetTitleOffset(1.15);
	mg1->GetXaxis()->SetTitle("s [m]");
	leg->SetBorderSize(1);
	leg->SetFillColor(0);
	leg->Draw();

	c1->Modified();
	c1->Update();
	
	return;
}
