  /* * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *                                                         *
*                   --<--<--  A fast simulator --<--<--     *
*                 / --<--<--     of particle   --<--<--     *
*  ----HECTOR----<                                          *
*                 \ -->-->-- transport through -->-->--     *
*                   -->-->-- generic beamlines -->-->--     *
*                                                           *
* JINST 2:P09005 (2007)                                     *
*      X Rouby, J de Favereau, K Piotrzkowski (CP3)         *
*       http://www.fynu.ucl.ac.be/hector.html               *
*                                                           *
* Center for Cosmology, Particle Physics and Phenomenology  *
*              Universite catholique de Louvain             *
*                 Louvain-la-Neuve, Belgium                 *
 *                                                         *
   * * * * * * * * * * * * * * * * * * * * * * * * * * * */

// ROOT includes
#include "TFile.h"
#include "TTree.h"
#include "TF1.h"
#include "TLorentzVector.h"
#include "TCanvas.h"

// Hector
#include "H_BeamLine.h"
#include "H_BeamParticle.h"
#include "H_RecRPObject.h"
#include "H_RomanPot.h"

//C++
#include <string>
using namespace std;

TF1* calibration_plot(float displacement = 0., int color = 1, char el_name[50] = "\"MQM.9R5.B1\"", string file  = "data/lpair_mumu_2gev.root") {

	vector<float> rec_e;
	vector<float> x_rp;

	H_BeamLine* beam1 = new H_BeamLine( 1,500);
	beam1->fill("data/LHCB1IR5_v6.500.tfs");
	H_RomanPot* rp = new H_RomanPot("RP",420.,3500.);
	beam1->add(rp);
	beam1->alignElement(el_name,displacement,0);

	
	// reading yanwen file
	TF1 * fake = new TF1("fake","1",0,1); // in case of trouble
	TFile* calfile = new TFile(file.c_str());
		if(calfile->IsZombie()) {cout << "Can not open the file : " << file << endl; return fake;}
	TTree* caltree = (TTree*) calfile->Get("h101");
		if(!caltree) {cout << "Can not open the tree.\n"; return fake;}
	delete fake; // forget about it
	

	float phep[4][5];
	caltree->SetBranchAddress("Phep",phep);
	/* phep[0] & phep[3] = protons
	 * phep[1] & phep[2] = muons */

	TLorentzVector mu1, mu2, cms;
	float m_mumu, pz_mumu;
	float gammae_1_central, gammae_2_central;
	float rho;
	int nentries = caltree->GetEntriesFast();
	for (int jentry=0; jentry<nentries;jentry++) {
		caltree->GetEntry(jentry);
		mu1.SetXYZM(phep[1][0],phep[1][1],phep[1][2],0);
		mu2.SetXYZM(phep[2][0],phep[2][1],phep[2][2],0);
		cms = mu1;
		cms += mu2;
		m_mumu = cms.M();
		pz_mumu = cms.Pz();
		rho = sqrt(m_mumu*m_mumu + pz_mumu*pz_mumu);
		gammae_1_central = (  pz_mumu + rho)/2.;
		gammae_2_central = ( -pz_mumu + rho)/2.;

		H_BeamParticle p1;
		p1.smearPos();
		p1.smearAng();
		p1.setE(phep[0][3]);
		p1.computePath(beam1,1);
		if(p1.stopped(beam1)) {
			if(p1.getStoppingElement()->getName()>="RP") {
				p1.propagate(420);
				rec_e.push_back(gammae_1_central);
				x_rp.push_back(p1.getX()/1000.);
			}
		}
	}
	const int npoints = rec_e.size();
	float rece[npoints];
	float xrp[npoints];
	for(int i = 0; i < npoints; i++) {
		rece[i] = rec_e[i];
		xrp[i]  = x_rp[i];
	}
	TGraph* correl = new TGraph(npoints,xrp,rece);
	char title[100];
	sprintf(title,"Calibration of roman pots from central objects variables - %d events",npoints);
	correl->SetTitle(title);
	correl->SetMarkerColor(color);

	correl->Draw("ap");

	sort(x_rp.begin(),x_rp.end());
	TF1* fit1 = new TF1("fit1","[0] + [1] * x + [2] * x * x",x_rp[0],x_rp[npoints-1]);
	correl->Fit("fit1","Q");
//	correl->GetYaxis()->SetTitle("Photon energy from central state (GeV)");
//	correl->GetXaxis()->SetTitle("X position of proton at RP (mm)");

	TCanvas * raoul = new TCanvas();
	correl->Draw("ap");
	return fit1;
}
