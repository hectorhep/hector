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

/// \file H_Misalignment.cpp
/// \brief Collection of routines about the misalignment impact on reconstruction (see also H_Reconstruction.cpp)

// todo : verifier l'efficacite pour la méthode de reconstructrion avec p0 et p1 
// l 211

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
#include "H_Beam.h"
#include "H_BeamParticle.h"
#include "H_RecRPObject.h"
#include "H_RomanPot.h"
//#include "H_Acceptance.h"

/// choice between shift and tilt
const bool shift = true; // if false, it is then tilt

const float length= 500.; // beam length
const float rp1_s = 420.; // RP420 pot 1
const float rp2_s = 428.; // RP420 pot 2
const float rp_x = 4000.; // RP420 x (µm)
static int hN =0; // histogram counter -- for a unique TH1 histo name

/// Draws the 1D-acceptance for any optical element, for a given energy range and a given virtuality
TH1F* acceptance_element(H_AbstractBeamLine * beamline, const string element_name, const unsigned int number_of_particles=10, const float E_min=0., const float E_max=150., const unsigned int Ne = 10, const float Q2 =-0.1, const bool NonLinear=true, const int crang_sign  = -1, bool drawlines = false) {
/// @param beamline : beamline to test
/// @param number_of_particles : in the beam
/// @param E_min : min value for the beam particle energy loss
/// @param E_max : max value for  the beam particle energy loss
/// @param Ne : number of steps in energy
/// @param Q2 : virtuality value (t) of the emitted gamma
/// @param NonLinear : (des)activates the NonLinear effects for the beam propagation. Should be always on.
	const bool show_path = false;

	vector<TVectorD> stop_pos;
	vector<TVectorD>::iterator spi_i;

	extern bool relative_energy;
	relative_energy = true;
	if(!relative_energy) {	
		cout << "ERROR : You should be in relative energy" << endl;
		return new TH1F();
	}
	if(!beamline) {
		cout << "ERROR : No beamline ! Empty pointer.\n";
		return new TH1F();
	}
	if(!element_name.size()) {
		cout << "ERROR : Empty element name !\n";
		return new TH1F();
	}
	if(!NonLinear) cout << "WARNING : Linear mode (NonLinear off).\n";

        float energies[Ne]; //, * log10x = new float[Ne];
	for(unsigned int e_i=0; e_i<Ne; e_i++) {energies[e_i] = E_min + (E_max - E_min)/(float)Ne * (e_i+0.5);}

	char can_name[100], title[20];
	sprintf(can_name,"Acceptance of %s",element_name.c_str());
	sprintf(title,"h1D_%d",hN++);
	TH1F * h1D = new TH1F(title,can_name,Ne,E_min,E_max); // 1D-histo at fixed Q^2 = 0
	TMultiGraph* hitx_m = new TMultiGraph("hitx_m","hitx_m");
	TMultiGraph* hity_m = new TMultiGraph("hity_m","hity_m");

	cout << "Generating particles " << endl; 
	{
		vector<H_OpticalElement> stopping_elements;
		vector<H_OpticalElement>::iterator el_i;
		vector<int> number_of_stopped_particles;

		extern int kickers_on;
		kickers_on =1;

		for (unsigned int e_i=0; e_i<Ne; e_i++) {
			H_Beam mybeam;
			mybeam.setPosition(PX,PY,crang_sign*CRANG,0.,0.);
			mybeam.setE(BE);
			if(drawlines) {
				mybeam.setDE(0);
				mybeam.setDX(0);
				mybeam.setDY(0);
				mybeam.setDTX(0);
				mybeam.setDTY(0);
			}
			mybeam.createBeamParticles(number_of_particles); // smearings

			for(unsigned int particle_i=0; particle_i<number_of_particles; particle_i++) {
				mybeam.getBeamParticle(particle_i)->emitGamma(energies[e_i],Q2); 
			}
			mybeam.computePath(beamline,NonLinear);
			if(drawlines) {
//				hitx_m->Add(mybeam.drawBeamX(1),"l");
//				hity_m->Add(mybeam.drawBeamY(1),"l");
			}

			vector<TVectorD> stop_pos_i = mybeam.getStoppingElements(beamline,stopping_elements,number_of_stopped_particles);
			for(spi_i = stop_pos_i.begin(); spi_i<stop_pos_i.end(); spi_i++) {
	//			if((*spi_i)[4]>=beamline->getElement(element_name)->getS() && (*spi_i)[4]<=beamline->getElement(element_name)->getS()+beamline->getElement(element_name)->getLength()) {
					stop_pos.push_back(*spi_i);
	//			}
			}

			for(el_i = stopping_elements.begin(); el_i < stopping_elements.end(); el_i++) {
				if((el_i->getName().c_str() == element_name) &&  (mybeam.getBeamParticle(1)->isPhysical())) 
					h1D->Fill(energies[e_i],number_of_stopped_particles[el_i-stopping_elements.begin()]/(float)(number_of_particles));
			}
		} // for e_i
	
        }// kills stopping_elements and other vector<...>'s
	cout << "\t done"<<endl;

	TGraph* svse = new TGraph();
	TGraph2D* hits = new TGraph2D();
	TGraph* hitx = new TGraph();
	if(drawlines) hitx->SetMarkerColor(2);
	TGraph* hity = new TGraph();
	if(drawlines) hity->SetMarkerColor(2);
	TGraph* hitp = new TGraph();
	for(spi_i = stop_pos.begin(); spi_i<stop_pos.end(); spi_i++) {
		svse->SetPoint(spi_i - stop_pos.begin(),(*spi_i)[4],energies[(spi_i - stop_pos.begin())/(int)number_of_particles]);
		hits->SetPoint(spi_i - stop_pos.begin(),(*spi_i)[0],(*spi_i)[4],(*spi_i)[2]);
		hitx->SetPoint(spi_i - stop_pos.begin(),(*spi_i)[4],(*spi_i)[0]);
		hity->SetPoint(spi_i - stop_pos.begin(),(*spi_i)[4],(*spi_i)[2]);
		hitp->SetPoint(spi_i - stop_pos.begin(),(*spi_i)[0],(*spi_i)[2]);
	}
	if(drawlines) {
		hitx->SetMarkerStyle(7);
		hitx->SetMarkerColor(2);
		hity->SetMarkerStyle(7);
		hity->SetMarkerColor(2);
	}
	hitx_m->Add(hitx,"p");
	hity_m->Add(hity,"p");
	TCanvas* can = new TCanvas("can","can",1);
	can->Divide(2,2);
	can->cd(1);
	hits->Draw("P");
	can->cd(2);
	hitx_m->Draw("ap");
	can->cd(3);
	hity_m->Draw("ap");
	can->cd(4);
	hitp->Draw("ap");
	TCanvas* ca2 = new TCanvas("ca2","ca2",1);
	ca2->SetGrid();
	svse->Draw("ap+");
	double height = (E_max - E_min)/10.;
	relative_energy = false;
	beamline->offsetElements(0,5*height/URAD);
	beamline->drawX(1*height,1*height,1);
	beamline->draw();
	svse->SetMarkerStyle(6);
	svse->GetXaxis()->SetTitle("s (m)");
	svse->GetYaxis()->SetTitle("Eloss (GeV)");
	svse->Draw("p+");

	// shows the beam path, just in case
	if (show_path) {
		H_BeamParticle p1;
		p1.setE(BE);
		p1.setPosition(p1.getX()+PX,p1.getY()+PY,p1.getTX()+crang_sign*CRANG,p1.getTY(),0.);
		p1.computePath(beamline,NonLinear);
		TGraph * ppath = p1.getPath(0,kGreen);
		TCanvas * path = new TCanvas("path","path");
		path->cd();
		ppath->Draw("ALP");
	}
	return h1D;
} // acceptance_element


/// Performs the RP calibration from a \f$ \gamma \gamma \rightarrow \mu \mu \f$ sample, in order to correct for the misalignement
TF1* misalignment_calibration_plot(const float displacement = 0., const int color = 1, const char el_name[50] = "\"MQM.9R5.B1\"",const int beam = 1, const int events = -1, const char file[50] = "data/lpair_mumu_2gev.root") {
	/// @param displacement : \f$ \mu m\f$
	/// @param color  for display
	/// @param el_name  : displaced optical element name
	/// @param beam : beam ID (1 or 2)
	/// @param events : max number of events for the calibration run. (-1 means all the sample)
	/// @param file : \f$ \gamma \gamma \rightarrow \mu \mu \f$ sample
	vector<float> rec_e;
	vector<float> x_rp;

	H_BeamLine* beam1;
	if(beam==1) { 
		beam1 = new H_BeamLine( 1,length);
		beam1->fill("data/LHCB1IR5_v6.500.tfs");
	} else {
		beam1 = new H_BeamLine( -1,length);
		beam1->fill("data/LHCB2IR5_v6.500.tfs");
	}
	H_RomanPot* rp = new H_RomanPot("RP",rp1_s,rp_x);
	beam1->add(rp);
	if (shift) beam1->alignElement(el_name,displacement,0);
        else    beam1->tiltElement(el_name,displacement,0);



	// reading yanwen file
	TFile* calfile = new TFile(file);
	if(calfile->IsZombie()) {cout << "Can not read the file " << calfile << endl; return 0;}
	TTree* caltree = (TTree*) calfile->Get("h101");
	if(!caltree) {cout << "Can not read the tree.\n"; return 0;}

	float phep[4][5];
	caltree->SetBranchAddress("Phep",phep);
	/* phep[0] & phep[3] = protons
	 * phep[1] & phep[2] = muons */

	TLorentzVector mu1, mu2, cms;
	float m_mumu, pz_mumu;
	float gammae_1_central, gammae_2_central;
	float rho;
	int nentries = caltree->GetEntriesFast();
	if(events>0) {
		if(events > nentries) {
			cout<<"file only contains "<<nentries<<" events, using all of them"<<endl;
		} else {
			nentries = events;
		}
	}
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
				p1.propagate(rp1_s);
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
	correl->GetYaxis()->SetTitle("Photon energy from central state (GeV)");
	correl->GetXaxis()->SetTitle("X position of proton at RP (mm)");
	sort(x_rp.begin(),x_rp.end());
	TF1* fit1 = new TF1("fit1","[0] + [1] * x + [2] * x * x + [3] * x * x * x",x_rp[0],x_rp[npoints-1]);
	correl->Fit("fit1","Q");
	Double_t* params = fit1->GetParameters();
	Double_t* errors = fit1->GetParErrors();
	cout << "Calibration based on " << npoints << " events.\n";
	cout<<"Parameters : "<<params[0]<<"\t"<<params[1]<<"\t"<<params[2]<<endl;
	cout<<"Errors     : "<<errors[0]<<"\t"<<errors[1]<<"\t"<<errors[2]<<endl;
	//TCanvas *robert = new TCanvas();
	//robert->cd();
	//correl->Draw("ap");
	return fit1;
}

/// Impact on E reconstruction of the misalignment of each quadrupole (separately), for a RP at 420.
void misalignment_quad_bias(float energy = 100., float displacement = 0.0005, int corrected = 1, int method = TM, bool testtag = false, bool acceptance = false, const float seuil = 0.9) {
	/// @param energy : gamma enegy
	/// @param displacement : amplitude of the misalignment (\f$\mu m\f$ for shifts, \f$\mu rad\f$ for tilts)
	/// @param method : reconstruction method (TM, ACM)
	/// @param corrected : 0 = no ; 1 = BPM ; 2 = calibration
	/// @param testtag : checks if the particle is still in the RP acceptance
	extern bool relative_energy;
	relative_energy = true;
	if(!relative_energy) {
		cout << "You should be in relative energy" << endl;
		return;
	}
	int inc = 0;
	char title[500], filename[500];

	float x1,y1,x2=0,y2=0;
	float x01,y01,x02,y02;
	const int number_of_particles = 100;
	const float E_min = 0;
	const float E_max = 100;
	const int Ne = 50;
	const float Q2 = -1E-4;

	H_BeamLine* beam1 = new H_BeamLine( 1,rp2_s+0.1);
	beam1->fill("data/LHCB1IR5_v6.500.tfs");
	H_RecRPObject recp(rp1_s,rp2_s,beam1);
	string rp_name = "RP";
	H_RomanPot * rp = new H_RomanPot(rp_name.c_str(),rp1_s,rp_x);
	beam1->add(rp);

	const int nelements_1 = beam1->getNumberOfElements();

	int qi = 0, nquadrupoles_1 = 0; // iterator on quadrupoles ; number of quadrupoles 
	for(int i = 0; i < nelements_1; i++) {
		if(beam1->getElement(i)->getType() == VQUADRUPOLE || beam1->getElement(i)->getType() == HQUADRUPOLE) {
			nquadrupoles_1++;
		}
	}
	cout<<"Testing "<<nquadrupoles_1<<" quadrupoles of the beam"<<endl;
	float s1[nquadrupoles_1], e1p[nquadrupoles_1], e1m[nquadrupoles_1];
	float acc_meanm[nquadrupoles_1], acc_sigmam[nquadrupoles_1];
	float acc_meanp[nquadrupoles_1], acc_sigmap[nquadrupoles_1];

	for(int i = 0; i < nelements_1; i++) {
		if(beam1->getElement(i)->getType() == VQUADRUPOLE || beam1->getElement(i)->getType() == HQUADRUPOLE) {
				
			cout<<beam1->getElement(i)->getName();
			s1[qi] = beam1->getElement(i)->getS();
			H_BeamParticle p1; // * test proton * /
			p1.setPosition(PX,PY,-CRANG,0,0);
			p1.setE(7000. - energy);
			if (shift) beam1->alignElement(beam1->getElement(i)->getName(), displacement,0); // shift or tilt
			else 	beam1->tiltElement(beam1->getElement(i)->getName(), displacement,0);     // ...
			sprintf(title,"%s at %.1f mm",beam1->getElement(i)->getName().c_str(),displacement*1000.);
			sprintf(filename,"acc_%d.eps",++inc);

			if(acceptance) {
				float acc_min = 0, acc_max = 0;
				TH1F * h1D_0 = acceptance_element(beam1,rp_name,number_of_particles,E_min,E_max,Ne,Q2);
				for (int j=0; j<h1D_0->GetNbinsX()-1; j++) {
					float value1 = h1D_0->GetBinContent(j);
					float value2 = h1D_0->GetBinContent(j+1);
					if ((value1 < seuil) && (value2 > seuil)) { // rising edge
						acc_min = E_min + (j+1)*(E_max-E_min)/Ne; 
						cout << "Acceptance (GeV)= [" << acc_min << "\t;";
					}
					else if ((value1 > seuil) && (value2 < seuil)) { // decreasing edge
						acc_max =  E_min + j*(E_max-E_min)/Ne;
						cout << acc_max << "] -- positive displacement \n";
					}
				}
				acc_meanp[qi] = (acc_min + acc_max)/2.;
				acc_sigmap[qi] = acc_meanp[qi] - acc_min;
				if(!h1D_0) return; // in case of troubles
				TCanvas canh1D_0("misaligned_acceptance","");
				canh1D_0.cd();
				h1D_0->SetTitle(title);
				h1D_0->Draw();
				h1D_0->GetXaxis()->SetTitle("GeV");
				h1D_0->GetYaxis()->SetTitle("Acceptance");
				canh1D_0.Print(filename);
				delete h1D_0;
			} // if acceptance

			p1.computePath(beam1,1);
			if(p1.stopped(beam1) || !testtag) {
				if(p1.getStoppingElement()->getName()>=rp_name || !testtag) {
					cout<<"\t d+ check";
					p1.propagate(rp1_s); x1 = p1.getX(); y1 = p1.getY();
					if(corrected==2) { // correction based on the calibration
						TF1* calib_1 = misalignment_calibration_plot(displacement,1,beam1->getElement(i)->getName().c_str());
						if(!calib_1){return;}
						e1p[qi] = calib_1->Eval(p1.getX()/1000.)/energy;
						delete calib_1;
					} else if (corrected==1) { //  correction based on bpm
						H_BeamParticle p0; // * reference proton, no energy loss * /
						p0.setPosition(PX,PY,-CRANG,0,0);
						p0.computePath(beam1,1);
						p0.propagate(rp1_s); x01= p0.getX();y01 = p0.getY();
						p1.propagate(rp2_s); x2 = p1.getX(); y2 = p1.getY();
						p0.propagate(rp2_s); x02= p0.getX();y02 = p0.getY();
						recp.setPositions(x1-x01,y1-y01,x2-x02,y2-y02); 
						e1p[qi] = recp.getE(method)/energy;
					} else {
						recp.setPositions(x1,y1,x2,y2); // no correction
						e1p[qi] = recp.getE(method)/energy;
					}
				} else {
					e1p[qi] = 0;
					cout<<"\t d+ no hit";
				}
			}



			if (shift) beam1->alignElement(beam1->getElement(i)->getName(),-2*displacement,0);
			else beam1->tiltElement(beam1->getElement(i)->getName(),-2*displacement,0);
                        sprintf(title,"%s at %.1f mm",beam1->getElement(i)->getName().c_str(),-displacement*1000.);
                        sprintf(filename,"acc_%d.eps",++inc);
			
			if(acceptance) {
				float acc_min = 0, acc_max = 0;
				TH1F *h1D_0 = acceptance_element(beam1,rp_name,number_of_particles,E_min,E_max,Ne,Q2);
				for (int j=0; j<h1D_0->GetNbinsX()-1; j++) {
                                        float value1 = h1D_0->GetBinContent(j);
                                        float value2 = h1D_0->GetBinContent(j+1);
                                        if ((value1 < seuil) && (value2 > seuil)) { // rising edge
                                                acc_min = E_min + (j+1)*(E_max-E_min)/Ne;
                                                cout << "Acceptance (GeV)= [" << acc_min << "\t";
                                        }
                                        else if ((value1 > seuil) && (value2 < seuil)) { // decreasing edge
                                                acc_max =  E_min + j*(E_max-E_min)/Ne;
                                                cout << acc_max << "] -- negative displacement \n";
                                        }
                                }
                                acc_meanm[qi] = (acc_min + acc_max)/2.;
                                acc_sigmam[qi] = acc_meanm[qi] - acc_min;
				if(!h1D_0) return; // in case of troubles
				TCanvas canh1D_0("misaligned_acceptance","");
				canh1D_0.cd();
				h1D_0->SetTitle(title);
				h1D_0->Draw();
				h1D_0->GetXaxis()->SetTitle("GeV");
				h1D_0->GetYaxis()->SetTitle("Acceptance");
				canh1D_0.Print(filename);
				delete h1D_0;
			}

			p1.resetPath();

			p1.computePath(beam1,1);
			if(p1.stopped(beam1) || !testtag) {
				if(p1.getStoppingElement()->getName()>=rp_name || !testtag) {
					cout<<"\t d- check"<<endl;;
					p1.propagate(rp1_s); x1 = p1.getX(); y1 = p1.getY();

                                        if(corrected==2) { // correction based on the calibration
                                                TF1* calib_1 = misalignment_calibration_plot(-2*displacement,1,beam1->getElement(i)->getName().c_str());
                                                if(!calib_1){return;}
                                                e1m[qi] = calib_1->Eval(p1.getX()/1000.)/energy;
                                                delete calib_1;
                                        } else if(corrected==1) {
						H_BeamParticle p0; // * reference proton, no energy loss * /
                                                p0.setPosition(PX,PY,-CRANG,0,0);
						p0.computePath(beam1,1);
						p0.propagate(rp1_s);x01 = p0.getX();y01 = p0.getY();
						p1.propagate(rp2_s); x2 = p1.getX(); y2 = p1.getY();
						p0.propagate(rp2_s);x02 = p0.getX();y02 = p0.getY();
						recp.setPositions(x1-x01,y1-y01,x2-x02,y2-y02);
						e1m[qi] = recp.getE(method)/energy;
					} else {
						recp.setPositions(x1,y1,x2,y2);
						e1m[qi] = recp.getE(method)/energy;
					}
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
	
	// acceptance graph
	const int N = nquadrupoles_1;
	float es[N], sp[N], sm[N];
	for (int i=0; i<N; i++) { 
		es[i] = 1.;
		sp[i] = s1[i]-es[i];
		sm[i] = s1[i]+es[i];
	}

	// drawing
	char mytitle[500];
	if (acceptance) {
		TCanvas* cAcc = new TCanvas("cAcc","Changes in acceptance");
		cAcc->SetGrid();
		cAcc->cd();
		TGraphErrors * tp = new TGraphErrors(N,sp,acc_meanp,es,acc_sigmap);
		TGraphErrors * tm = new TGraphErrors(N,sm,acc_meanm,es,acc_sigmam);
		tp->Draw("APE2");
		sprintf(mytitle,"Impact of quadrupole displacement (%.0f #mum) on RP acceptance ",displacement*1000000.);
		tp->SetTitle(mytitle);
		tp->GetXaxis()->SetTitle("s (m)");
		tp->GetYaxis()->SetTitle("RP acceptance (> 90%) [GeV]");
		tm->Draw("PE2same");
		tm->SetFillColor(kRed);
		TLegend * lef = new TLegend(0.68,0.10,0.99,0.30);
		lef->AddEntry(tp,"positive displacement","f");
		lef->AddEntry(tm,"negative displacement","f");
		lef->Draw();
	        lef->SetBorderSize(1);
        	lef->SetFillColor(0);
	}
	
	TLegend* leg=new TLegend(0.68,0.75,0.99,0.88);
	TCanvas* c1 = new TCanvas("misalignment","Misalignment of quadrupoles",1);
	c1->SetGrid();
	c1->cd();
	TGraph* g1p = new TGraph(nquadrupoles_1,s1,e1p);
	g1p->SetLineStyle(1);
	g1p->SetLineColor(1);
	g1p->SetMarkerColor(1);
	g1p->SetMarkerStyle(22);

	leg->AddEntry(g1p,"positive displacement","pl");
	TGraph* g1m = new TGraph(nquadrupoles_1,s1,e1m);
	g1m->SetLineStyle(2);
	g1m->SetLineColor(2);
	g1m->SetMarkerColor(2);
	g1m->SetMarkerStyle(23);

	leg->AddEntry(g1m,"negative displacement","pl");
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


/// Outputs a graph for the misalignment study : see misalignment_multigraph
TGraph * misalignment_quad_graph(float energy = 100., float displacement = 0.0005, bool corrected = true, int markerstyle = 22, int linestyle = 1, int color = 1) {

	extern bool relative_energy;
	relative_energy = true;
	if(!relative_energy) {
		cout << "You should be in relative energy" << endl;
		return 0;
	}

	float x1,y1,x2=0,y2=0;
	float x01,y01,x02,y02;
	
	H_BeamLine* beam1 = new H_BeamLine( 1,rp2_s+0.1);
	beam1->fill("data/LHCB1IR5_v6.500.tfs");
	H_RecRPObject recp1(rp1_s,rp2_s,beam1);
	H_RomanPot * rp = new H_RomanPot("RP",rp1_s,rp_x);
	beam1->add(rp);

	const int nelements = beam1->getNumberOfElements();

	int qi = 0;
	int nquadrupoles = 0;

	for(int i = 0; i < nelements; i++) {
		if(beam1->getElement(i)->getType() == VQUADRUPOLE || beam1->getElement(i)->getType() == HQUADRUPOLE) {
			nquadrupoles++;
		}
	}

	float s1[nquadrupoles];
	float e1[nquadrupoles];

	for(int i = 0; i < nelements; i++) {
		if(beam1->getElement(i)->getType() == VQUADRUPOLE || beam1->getElement(i)->getType() == HQUADRUPOLE) {
			s1[qi] = beam1->getElement(i)->getS();
			H_BeamParticle p1; /* test proton */
			H_BeamParticle p0; /* reference proton, no energy loss */
			p1.setE(7000. - energy);
			if (shift) beam1->alignElement(beam1->getElement(i)->getName(), displacement,0);
                        else    beam1->tiltElement(beam1->getElement(i)->getName(), displacement,0);
			
			p1.computePath(beam1,1);
			p0.computePath(beam1,1);
			p1.propagate(rp1_s);
			p0.propagate(rp1_s);
			x1 = p1.getX();
			x01 = p0.getX();
			y1 = p1.getY();
			y01 = p0.getY();
			p1.propagate(rp2_s);
			p0.propagate(rp2_s);
			x2 = p1.getX();
			x02 = p0.getX();
			y2 = p1.getY();
			y02 = p0.getY();
			if(corrected) {
				recp1.setPositions(x1-x01,y1-y01,x2-x02,y2-y02);
			} else {
				recp1.setPositions(x1,y1,x2,y2);
			}
			e1[qi] = recp1.getE(TM)/energy;
			
			if (shift) beam1->alignElement(beam1->getElement(i)->getName(),-displacement,0);
                        else    beam1->tiltElement(beam1->getElement(i)->getName(), displacement,0);
			qi++;
		}
	}

	char name[50];
	char title[150];
	if(corrected) { 
		sprintf(name,"E = %.0fGeV, corrected",energy);
	} else {
		sprintf(name,"E = %.0fGeV, uncorrected",energy);
	}
	sprintf(title,"Effect of %.0f #mum displacement of quadrupoles on photon energy reconstruction (E = %.0f GeV)",displacement*1000000,energy);
	TGraph* gr = new TGraph(nquadrupoles,s1,e1);
	gr->SetMarkerStyle(markerstyle);
	gr->SetMarkerColor(color);
	gr->SetLineStyle(linestyle);
	gr->SetLineColor(color);
	gr->SetName(name);
	gr->SetTitle(title);
	return gr;
}
/// Impact on E reconstruction of the misalignment of each quadrupole (separately), for a RP at 420, with E scan. 
void misalignment_multigraph(float emin = 20, float emax = 100, int ne = 3, int docorrected = 1, bool donegative = false, float displacement = 0.0005) {
	/// @param emin : minimal proton energy loss 
	/// @param emax : maximal proton energy loss 
	/// @param ne : number of steps in the scan in energy
        /// @param docorrected : 0 = only uncorrected, 1 = only corrected, 2 = both
        /// @param donegative : plots the results also for the opposite displacement.
	/// @param displacement : amplitude of the misalignment (\f$\mu m\f$ for shifts)
	char title[150];
	sprintf(title,"Effect of a transverse displacement of quadrupoles of %.0f #mum on energy reconstruction",1000000.*displacement);
	TMultiGraph* mg = new TMultiGraph("mg",title);
	TLegend* lg=new TLegend(0.75,0.7,0.95,0.8);
	TGraph* tg = 0;
	float ei = emin;
	for(int i = 0; i < ne; i++) {
		ei = (ne==1)?emin:emin + i*(emax-emin)/((float)ne-1);
		if(docorrected==0 || docorrected==2) { 
			tg = misalignment_quad_graph(ei,displacement,0,22,1,i+1);
			mg->Add(tg,"lp");
			lg->AddEntry(tg,tg->GetName(),"lp");
			if(donegative) {
				tg = misalignment_quad_graph(ei,-displacement,0,26,1,i+1);
				mg->Add(tg,"lp");
			}
		}
		if(docorrected==1 || docorrected==2) {
			tg = misalignment_quad_graph(ei,displacement,1,22,2,i+1);
			mg->Add(tg,"lp");
			lg->AddEntry(tg,tg->GetName(),"lp");
			if(donegative) {
				tg = misalignment_quad_graph(ei,-displacement,1,26,2,i+1);
				mg->Add(tg,"lp");
			}
		}
	}
	mg->Draw("a");
	lg->Draw();
	return;
}


/// Misalignment impact on the reconstruction of the Higgs mass (exclusive photoproduction).
void misalignment_compare_higgs(char* name = "\"MQM.9R5.B1\"", float dx = 0.0001, float dy = 0.) {
	/// @param name : displaced element name
	/// @param dx : displacement amplitude, in horizontal plane (shift in \f$ \mu m \f$)
	/// @param dy : displacement amplitude, in vertical plane (shift in \f$ \mu m \f$)

	char* calfile = "data/lpair_mumu_2gev.root";

	cout<<"Computing calibration for beam 1 ..."<<endl;
	TF1* calib_1 = misalignment_calibration_plot(dx,1,name,1,-1,calfile);
	cout<<"Computing calibration for beam 2 ..."<<endl;
	TF1* calib_2 = misalignment_calibration_plot(0 ,1,"\"MQXA.1L5\"",2,-1,calfile);
	cout<<"Calibration ok, reading Higgs file ..."<<endl;

	TFile* gammafile = new TFile("data/aah_gammae.root");
	TTree* gammatree = (TTree*) gammafile->Get("newtree");
	float gamma1_e, gamma2_e;
	gammatree->SetBranchAddress("gamma1_e",&gamma1_e);
	gammatree->SetBranchAddress("gamma2_e",&gamma2_e);

	extern bool relative_energy;
	relative_energy = true;
	if(!relative_energy) {
		cout << "You should be in relative energy" << endl;
		return;
	}

    H_BeamLine* beam1 = new H_BeamLine( 1,length);
	H_BeamLine* beam2 = new H_BeamLine(-1,length);
	H_BeamLine* beam3 = new H_BeamLine( 1,length);
	H_BeamLine* beam4 = new H_BeamLine(-1,length);
    beam1->fill("data/LHCB1IR5_v6.500.tfs");
	beam2->fill("data/LHCB2IR5_v6.500.tfs");
	beam3->fill("data/LHCB1IR5_v6.500.tfs");
	beam4->fill("data/LHCB2IR5_v6.500.tfs");
    H_RecRPObject recp1_2(rp1_s,rp2_s,beam1);
	H_RecRPObject recp2_2(rp1_s,rp2_s,beam2);
/*	H_RecRPObject recp1_1(220,225,beam1);
	H_RecRPObject recp2_1(220,225,beam2);

	H_RomanPot * rp1 = new H_RomanPot("RP1",220.,2000.);
	beam1->add(rp1);
	beam2->add(rp1);
	beam3->add(rp1);
	beam4->add(rp1);
	beam3->showElements();
*/
	H_RomanPot * rp2 = new H_RomanPot("RP2",rp1_s,rp_x);
	beam1->add(rp2);
	beam2->add(rp2);
	beam3->add(rp2);
	beam4->add(rp2);
	           
	if (shift) beam3->alignElement(name,dx,dy);
        else    beam3->tiltElement(name,dx,dy);

	float x1_1,x1_2,x2_1,x2_2,y1_1,y1_2,y2_1,y2_2;
	float x3_1,x3_2,x4_1,x4_2,y3_1,y3_2,y4_1,y4_2;
	float x03_1,x03_2,x04_1,x04_2,y03_1,y03_2,y04_1,y04_2;
	float rece_1=0, rece_2=0, rece_3=0, rece_4=0, rece_3corr=0, rece_4corr=0, rece_3cor2=0, rece_4cor2=0;
	bool p1hit, p2hit, p3hit, p4hit;

	float nacc_base = 0;
	float nacc_misa = 0;

	TH1F* higgs_theo = new TH1F("higgs_theo","MC higgs mass",50,95,135);
        higgs_theo->SetLineColor(kBlack);
        higgs_theo->SetLineWidth(2);
        higgs_theo->SetFillStyle(4000);
	TH1F* higgs_reco = new TH1F("higgs_reco","reco higgs   ",50,95,135);
	higgs_reco->SetLineColor(kRed);
	higgs_reco->SetLineStyle(9);
	higgs_reco->SetLineWidth(2);
	TH1F* higgs_rec2 = new TH1F("higgs_rec2","reco higgs 2 ",50,95,135);
	higgs_rec2->SetLineColor(kGreen);
	higgs_rec2->SetLineWidth(2);
	TH1F* higgs_rec3 = new TH1F("higgs_rec3","reco higgs 3 ",50,95,135);
	higgs_rec3->SetLineColor(kYellow);
	higgs_rec3->SetLineWidth(2);
	higgs_rec3->SetFillStyle(3004);
	higgs_rec3->SetFillColor(kYellow);
	TH1F* higgs_rec4 = new TH1F("higgs_rec4","reco higgs 4 ",50,95,135);
	higgs_rec4->SetLineColor(38);
        higgs_rec4->SetFillColor(38);

   int nentries = gammatree->GetEntriesFast();
   for (int jentry=0; jentry<nentries;jentry++) {
	   gammatree->GetEntry(jentry);

	  // start of loop
	
	    // MC higgs mass :  
	  	higgs_theo->Fill(2*sqrt(gamma1_e*gamma2_e));

	  	// Hector higgs mass : 
	  	H_BeamParticle p1;
	 	H_BeamParticle p2;
		p1.smearPos();
		p1.smearAng();
		p2.smearPos();
		p2.smearAng();
	  	p1.setE(7000. - gamma1_e);
	  	p2.setE(7000. - gamma2_e);
	  	p1.computePath(beam1,1);
	  	p2.computePath(beam2,1);

		p1hit = false;
		p2hit = false;
		p3hit = false;
		p4hit = false;

		if(p1.stopped(beam1)) {
			/*if(p1.getStoppingElement()->getName()>="RP1") {
				p1.propagate(220);
				x1_1 = p1.getX();
				y1_1 = p1.getY();
				p1.propagate(225);
				x1_2 = p1.getX();
				y1_2 = p1.getY();
				//recp1_1.setPositions(x1_1,y1_1,x1_2,y1_2);
				//rece_1 = recp1_1.getE(TM);
				//p1hit = true;
			}*/
			if(p1.getStoppingElement()->getName()>="RP2") {
				p1.propagate(rp1_s);
				x1_1 = p1.getX();
				y1_1 = p1.getY();
				p1.propagate(rp2_s);
				x1_2 = p1.getX();
				y1_2 = p1.getY();
				recp1_2.setPositions(x1_1,y1_1,x1_2,y1_2);
				rece_1 = recp1_2.getE(TM);
				p1hit = true;
			}
		}

		if(p2.stopped(beam2)) {
			/*if(p2.getStoppingElement()->getName()>="RP1") {
				p2.propagate(220);
				x2_1 = p2.getX();
				y2_1 = p2.getY();
				p2.propagate(225);
				x2_2 = p2.getX();
				y2_2 = p2.getY();
	  			//recp2_1.setPositions(x2_1,y2_1,x2_2,y2_2);
				//rece_2 = recp2_1.getE(TM);
				//p2hit = true;
			}*/
			if(p2.getStoppingElement()->getName()>="RP2") {
				p2.propagate(rp1_s);
				x2_1 = p2.getX();
				y2_1 = p2.getY();
				p2.propagate(rp2_s);
				x2_2 = p2.getX();
				y2_2 = p2.getY();
				recp2_2.setPositions(x2_1,y2_1,x2_2,y2_2);
				rece_2 = recp2_2.getE(TM);
				p2hit = true;
			}
		}

		if(p1hit && p2hit) {
			higgs_reco->Fill(2*sqrt(rece_1*rece_2));
			nacc_base++;
		}

	    // Hector Higgs mass with misalignement :

		// reference protons :

		H_BeamParticle p03;
		H_BeamParticle p04;
		p03.computePath(beam3,1);
		p04.computePath(beam4,1);
		p03.propagate(rp1_s);
		p04.propagate(rp1_s);

		x03_1 = p03.getX();
		x04_1 = p04.getX();
		y03_1 = p03.getY();
		y04_1 = p04.getY();

		p03.propagate(rp2_s);
		p04.propagate(rp2_s);

		x03_2 = p03.getX();
		x04_2 = p04.getX();
		y03_2 = p03.getY();
		y04_2 = p04.getY();

		// MC protons :

		H_BeamParticle p3;
	 	H_BeamParticle p4;
		p3.smearPos();
		p3.smearAng();
		p4.smearPos();
		p4.smearAng();
		p3.setE(7000. - gamma1_e);
		p4.setE(7000. - gamma2_e);
		p3.computePath(beam3,1);
		p4.computePath(beam4,1);

		if(p3.stopped(beam3)) {
			if(p3.getStoppingElement()->getName()>="RP1") {
				p3.propagate(220);
				p03.propagate(220);
				x3_1 = p3.getX();
				x03_1 = p03.getX();
				y3_1 = p4.getY();
				y03_1 = p04.getY();
				p3.propagate(225);
				p03.propagate(225);
				x3_2 = p3.getX();
				x03_2 = p03.getX();
				y3_2 = p3.getY();
				y03_2 = p03.getY();
//				recp1_1.setPositions(x3_1,y3_1,x3_2,y3_2);
//				rece_3 = recp1_1.getE(TM);
//				recp1_1.setPositions(x3_1 - x03_1,y3_1 - y03_1,x3_2 - x03_2,y3_2 - y03_2);
//				rece_3corr = recp1_1.getE(TM);
//				p3hit = true;
			}
			if(p3.getStoppingElement()->getName()>="RP2") {
				p3.propagate(rp1_s);
				p03.propagate(rp1_s);
				x3_1 = p3.getX();
				rece_3cor2 = calib_1->Eval(p3.getX()/1000.);
				x03_1 = p03.getX();
				y3_1 = p3.getY();
				y03_1 = p03.getY();
				p3.propagate(rp2_s);
				p03.propagate(rp2_s);
				x3_2 = p3.getX();
				x03_2 = p03.getX();
				y3_2 = p3.getY();
				y03_2 = p03.getY();
				recp1_2.setPositions(x3_1,y3_1,x3_2,y3_2);
				rece_3 = recp1_2.getE(TM);
				recp1_2.setPositions(x3_1 - x03_1,y3_1 - y03_1,x3_2 - x03_2,y3_2 - y03_2);
				rece_3corr = recp1_2.getE(TM);
				p3hit = true;
			}
		}

		if(p4.stopped(beam4)) {
			if(p4.getStoppingElement()->getName()>="RP1") {
				p4.propagate(220);
				p04.propagate(220);
				x4_1 = p4.getX();
				x04_1 = p04.getX();
				y4_1 = p4.getY();
				y04_1 = p04.getY();
				p4.propagate(225);
				p04.propagate(225);
				x4_2 = p4.getX();
				x04_2 = p04.getX();
				y4_2 = p4.getY();
				y04_2 = p04.getY();
//				recp2_1.setPositions(x4_1,y4_1,x4_2,y4_2);
//				rece_4 = recp2_1.getE(TM);
//				recp2_1.setPositions(x4_1 - x04_1,y4_1 - y04_1,x4_2 - x04_2,y4_2 - y04_2);
//				rece_4corr = recp2_1.getE(TM);
//				p4hit = true;
			}
			if(p4.getStoppingElement()->getName()>="RP2") {
				p4.propagate(rp1_s);
				p04.propagate(rp1_s);
				x4_1 = p4.getX();
				rece_4cor2 = calib_2->Eval(p4.getX()/1000.);
				x04_1 = p04.getX();
				y4_1 = p4.getY();
				y04_1 = p04.getY();
				p4.propagate(rp2_s);
				p04.propagate(rp2_s);
				x4_2 = p4.getX();
				x04_2 = p04.getX();
				y4_2 = p4.getY();
				y04_2 = p04.getY();
				recp2_2.setPositions(x4_1,y4_1,x4_2,y4_2);
				rece_4 = recp2_2.getE(TM);
				recp2_2.setPositions(x4_1 - x04_1,y4_1 - y04_1,x4_2 - x04_2,y4_2 - y04_2);
				rece_4corr = recp2_2.getE(TM);
				p4hit = true;
			}
		}

		if(p3hit && p4hit) { 
			higgs_rec2->Fill(2*sqrt(rece_3*rece_4));
			higgs_rec3->Fill(2*sqrt(rece_3corr*rece_4corr));
			higgs_rec4->Fill(2*sqrt(rece_3cor2*rece_4cor2));
			nacc_misa++;
		}
	  // end of loop
   }

	TCanvas* can = new TCanvas("higgs","higgs reconstruction",1);
	THStack * stack = new THStack();
	stack->Add(higgs_rec4);
	stack->Add(higgs_reco);
	stack->Add(higgs_rec2);
	stack->Add(higgs_theo);
	stack->Add(higgs_rec3);
	stack->Draw("nostack");
	stack->GetXaxis()->SetTitle("Reconstructed Higgs Mass (GeV)");
	stack->GetYaxis()->SetTitle("Events");
	char title[500];
        sprintf(title,"Misalignment impact on Higgs mass reconstruction");
	stack->SetTitle(title);	

        TLegend* leg = new TLegend(0.53,0.61,0.999,0.88);
        leg->SetFillColor(0);
	leg->SetBorderSize(0);
	leg->SetMargin(0.12);
	sprintf(title,"E_{rec} with Hector -- %s shifted by %.1fmm",name, dx*1000.);
        leg->SetHeader(title);
        leg->AddEntry(higgs_theo,"MC particles","l");

	sprintf(title,"No misalignment (%.1f #pm %.1f GeV)",higgs_reco->GetMean(),higgs_reco->GetRMS());
        leg->AddEntry(higgs_reco,title,"l");

	sprintf(title,"Misalignment (%.1f #pm %.1f GeV)",higgs_rec2->GetMean(),higgs_rec2->GetRMS());
        leg->AddEntry(higgs_rec2,title,"l");

	sprintf(title,"Misalign + beam pos corr. (%.1f #pm %.1f GeV)",higgs_rec3->GetMean(),higgs_rec3->GetRMS());
        leg->AddEntry(higgs_rec3,title,"fl");

	sprintf(title,"Misalign + calibration (%.1f #pm %.1f GeV)",higgs_rec4->GetMean(),higgs_rec4->GetRMS());
        leg->AddEntry(higgs_rec4,title,"f");


	leg->Draw();
	can->SetLogy();

	cout<<"--------------------------"<<endl;
	cout<<"Hector reconstruction  : "<<endl;
	cout<<"\t Higgs mass : "<< higgs_reco->GetMean() <<"\t GeV"<<endl;
	cout<<"\t RMS        : "<< higgs_reco->GetRMS()  <<"\t GeV"<<endl;
	cout<<"\t Efficiency : "<< (float)nacc_base/(float)nentries <<endl;
	cout<<"Misalignement included : "<<endl;
	cout<<"\t Higgs mass : "<< higgs_rec2->GetMean() <<"\t GeV"<<endl;
	cout<<"\t RMS        : "<< higgs_rec2->GetRMS()  <<"\t GeV"<<endl;
	cout<<"\t Efficiency : "<< (float)nacc_misa/(float)nentries <<endl;
	cout<<"Misalignement and beam position correction included : "<<endl;
	cout<<"\t Higgs mass : "<< higgs_rec3->GetMean() <<"\t GeV"<<endl;
	cout<<"\t RMS        : "<< higgs_rec3->GetRMS()  <<"\t GeV"<<endl;
	cout<<"Misalignement and physics-based calibration included : "<<endl;
	cout<<"\t Higgs mass : "<< higgs_rec4->GetMean() <<"\t GeV"<<endl;
	cout<<"\t RMS        : "<< higgs_rec4->GetRMS()  <<"\t GeV"<<endl;

   return;
}

int main() {
	misalignment_quad_bias(100.,5E-4,1,TM,false);
	return 0;
}