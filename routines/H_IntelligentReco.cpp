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

#include <iostream>
#include "TF1.h"
#include "TH1F.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TRandom.h"
#include "TMultiGraph.h"

#include "H_BeamLine.h"
#include "H_BeamParticle.h"
#include "H_RecRPObject.h"

#define MDIM 6
#define MEGA 1000000.

/* The purpose here is to offer an intelligent alternative to all the naive reco methods used.
 * In order to achieve this goal, we will cease neglecting parts of the information
 * as a first try, we will only use the horizontal plane
 * after that, the vertical plane will be used too in order to obtain the final method
 * which will be general and, hopefully, ultimate.
 *
 * The X (horizontal) system is described as follows : 
 * 		X_1 = f_1(E) x* + g_1(E) x'* + D_1(E) E
 * 		X_2 = f_2(E) x* + g_2(E) x'* + D_2(E) E
 *		Y_1 = k_1(E) y* + l_1(E) y'*
 *		Y_2 = k_2(E) y* + l_2(E) y'*
 *
 * with X_i from detector, f_i, g_i and D_i from simulations
 * the latter are in fact polynomial fits
 *
 * 1) get x*(E) and x'*(E) : 
 * x (E) = (X_1*g_2 - X_2*g_1 + D_2*g_1*E - D_1*g_2*E) / (f_1*g_2 - f_2*g_1) 
 * x'(E) = (X_1*f_2 - X_2*f_1 + D_2*f_1*E - D_1*f_2*E) / (g_1*f_2 - g_2*f_1)
 * y (E) = Y_1*k_2 - Y_2*k_1 / (k_2*l_1 - k_1*l_2)
 * y'(E) = Y_1*l_2 - Y_2*l_1 / (l_2*k_1 - l_1*k_2)
 *
 * 2) get the likelihood of this solution L(E) = P(x(E))*P(x'(E))*P(E)
 * P being the probability density.
 * for x*,x*' , this is a gaussian centered on 0, of know width
 * for E, this is the EPA spectrum (oh, frak !)
 *
 * 3) take the derivative of this Likelyhood to get the min. 
 * One possible (easy) way to do this is to use the usual -2*log(L)
 * this Q' = -Q/2 = (log(P(x)) + log(P(x')) + log(P(E)))
 * and dQ'/dE = (1/P(x))*dPx/dE + ...
 *
 * 3') we could get the minimum from the minimalizing root function without using derivatives
 * the advantage of this is that we use L directly, and normaization doesn't matter anymore
 *
 * This minimum gives the most probable solution, which is our best pick
 *
 * Quicker and yet performant solution : 
 * 1) compute x(E)
 * 2) minimize x²(E) to get most probable E
 * 3) compute x'(E) from found E
 *
 * Add Y the same way : 
 * 1) compute y(E)
 * 2) minimize x²(E) + y²(E) to get the best pick (can it improve energy reco ?)
 * 3) compute x' and y' from this E
*/
/*
double intelligentreco_simple(double ene = 50., double pt = 0., double detres = 0, double vfdet1 = 420, double vfdet2 = 428, double emin = 20, double emax = 120, bool eorpt = 1, string filename = "data/LHCB1IR5_v6.500.tfs", int side = 1); 

void intelligentreco_graphs(bool eorpt = 1) {
	vfdet = 220, 420
	 * detector resolution = 0,5,30 µm
	 * pt = 0 -> 1
	 * E = 20 -> 120 / 120 -> 1200
	 

	// points : 
	const int N = 5;

	double vfdet1[2] = {420,220};
	double vfdet2[2] = {428,225};
	double enemin[2] = {20,120};
	double enemax[2] = {120,1200};
	double enmean[2] = {50,400};
	double detres[3] = {0,5,30};
	double ptmin     = 0.1;
	double ptmax     = 1.;
	double ptmean    = 0.5;

	double ee[N];
	double pt[N];
	double ys[N];
	double xs[N];

	TCanvas* rescan = new TCanvas("rescan","pt resolutions",1);
	rescan->Divide(2,2);

	char title[100];
	char addt[50];
	char prevt[20];

	for(int i=0;i<2;i++) { 							// i = loop on vfdet 
		for(int j=0;j<2;j++) { 						// j = Pt(E), Pt(Pt) graphs
			sprintf(prevt,(eorpt?"Pt":"Energy"));
			sprintf(addt,(!j?"Pt = %.1f GeV":"E = %.1f GeV"),(!j?ptmean:enmean[i]));
			sprintf(title,"%s reconstruction resolution @ %.1fm, %s",prevt,vfdet1[i],addt);
			TMultiGraph* mg = new TMultiGraph("mg",title); 
			for(int k=0;k<3;k++) { 					// k = detector resolution loop 
				for(int l=0;l<N;l++) {				// l = points/line 
					ee[l] = j?(enmean[i]):(enemin[i] + l*(enemax[i]-enemin[i])/((double)N-1.));
					pt[l] = j?(ptmin + l*(ptmax-ptmin)/((double)N-1.)):ptmean;
					xs[l] = j?pt[l]:ee[l];
					ys[l] = intelligentreco_simple(ee[l],pt[l],detres[k],vfdet1[i],vfdet2[i],enemin[i],enemax[i],eorpt);
				}
				TGraph* temp = new TGraph(N,xs,ys);
				temp->SetMarkerColor(k+1);
				temp->SetMarkerStyle(k+20);
				mg->Add(temp);
			}
			cout<<"Ending graph "<<j+2*i+1<<endl;
			rescan->cd(j+2*i+1);
			gPad->SetGrid();
			mg->Draw("alp");
			rescan->Modified();
			rescan->Update();
			mg->GetXaxis()->SetTitle(j?"Pt (GeV)":"E (GeV)");
			mg->GetYaxis()->SetTitle(eorpt?"#delta Pt (GeV)":"#delta E (GeV)");
			mg->GetXaxis()->SetTitleSize(0.05);
			mg->GetYaxis()->SetTitleSize(0.05);
			mg->GetXaxis()->SetLabelSize(0.05);
			mg->GetYaxis()->SetLabelSize(0.05);
		}
	}
	return;
}
*/
/*
double intelligentreco_simple(double ene , double pt , double detres , double vfdet1 , double vfdet2 , double emin , double emax , bool eorpt, string filename , int side) { 

	const int N_p = 1000;

	// caution : the following doesn't yet include kickers effects !
	// no crossing angle should be used yet
	extern int kickers_on;
	kickers_on = 0;

	extern bool relative_energy;
	relative_energy = true;

	double xp = 0;
	double yp = 0;

	H_BeamLine* beam1 = new H_BeamLine(side,vfdet1);
	H_BeamLine* beam2 = new H_BeamLine(side,vfdet2);
	H_BeamLine* beam = new H_BeamLine(side,vfdet2+5.);
	beam1->fill(filename);
	beam2->fill(filename);
	beam->fill(filename);

	// Parameters come from FIT
	// The fits should be run only once to set the parameters
	// In the final version, those should be protected members of H_RecRPOBject
	// 		and initialized at creation
	// X parameters
	TF1 f_1("f_1","[0] + [1]*x + [2]*x*x ",emin-20.,emax+20.);
	TF1 f_2("f_2","[0] + [1]*x + [2]*x*x ",emin-20.,emax+20.);
	TF1 g_1("g_1","[0] + [1]*x + [2]*x*x ",emin-20.,emax+20.);
	TF1 g_2("g_2","[0] + [1]*x + [2]*x*x ",emin-20.,emax+20.);
	TF1 d_1("d_1","[0] + [1]*x + [2]*x*x ",emin-20.,emax+20.);
	TF1 d_2("d_2","[0] + [1]*x + [2]*x*x ",emin-20.,emax+20.);

	TF1 h_1("h_1","[0] + [1]*x + [2]*x*x ",emin-20.,emax+20.);
	TF1 h_2("h_2","[0] + [1]*x + [2]*x*x ",emin-20.,emax+20.);

	// Y parameters
	TF1 k_1("k_1","[0] + [1]*x + [2]*x*x ",emin-20.,emax+20.);
	TF1 k_2("k_2","[0] + [1]*x + [2]*x*x ",emin-20.,emax+20.);
	TF1 l_1("l_1","[0] + [1]*x + [2]*x*x ",emin-20.,emax+20.);
	TF1 l_2("l_2","[0] + [1]*x + [2]*x*x ",emin-20.,emax+20.);

	// Fitting :
	// float * el_mat = new float[MDIM*MDIM];
	const int N = 20;
	double e_i[N]; 
	double f_1i[N], f_2i[N], g_1i[N], g_2i[N], d_1i[N], d_2i[N], h_1i[N], h_2i[N];
	double k_1i[N], k_2i[N], l_1i[N], l_2i[N];
	for(int i = 0; i < N; i++) {
		e_i[i] = emin + i * (emax - emin)/((double)N-1);
		const float* el_mat1 = (beam1->getBeamMatrix(e_i[i],MP,QP)->GetMatrixArray());
		f_1i[i] = el_mat1[0];
		g_1i[i] = el_mat1[1*MDIM];
		d_1i[i] = MEGA*el_mat1[4*MDIM];
		h_1i[i] = el_mat1[5*MDIM];
		k_1i[i] = el_mat1[2*MDIM+2];
		l_1i[i] = el_mat1[3*MDIM+2];
		const float* el_mat2 = (beam2->getBeamMatrix(e_i[i],MP,QP)->GetMatrixArray());
		f_2i[i] = el_mat2[0];
       	g_2i[i] = el_mat2[1*MDIM];
       	d_2i[i] = MEGA*el_mat2[4*MDIM];
		h_1i[i] = el_mat2[5*MDIM];
		k_2i[i] = el_mat2[2*MDIM+2];
		l_2i[i] = el_mat2[3*MDIM+2];
	}

	TGraph gf_1(N,e_i,f_1i);
	TGraph gg_1(N,e_i,g_1i);
	TGraph gd_1(N,e_i,d_1i);
	TGraph gh_1(N,e_i,h_1i);
	TGraph gf_2(N,e_i,f_2i);
	TGraph gg_2(N,e_i,g_2i);
	TGraph gd_2(N,e_i,d_2i);
	TGraph gh_2(N,e_i,h_2i);
	TGraph gk_1(N,e_i,k_1i);
	TGraph gl_1(N,e_i,l_1i);
	TGraph gk_2(N,e_i,k_2i);
	TGraph gl_2(N,e_i,l_2i);

	// and fitting...
	gf_1.Fit("f_1","Q");
	gg_1.Fit("g_1","Q");
	gd_1.Fit("d_1","Q");
	gh_1.Fit("h_1","Q");
	gf_2.Fit("f_2","Q");
	gg_2.Fit("g_2","Q");
	gd_2.Fit("d_2","Q");
	gh_2.Fit("h_2","Q");
	gk_1.Fit("k_1","Q");
	gl_1.Fit("l_1","Q");
	gk_2.Fit("k_2","Q");
	gl_2.Fit("l_2","Q");

//	cout<<"end of fits, starting reco"<<endl;

	// testing the method 
	TH1F rec_ee("rec_ee","reconstructed energy  ",50,ene-20.,ene+20.);
	TH1F rec_xp("rec_xp","reconstructed x angle ",50,-200,200);
	TH1F rec_yp("rec_yp","reconstructed y angle ",50,-200,200);
	TH1F rec_pt("rec_pt","reconstructed Pt      ",20,0,2);
	

	double eformax, xp_likely, yp_likely;

	for(int i = 0; i < N_p; i++) {
		H_BeamParticle p;
//		p.setPosition(PX,PY,xp-CRANG,yp,0);
		p.setPosition(0,0,xp,yp,0);
		if(i) {
			p.smearPos(SX/sqrt(2.),SY/sqrt(2.));
			p.smearAng();
			p.smearE();
		}
		p.emitGamma(ene,-pt*pt);
		p.computePath(beam,1);
		p.propagate(vfdet1);
		TF1 par0("par0","[0]",emin-20.,emax+20.);
		par0.SetParameter(0,gRandom->Gaus(-p.getX(),detres));
		TF1 par2("par2","[0]",emin-20.,emax+20.);
		par2.SetParameter(0,gRandom->Gaus(-p.getY(),detres));
		p.propagate(vfdet2);
		TF1 par1("par1","[0]",emin-20.,emax+20.);
		par1.SetParameter(0,gRandom->Gaus(-p.getX(),detres));
		TF1 par3("par3","[0]",emin-20.,emax+20.);
		par3.SetParameter(0,gRandom->Gaus(-p.getY(),detres));

		TF1 xx_E("xx_E","(g_2 * (par0 - d_1 * x) - g_1 * (par1 - d_2 * x))/(f_2 * g_1 - f_1 * g_2)",emin-20.,emax+20.);
		TF1 xp_E("xp_E","(f_2 * (par0 - d_1 * x) - f_1 * (par1 - d_2 * x))/(g_2 * f_1 - g_1 * f_2)",emin-20.,emax+20.);
		TF1 yy_E("yy_E","(par2*l_2 - par3*l_1) / (k_2*l_1 - k_1*l_2)",emin-20.,emax+20.);
		TF1 yp_E("yp_E","(par2*k_2 - par3*k_1) / (l_2*k_1 - l_1*k_2)",emin-20.,emax+20.);
		TF1 p_xy_E("p_xy_E","(-xx_E*xx_E - yy_E*yy_E)",emin-20.,emax+20.);
		
		eformax = p_xy_E.GetMaximumX(emin-20.,emax+20.);
		xp_likely = xp_E.Eval(eformax);
		yp_likely = yp_E.Eval(eformax);

		rec_ee.Fill(eformax);
//		rec_xp->Fill(xp_likely);
//		rec_yp->Fill(yp_likely);
		rec_pt.Fill(sqrt(7000*(7000-eformax)*(xp_likely*xp_likely + yp_likely*yp_likely)/(MEGA*MEGA)));
	}

	cout<<"digest : "<<endl;
	cout<<"\t detector @ "<<vfdet1<<" m resolution : "<<detres<<" µm"<<endl;
	cout<<"\t IP   energy : "<< ene              <<"\t IP   Pt : "<<pt<<endl;
    cout<<"\t reco energy : "<<	rec_ee->GetMean()<<"\t reco Pt : "<<rec_pt->GetMean()<<endl;
    cout<<"\t energy rms  : "<< rec_ee->GetRMS() <<"\t Pt RMS  : "<<rec_pt->GetRMS()<<endl;	

	TCanvas* reco = new TCanvas("reco","IP variables reconstruction");
	reco->Divide(2);
	reco->cd(1);
	rec_ee->Draw();
//	reco->cd(3);
//	rec_xp->Draw();
//	reco->cd(4);
//	rec_yp->Draw();
	reco->cd(2);
	rec_pt->Draw();
	return (eorpt?(rec_pt.GetRMS()):(rec_ee.GetRMS()));
}
*/
void intelligentreco_rpo(double ene , double pt , double vfdet1 , double vfdet2 , string filename , int side) { 

	const int N_p = 1000;

	// caution : the following doesn't yet include kickers effects !
	// no crossing angle should be used yet !!!
	extern int kickers_on;
	kickers_on = 0;

	extern bool relative_energy;
	relative_energy = true;

	H_BeamLine* beam = new H_BeamLine(side,vfdet2+5.);
	beam->fill(filename);

	// creating reco object
	H_RecRPObject rec1(vfdet1,vfdet2,beam);
	rec1.computeERange();
//	rec1.setERange(emin,emax);
	// computing reco constants
	rec1.initialize();

	TH1F* rec_ee = new TH1F("rec_ee","reconstructed energy  ",50,ene-20.,ene+20.);
	TH1F* rec_pt = new TH1F("rec_pt","reconstructed Pt      ",50,0,(pt<1)?2:(pt*10));
	
	for(int i = 0; i < N_p; i++) {
		H_BeamParticle p;
		p.smearPos(SX/sqrt(2.),SY/sqrt(2.));
		p.smearAng();
		p.smearE();
		p.emitGamma(ene,-pt*pt);
		p.computePath(beam,1);
		p.propagate(vfdet1);
		rec1.setPosition_det1(p.getX(),p.getY());
		p.propagate(vfdet2);
		rec1.setPosition_det2(p.getX(),p.getY());
		rec1.computeAll();
		rec_ee->Fill(rec1.getE());
		rec_pt->Fill(rec1.getPt());
	}

	cout<<"digest : "<<endl;
	cout<<"\t detector @ "<<vfdet1<<" and "<<vfdet2<<" m"<<endl;
	cout<<"\t IP   energy : "<< ene              <<"\t IP   Pt : "<<pt<<endl;
	cout<<"\t reco energy : "<<	rec_ee->GetMean()<<"\t reco Pt : "<<rec_pt->GetMean()<<endl;
	cout<<"\t energy rms  : "<< rec_ee->GetRMS() <<"\t Pt RMS  : "<<rec_pt->GetRMS()<<endl;	

	TCanvas* reco = new TCanvas("reco","IP variables reconstruction");
	reco->Divide(2);
	reco->cd(1);
	rec_ee->Draw();
	reco->cd(2);
	rec_pt->Draw();
	return;
}

/*
void intelligentreco_xonly_test(double ene = 50., double xx = 0., double xp = 0., string filename = "data/LHCB1IR5_v6.500.tfs", int side = 1) {

	double emin = 20;
	double emax = 120;

	double vfdet1 = 420;
	double vfdet2 = 428;

	// caution : the following doesn't yet include kickers effects !
	// no crossing angle should be used yet
	extern int kickers_on;
	kickers_on = 0;

	extern bool relative_energy;
	relative_energy = true;

	H_BeamLine* beam1 = new H_BeamLine(side,vfdet1);
	H_BeamLine* beam2 = new H_BeamLine(side,vfdet2);
	H_BeamLine* beam = new H_BeamLine(side,vfdet2+5.);
	beam1->fill(filename);
	beam2->fill(filename);
	beam->fill(filename);

	// Parameters come from FIT
	// The fits should be run only once to set the parameters
	// In the final version, those should be protected members of H_RecRPOBject
	// 		and initialized at creation
	TF1 f_1("f_1","[0] + [1]*x + [2]*x*x + [3]*x*x*x",emin,emax);
	TF1 f_2("f_2","[0] + [1]*x + [2]*x*x + [3]*x*x*x",emin,emax);
	TF1 g_1("g_1","[0] + [1]*x + [2]*x*x + [3]*x*x*x",emin,emax);
	TF1 g_2("g_2","[0] + [1]*x + [2]*x*x + [3]*x*x*x",emin,emax);
	TF1 d_1("d_1","[0] + [1]*x + [2]*x*x + [3]*x*x*x",emin,emax);
	TF1 d_2("d_2","[0] + [1]*x + [2]*x*x + [3]*x*x*x",emin,emax);

	// Fitting :
	// float * el_mat = new float[MDIM*MDIM];
	const int N = 200;
	double e_i[N], f_1i[N], f_2i[N], g_1i[N], g_2i[N], d_1i[N], d_2i[N];
	for(int i = 0; i < N; i++) {
		e_i[i] = emin + i * (emax - emin)/((double)N-1);
		const float* el_mat1 = (beam1->getBeamMatrix(e_i[i],MP,QP)->GetMatrixArray());
		f_1i[i] = el_mat1[0];
		g_1i[i] = el_mat1[1*MDIM];
		d_1i[i] = MEGA*el_mat1[4*MDIM];
		const float* el_mat2 = (beam2->getBeamMatrix(e_i[i],MP,QP)->GetMatrixArray());
		f_2i[i] = el_mat2[0];
        g_2i[i] = el_mat2[1*MDIM];
        d_2i[i] = MEGA*el_mat2[4*MDIM];
	}
	TGraph gf_1(N,e_i,f_1i);
	TGraph gg_1(N,e_i,g_1i);
	TGraph gd_1(N,e_i,d_1i);
	TGraph gf_2(N,e_i,f_2i);
	TGraph gg_2(N,e_i,g_2i);
	TGraph gd_2(N,e_i,d_2i);

	// and fitting...
	gf_1.Fit("f_1","Q");
	gg_1.Fit("g_1","Q");
	gd_1.Fit("d_1","Q");
	gf_2.Fit("f_2","Q");
	gg_2.Fit("g_2","Q");
	gd_2.Fit("d_2","Q");

	cout<<"end of fits, starting reco"<<endl;

	// Solutions for x, x' as a function of E :
	// parameters : [0] = X_1, [1] = X_2. SetParameter should be done for each event
//	TF1 par0("par0","[0]",emin,emax);
//	par0.SetParameter(0,6183.21);
//	TF1 par1("par1","[0]",emin,emax);
//	par1.SetParameter(0,7008.86);
//	TF1 xx_E("xx_E","(g_2 * (par0 - d_1 * x) + g_1 * (par1 - d_2 * x))/(f_2 * g_1 - f_1 * g_2)",emin,emax);
//	TF1 xp_E("xp_E","(f_2 * (par0 - d_1 * x) + f_1 * (par1 - d_2 * x))/(g_2 * f_1 - g_1 * f_2)",emin,emax);

//	TF1 sx("sx","[0]",emin,emax);
//	sx.SetParameter(0,SX);
//	TF1 stx("stx","[0]",emin,emax);
//	stx.SetParameter(0,STX);


	// tests !
	TF1* xx_E2 = new TF1("xx_E2","((g_2 * (par0 - d_1 * x) + g_1 * (par1 - d_2 * x))/(f_2 * g_1 - f_1 * g_2))",emin,emax);
	TF1* xp_E2 = new TF1("xp_E2","((f_2 * (par0 - d_1 * x) + f_1 * (par1 - d_2 * x))/(g_2 * f_1 - g_1 * f_2))",emin,emax);

	TF1* test = new TF1("test","(d_1 * x)",emin,emax);
	TCanvas* test1 = new TCanvas("test1","test1",1);
	test->Draw();

//	TF1* p_xx_E2 = new TF1("p_xx_E2","1/sqrt(2*[1]) * 1/[0] * exp(-xx_E2*xx_E2/(2*[0]*[0]))",emin,emax);
	TF1* p_xx_E2 = new TF1("p_xx_E2","(-xx_E2*xx_E2)/(2*sx*sx)",emin,emax);
//	TF1* p_xp_E2 = new TF1("p_xp_E2","1/sqrt(2*[1]) * 1/[0] * exp(-xx_E2*xx_E2/(2*[0]*[0]))",emin,emax);
	TF1* p_xp_E2 = new TF1("p_xp_E2","(-xp_E2*xp_E2)/(2*stx*stx)",emin,emax);

	TF1* L_E2 = new TF1("L_E2","(p_xx_E2 + p_xp_E2)",emin,emax);

	TCanvas* b1 = new TCanvas("b1","b1",1);
	xx_E2->Draw();
	TCanvas* b2 = new TCanvas("b2","b2",1);
	xp_E2->Draw();
	TCanvas* b3 = new TCanvas("b3","b3",1);
	p_xx_E2->Draw();
	TCanvas* b4 = new TCanvas("b4","b4",1);
	p_xp_E2->Draw();
	TCanvas* b5 = new TCanvas("b5","b5",1);
	L_E2->Draw();
    double eformax2 = L_E2->GetMaximumX(emin,emax);
    cout<<eformax2<<" "<<xx_E2->Eval(eformax2)<<" "<<xp_E2->Eval(eformax2)<<endl;


	// probability densities for x,x' as a function of E :
	// parameters : [0] = sigma
	// caution : xp distribution ignores physics !
//	TF1 p_xx_E("p_xx_E","((-xx_E*xx_E)/(2*sx*sx))",emin,emax);
//	TF1 p_xp_E("p_xp_E","((-xp_E*xp_E)/(2*stx*stx))",emin,emax);

	// likelihood function vs E : 
	// caution : this lacks p_E !!! solutions won't favour low E's
//	TF1* L_E = new TF1("L_E","(p_xx_E + p_xp_E)",emin,emax);

	// testing the method 
	TH1F* rec_ee = new TH1F("rec_ee","reconstructed energy  ",21,ene-20.,ene+20.);
	TH1F* rec_xx = new TH1F("rec_xx","reconstructed position",101,-50,50);
	TH1F* rec_xp = new TH1F("rec_xp","reconstructed angle   ",21,-200,200);
	
	const int N_p = 100;

	double eformax, xx_likely, xp_likely;

	for(int i = 0; i < N_p; i++) {
		H_BeamParticle p;
		p.setPosition(xx,0,xp,0,0);
		p.smearPos(SX/sqrt(2.),SY/sqrt(2.));
		p.smearAng();
		p.emitGamma(ene,0);
		p.computePath(beam,1);
		p.propagate(vfdet1);
		TF1 par0("par0","[0]",emin,emax);
		par0.SetParameter(0,-p.getX());
		p.propagate(vfdet2);
		TF1 par1("par1","[0]",emin,emax);
		par1.SetParameter(0,-p.getX());

		TF1* xx_E = new TF1("xx_E","(g_2 * (par0 - d_1 * x) - g_1 * (par1 - d_2 * x))/(f_2 * g_1 - f_1 * g_2)",emin,emax);
		TF1* xp_E = new TF1("xp_E","(f_2 * (par0 - d_1 * x) - f_1 * (par1 - d_2 * x))/(g_2 * f_1 - g_1 * f_2)",emin,emax);

		TF1 sx("sx","[0]",emin,emax);
		sx.SetParameter(0,SX/sqrt(2.));
		TF1 stx("stx","[0]",emin,emax);
		stx.SetParameter(0,STX);

		TF1 p_xx_E("p_xx_E","((-xx_E*xx_E)/(2*sx*sx))",emin,emax);
//		TF1 p_xp_E("p_xp_E","-2*log(xp_E)",emin,emax);
		TF1 p_xp_E("p_xp_E","((-xp_E*xp_E)/(2*stx*stx))",emin,emax);
		TF1 p_E("p_E","(log(7000.-x) - log(x))",emin,emax);
		TF1* L_E = new TF1("L_E","(p_xx_E)",emin,emax);

		if(!i) {
			TCanvas* lecan = new TCanvas("lecan","optimisation canvas",1);
			lecan->SetGrid();
			xp_E->Draw();
//			L_E->Draw();
		}
		
		eformax = L_E->GetMaximumX(emin,emax);
		xx_likely = xx_E->Eval(eformax);
		xp_likely = xp_E->Eval(eformax);
		cout<<eformax<<" "<<xp_E->Eval(eformax)<<endl;

		rec_ee->Fill(eformax);
		rec_xx->Fill(xx_likely);
		rec_xp->Fill(xp_likely);
//		cout<<eformax<<"\t"<<xx_likely<<"\t"<<xp_likely<<endl;
	}

	TCanvas* reco = new TCanvas("reco","IP variables reconstruction");
	reco->Divide(3,1);
	reco->cd(1);
	rec_ee->Draw();
	reco->cd(2);
	rec_xx->Draw();
	reco->cd(3);
	rec_xp->Draw();
	return;
}
*/
