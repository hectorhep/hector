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

{
	ifstream linuxdistrib("/etc/issue");
	if (linuxdistrib.good()) { // be careful if using Fedora ! see your SELinux TroubleShooter
		string distrib;
		linuxdistrib >> distrib;
		if (distrib.find("edora")) cout << "Using Hector on Fedora ? Check README and your SELinux TroubleShooter about libHector_routines" << endl;
	}

	gSystem->Load("lib/libHector");
	//gSystem->Load("lib/libHector_routines"); // be careful if using Fedora ! see your SELinux TroubleShooter

	gROOT->SetStyle("Plain");
	gSystem->AddIncludePath("-I\"./include\"");
	gSystem->AddIncludePath("-I\"./routines\"");
	ifstream fversion("VERSION");
	string ver;
	fversion >> ver;
	fversion.close();
	cout << "Ready for Hector (v"<<ver<<") -- enjoy !" << endl;
}
