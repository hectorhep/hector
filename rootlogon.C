{
	ifstream linuxdistrib("/etc/issue");
	if (linuxdistrib.good()) { // be careful if using Fedora ! see your SELinux TroubleShooter
		string distrib;
		linuxdistrib >> distrib;
		if (distrib.find("edora")) cout << "Using Hector on Fedora ? Check README and your SELinux TroubleShooter about libHector_routines" << endl;
	}

	gSystem->Load("lib/libHector.so");
	gSystem->Load("lib/libHector_routines"); // be careful if using Fedora ! see your SELinux TroubleShooter

	gROOT->SetStyle("Plain");
	gSystem->AddIncludePath("-I\"./include\"");
	gSystem->AddIncludePath("-I\"./routines\"");
	ifstream fversion("VERSION");
	string ver;
	fversion >> ver;
	fversion.close();
	cout << "Ready for Hector (v"<<ver<<") -- enjoy !" << endl;
}
