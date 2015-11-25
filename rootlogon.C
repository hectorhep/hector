{
	gSystem->Load("lib/libHector");
	gROOT->SetStyle("Plain");
	gSystem->AddIncludePath("-I\"./include\"");
	ifstream fversion("VERSION");
	string ver;
	fversion >> ver;
	fversion.close();
	cout << "Ready for Hector (v"<<ver<<") -- enjoy !" << endl;
}
