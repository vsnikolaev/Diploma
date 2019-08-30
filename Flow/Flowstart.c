{
	gROOT->ProcessLine(".L Qvector.h");
	gROOT->ProcessLine(".L Qvector.cxx");
	gROOT->ProcessLine(".L particles.C");
	gROOT->ProcessLine("particles t");
	gROOT->ProcessLine("t.Loop()");
	gROOT->ProcessLine(".q");
}

