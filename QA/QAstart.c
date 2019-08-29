{
	gROOT->ProcessLine(".L particles.C");
	gROOT->ProcessLine("particles t");
	gROOT->ProcessLine("t.Loop()");
	gROOT->ProcessLine(".q");
}

