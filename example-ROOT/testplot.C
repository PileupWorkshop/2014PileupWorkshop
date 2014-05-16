{
    TFile *tF = new TFile("example02.root", "OPEN");
    TTree *tree = (TTree*) tF->Get("Events");
    tree->Show(1);

    // make plot

    TH1D* h = new TH1D("h_ptreco_minus_pttruth", "", 100, -50, 50);
    tree->Draw("JetPt-JetPtTruth>>h_ptreco_minus_pttruth");
    h->SetXTitle("offset [GeV]");

   h->Draw(); 

}
