{

// Open the input file
TFile* infile = new TFile("poetoutput.root");

// Get some of the TTree objects
TDirectoryFile* electrons = (TDirectoryFile*)infile->Get("myelectrons");
TTree* tree_el = (TTree*)electrons->Get("Events");

TDirectoryFile* muons = (TDirectoryFile*)infile->Get("mymuons");
TTree* tree_mu = (TTree*)muons->Get("Events");

tree_el->Print();
tree_mu->Print();

// Let's make the output file
TFile fout("analysis_output.root","recreate");

// And we will declare two histograms 
TH1F h1("h1","electron pT (GeV/c)",50,0,150);
TH1F h2("h2","muon pT (GeV/c)",50,0,150);

// Now we need a few variables to hold the data we'll extract from 
// the input ROOT file.

Int_t nelectron; // Necessary to keep track of the number of electrons
Int_t nmuon; // Necessary to keep track of the number of muons

// We'll define these assuming we will not write information for
// more than 16 jets. We'll have to check for this in the code otherwise
// it could crash!
//Float_t el_pt[16];
//Float_t mu_pt[16];
std::vector<float> *el_pt=0;
std::vector<float> *mu_pt=0;

// Assign these variables to specific branch addresses
tree_el->SetBranchAddress("numberelectron",&nelectron);
tree_el->SetBranchAddress("electron_pt",&el_pt);

tree_mu->SetBranchAddress("numbermuon",&nmuon);
tree_mu->SetBranchAddress("muon_pt",&mu_pt);

// Get the number of events in the file
Int_t nevents_el = tree_el->GetEntries();
Int_t nevents_mu = tree_el->GetEntries();

printf("%d events in electron tree\n",nevents_el);
printf("%d events in muon tree\n",nevents_mu);

Int_t nevents = -1;

// Check to make sure there are the same number of events in each TTree
if (nevents_el != nevents_mu) {
    printf("The number of events in the elctron tree is not equal to the number of events in the muon tree!");
    exit;
}
else { nevents = nevents_el; }

for (Int_t i=0;i<nevents;i++) {

    if (i%1000==0) {
        printf("Event %d out of %d",i, nevents);
    }

    // Get the values for the i`th event and fill all our local variables
    // that were assigned to TBranches

    /*-----------------------------------------------------------------------*/
    // First do the electron TTree for the i'th event
    tree_el->GetEntry(i);

    // Print the number of electrons in this event
    //printf("%d\n",nelectron);

    // This is because we declared our variable to only hold 16 values.
    if (nelectron > 16) { nelectron = 16; }

    // Print out the transverse momentum for each electron in this event
    for (Int_t j=0;j<nelectron;j++) {
        //printf("Event %d: electron %d  %f\n",i, j, el_pt->at(j)); // For debugging
        // Fill the histogram with each value of pT
        h1.Fill(el_pt->at(j));
    }

    /*-----------------------------------------------------------------------*/
    // Now do the muon TTree for the i'th event
    tree_mu->GetEntry(i);

    // Print the number of muons in this event
    //printf("%d\n",nmuon);

    // This is because we declared our variable to only hold 16 values.
    if (nmuon > 16) { nmuon = 16; }

    // Print out the transverse momentum for each muon in this event
    for (Int_t j=0;j<nmuon;j++) {
        //printf("Event %d: muon %d     %f\n", i, j,mu_pt->at(j)); // For debugging
        // Fill the histogram with each value of pT
        h2.Fill(mu_pt->at(j));
    }

}
// Done looping over the events.
  // Declare a TCanvas
  TCanvas *c1 = new TCanvas("c1", "Electrons", 800, 400);
  TCanvas *c2 = new TCanvas("c2", "Muons", 800, 400);

  c1->cd(0);
  h1.Draw();
  c1->SaveAs("h_el_pt.png");

  c2->cd(0);
  h2.Draw();
  c2->SaveAs("h_mu_pt.png");

  fout.cd();
  h1.Write();
  h2.Write();
  fout.Close();

}
