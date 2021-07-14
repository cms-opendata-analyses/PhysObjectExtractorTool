void plotBeff(){

  TFile *_file0 = TFile::Open("flavortagefficiencies.root"); // your file name goes here "Hadd2016.root"
  TH1D* num = (TH1D*)_file0->Get("mcweightanalyzer/BEffMed_Nptbins_b");
  TH1D* den = (TH1D*)_file0->Get("mcweightanalyzer/BEff_Dptbins_b");
  num->Divide(num,den,1.0,1.0,"B");

  cout << "B FLAVOR Med 2011" << endl;
  for(int i = 1; i < num->GetNbinsX()+1; i++){printf("pT > %f, eff = %f\n",num->GetBinLowEdge(i),num->GetBinContent(i));}

  TH1D* numC = (TH1D*)_file0->Get("mcweightanalyzer/BEffMed_Nptbins_c");
  TH1D* denC = (TH1D*)_file0->Get("mcweightanalyzer/BEff_Dptbins_c");
  numC->Divide(numC,denC,1.0,1.0,"B");

  cout << "C FLAVOR Med 2011" << endl;
  for(int i = 1; i < numC->GetNbinsX()+1; i++){printf("pT > %f, fake = %f\n",numC->GetBinLowEdge(i),numC->GetBinContent(i));}

  TH1D* numL = (TH1D*)_file0->Get("mcweightanalyzer/BEffMed_Nptbins_udsg");
  TH1D* denL = (TH1D*)_file0->Get("mcweightanalyzer/BEff_Dptbins_udsg");
  numL->Divide(numL,denL,1.0,1.0,"B");

  cout << "L FLAVOR Med 2011" << endl;
  for(int i = 1; i < numL->GetNbinsX()+1; i++){printf("pT > %f, fake = %f\n",numL->GetBinLowEdge(i),numL->GetBinContent(i));}

  TCanvas *c1 = new TCanvas("c1","c1",800,600);
  gStyle->SetOptStat(0);
  num->GetYaxis()->SetRangeUser(0.0,1.0);
  num->GetYaxis()->SetTitle("Efficiency of CSV Medium");
  num->GetXaxis()->SetTitle("AK5 jet p_{T} [GeV]");
  num->SetLineWidth(2);
  num->SetMarkerStyle(20);
  num->SetMarkerColor(kBlack);
  num->SetLineColor(kBlack);
  numC->SetLineWidth(2);
  numC->SetMarkerStyle(20);
  numC->SetMarkerColor(kBlue);
  numC->SetLineColor(kBlue);
  numL->SetLineWidth(2);
  numL->SetMarkerStyle(20);
  numL->SetMarkerColor(kRed);
  numL->SetLineColor(kRed);
  num->Draw();
  numC->Draw("same");
  numL->Draw("same");
  TLegend *leg = new TLegend(0.7,0.7,0.98,0.98);
  leg->AddEntry(num,"b quarks","pl");                                        
  leg->AddEntry(numC,"c quarks","pl");                                       
  leg->AddEntry(numL,"udsg","pl");                                           
  leg->Draw("same");
  c1->SaveAs("CSVMedium_2011.png");
  c1->SaveAs("CSVMedium_2011.pdf");
  c1->SaveAs("CSVMedium_2011.root");


}
