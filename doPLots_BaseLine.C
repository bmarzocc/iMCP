
void compareScan(){
  int colorsExt[9] = {kRed+1, kBlue, kPink, kGreen+2, kAzure+7, kOrange+7, kCyan, kWhite, kBlack};
  std::vector<int> colors;
  for(int posVec = 0; posVec<6; ++posVec){
    colors.push_back(colorsExt[posVec]);
  }

  std::vector<std::string> nameMCP;
  nameMCP.push_back("MiB1");
  nameMCP.push_back("MiB2");
  nameMCP.push_back("ScB");
  nameMCP.push_back("Planacon");
  nameMCP.push_back("MiB3");
  nameMCP.push_back("RomaX");

  std::string PC_CH0[4] = {"ON", "ON", "ON", "ON"};
  std::string PC_CH1[4] = {"ON", "OFF", "i", "OFF"};
  std::string PC_CH3[4] = {"OFF", "ON", "OFF", "OFF"};
  std::string PC_CH4[4] = {"OFF", "ON", "ON", "OFF"};
  std::string PC_CH5[4] = {"ON", "ON", "ON", "ON"};



  TCanvas** cBas = new TCanvas*[6];  
  TCanvas** cSig = new TCanvas*[6];  
  
  TH1F* chHistoBase_All[6][4];
  TH1F* chHistoSignal_All[6][4];
  TFile* inF;
  
  for(int iw=0; iw<6; ++iw){
      if(iw == 2) continue;
      cBas[iw] = new TCanvas();
      //      cSig[iw] = new TCanvas();
  }

  for(int iS=0; iS<4; ++iS){
    inF = TFile::Open(Form("Scan%d_outHistos.root",iS+1), "read");
    for(int iw=0; iw<6; ++iw){
      if(iw == 2) continue;

      chHistoBase_All[iw][iS] = (TH1F*)inF->Get(Form("histoBase_All_Ch%d",iw));
      chHistoBase_All[iw][iS]->SetDirectory(0);
      chHistoSignal_All[iw][iS] = (TH1F*)inF->Get(Form("histoSignal_All_Ch%d",iw));
      chHistoSignal_All[iw][iS]->SetDirectory(0);
      //      inF->Close();

      chHistoBase_All[iw][iS]->Rebin(10);
      chHistoBase_All[iw][iS]->Scale(1./chHistoBase_All[iw][iS]->Integral());
      chHistoBase_All[iw][iS]->SetLineColor(colors.at(iw)+2.*iS);
      chHistoBase_All[iw][iS]->SetLineWidth(2);  
      chHistoBase_All[iw][iS]->SetLineStyle(1);  

      chHistoSignal_All[iw][iS]->Rebin(10);
      chHistoSignal_All[iw][iS]->SetLineColor(colors.at(iw));
      chHistoSignal_All[iw][iS]->SetLineWidth(2);  
      chHistoSignal_All[iw][iS]->SetLineStyle(1+3*iS);  

      cBas[iw]->cd();
      if(iS == 0) chHistoBase_All[iw][iS]->DrawNormalized();
      else if(iw == 1 && iS == 2) {cBas[5]->cd(); chHistoBase_All[iw][iS]->Draw("same");}
      else chHistoBase_All[iw][iS]->DrawNormalized("same");

//       cSig[iw]->cd();
//       if(iS == 0) chHistoSignal_All[iw][iS]->DrawNormalized();
//       else chHistoSignal_All[iw][iS]->DrawNormalized("same");
    }
  }

  //  return;


  TLegend** leg = new TLegend[6];
  for (int iw=0;iw<6;iw++) {
    if(iw == 2) continue;
    leg[iw] = new TLegend(0.88,0.65,0.98,0.85);
    leg[iw]->SetFillColor(0);
    leg[iw]->SetTextFont(41);
  }

    leg[0]->AddEntry(chHistoBase_All[0][0],(nameMCP.at(0)+PC_CH0[0]).c_str(),"l");
    leg[0]->AddEntry(chHistoBase_All[0][1],(nameMCP.at(0)+PC_CH0[1]).c_str(),"l");
    leg[0]->AddEntry(chHistoBase_All[0][2],(nameMCP.at(0)+PC_CH0[2]).c_str(),"l");
    leg[0]->AddEntry(chHistoBase_All[0][3],(nameMCP.at(0)+PC_CH0[3]).c_str(),"l");

    leg[1]->AddEntry(chHistoBase_All[1][0],(nameMCP.at(1)+PC_CH1[0]).c_str(),"l");
    leg[1]->AddEntry(chHistoBase_All[1][1],(nameMCP.at(1)+PC_CH1[1]).c_str(),"l");
    leg[1]->AddEntry(chHistoBase_All[1][3],(nameMCP.at(1)+PC_CH1[3]).c_str(),"l");
   
    leg[3]->AddEntry(chHistoBase_All[3][0],(nameMCP.at(3)+PC_CH3[0]).c_str(),"l");
    leg[3]->AddEntry(chHistoBase_All[3][1],(nameMCP.at(3)+PC_CH3[1]).c_str(),"l");
    leg[3]->AddEntry(chHistoBase_All[3][2],(nameMCP.at(3)+PC_CH3[2]).c_str(),"l");
    leg[3]->AddEntry(chHistoBase_All[3][3],(nameMCP.at(3)+PC_CH3[3]).c_str(),"l");

    leg[4]->AddEntry(chHistoBase_All[4][0],(nameMCP.at(4)+PC_CH4[0]).c_str(),"l");
    leg[4]->AddEntry(chHistoBase_All[4][1],(nameMCP.at(4)+PC_CH4[1]).c_str(),"l");
    leg[4]->AddEntry(chHistoBase_All[4][2],(nameMCP.at(4)+PC_CH4[2]).c_str(),"l");
    leg[4]->AddEntry(chHistoBase_All[4][3],(nameMCP.at(4)+PC_CH4[3]).c_str(),"l");

    leg[5]->AddEntry(chHistoBase_All[5][0],(nameMCP.at(5)+PC_CH5[0]).c_str(),"l");
    leg[5]->AddEntry(chHistoBase_All[5][1],(nameMCP.at(5)+PC_CH5[1]).c_str(),"l");
    leg[5]->AddEntry(chHistoBase_All[1][2],(nameMCP.at(5)+PC_CH1[2]).c_str(),"l"); 
    leg[5]->AddEntry(chHistoBase_All[5][2],(nameMCP.at(5)+PC_CH5[2]).c_str(),"l");
    leg[5]->AddEntry(chHistoBase_All[5][3],(nameMCP.at(5)+PC_CH5[3]).c_str(),"l");


  //  return;

  for (int iw=0;iw<6;++iw) {
    if(iw == 2) continue;
    cBas[iw]->cd();
    leg[iw]->Draw("same");
//     cSig[iw]->cd();
//     leg[iw]->Draw("same");
  }

}



void doPlots_BaseLine(std::string fileName){
  TFile* inF = TFile::Open(fileName.c_str(), "read");

  int colorsExt[9] = {kRed+1, kBlue, kPink, kGreen+2, kAzure+7, kOrange+7, kCyan, kWhite, kBlack};
  std::vector<int> colors;
  for(int posVec = 0; posVec<6; ++posVec){
    colors.push_back(colorsExt[posVec]);
  }

  std::vector<std::string> nameMCP;
  nameMCP.push_back("MiB1");
  nameMCP.push_back("MiB2");
  nameMCP.push_back("ScB");
  nameMCP.push_back("Planacon");
  nameMCP.push_back("MiB3");
  nameMCP.push_back("RomaX");

  TGraph* gCont[6];

  TH1F* hBaseLineAll[9];
  char h0[10];
  for (int iw=0;iw<6;iw++){
    sprintf (h0,"BsACh_Ch%d",iw);
    hBaseLineAll[iw] = (TH1F*)inF->Get(h0);
    hBaseLineAll[iw]->SetLineColor(colors.at(iw));
    hBaseLineAll[iw]->SetLineWidth(2);

    gCont[iw] = new TGraph();
    gCont[iw]->SetMarkerColor(colors.at(iw)); 
    gCont[iw]->SetLineColor(colors.at(iw)); 
    gCont[iw]->SetMarkerStyle(21); 
    gCont[iw]->SetMarkerSize(0.6); 
 }

  TLegend* leg = new TLegend(0.88,0.65,0.98,0.85);
  leg->SetFillColor(0);
  leg->SetTextFont(41);
  for (int iw=0;iw<6;iw++) {
    if(iw == 2) continue;
    leg->AddEntry(hBaseLineAll[iw],nameMCP.at(iw).c_str(),"l");
  }

 
  TCanvas* c = new TCanvas();
  hBaseLineAll[0]->Draw();
  for (int iw=1;iw<6;iw++){
    hBaseLineAll[iw]->Draw("same");
  }
  leg->Draw("same");


  for(int iBin=0; iBin<80; ++iBin){
    for (int iw=0;iw<6;iw++){
      if(iw == 2) continue;
      double fake = -1;
      double num = hBaseLineAll[iw]->Integral(hBaseLineAll[iw]->FindBin(iBin), hBaseLineAll[iw]->GetNbinsX());
      double den = hBaseLineAll[iw]->Integral(1, hBaseLineAll[iw]->GetNbinsX());
      if(den != 0) fake = num/den;
      
      gCont[iw]->SetPoint(iBin, iBin, fake*100.);
    }
  }
}

///////////

void doPLots_BaseLine_Envelope(std::string fileName){
  TFile* inF = TFile::Open(fileName.c_str(), "read");

  int colorsExt[9] = {kRed+1, kPink, kBlue, kGreen+2, kAzure+7, kOrange+7, kCyan, kWhite, kBlack};
  std::vector<int> colors;
  for(int posVec = 0; posVec<6; ++posVec){
    colors.push_back(colorsExt[posVec]);
  }

  int nScan = 17; // First
  int nScan = 12; // Last

  std::vector<std::string> nameMCP;
  nameMCP.push_back("MiB1");
  nameMCP.push_back("MiB2");
  nameMCP.push_back("ScB");
  nameMCP.push_back("Planacon");
  nameMCP.push_back("MiB3");
  nameMCP.push_back("RomaX");

  //  TGraph* gCont[6];

  TCanvas** cBs = new TCanvas*[6];
  TCanvas** cCh = new TCanvas*[6];

  TH1F* hIntegral_Bs[9][20];
  TH1F* hIntegral_Ch[9][20];

  char h0[10];
  char ha[10];
  for (int iw=0;iw<6;iw++){
    if(iw == 2)continue;
  //  for (int iw=0;iw<1;iw++){
    cBs[iw] = new TCanvas();
    cCh[iw] = new TCanvas();

    //  for (int iScan=1;iScan<nScan;iw++){
    for (int iScan=1;iScan<10;iScan++){
      sprintf (h0,"IntBs_Ch%d_Scan_%d",iw, iScan);
      hIntegral_Bs[iw][iScan-1] = (TH1F*)inF->Get(h0);
      hIntegral_Bs[iw][iScan-1]->SetDirectory(0);
      hIntegral_Bs[iw][iScan-1]->SetLineColor(colors.at(iw));
      hIntegral_Bs[iw][iScan-1]->GetXaxis()->SetTitle("Bs integral");
      hIntegral_Bs[iw][iScan-1]->SetLineWidth(2);
      hIntegral_Bs[iw][iScan-1]->Rebin(50);

      cBs[iw]->cd();
      if(iScan == 1) hIntegral_Bs[iw][iScan-1]->Draw();
      else hIntegral_Bs[iw][iScan-1]->Draw("same");

      sprintf (ha,"IntCh_Ch%d_Scan_%d",iw, iScan);
      hIntegral_Ch[iw][iScan-1] = (TH1F*)inF->Get(ha);
      hIntegral_Ch[iw][iScan-1]->SetDirectory(0);
      hIntegral_Ch[iw][iScan-1]->SetLineColor(colors.at(iw));
      hIntegral_Ch[iw][iScan-1]->GetXaxis()->SetTitle("Ch integral");
      hIntegral_Ch[iw][iScan-1]->SetLineWidth(2);
      hIntegral_Ch[iw][iScan-1]->Rebin(50);

      cCh[iw]->cd();
      if(iScan == 1) hIntegral_Ch[iw][iScan-1]->Draw();
      else hIntegral_Ch[iw][iScan-1]->Draw("same");
    }
//     gCont[iw] = new TGraph();
//     gCont[iw]->SetMarkerColor(colors.at(iw)); 
//     gCont[iw]->SetLineColor(colors.at(iw)); 
//     gCont[iw]->SetMarkerStyle(21); 
//     gCont[iw]->SetMarkerSize(0.6); 
  }

  return;
  TLegend* leg = new TLegend(0.88,0.65,0.98,0.85);
  leg->SetFillColor(0);
  leg->SetTextFont(41);
  for (int iw=0;iw<6;iw++) {
    if(iw == 2) continue;
    leg->AddEntry(hIntegral_Bs[iw][0],nameMCP.at(iw).c_str(),"l");
  }

  for(int iw=0;iw<6;iw++){
    //for(int iw=0;iw<1;iw++){
    cBs[iw]->cd();
    leg->Draw("same");
    cCh[iw]->cd();
    leg->Draw("same");
  }

  /*
  for(int iBin=0; iBin<80; ++iBin){
    for (int iw=0;iw<6;iw++){
      if(iw == 2) continue;
      double fake = -1;
      double num = hBaseLineAll[iw]->Integral(hBaseLineAll[iw]->FindBin(iBin), hBaseLineAll[iw]->GetNbinsX());
      double den = hBaseLineAll[iw]->Integral(1, hBaseLineAll[iw]->GetNbinsX());
      if(den != 0) fake = num/den;
      
      gCont[iw]->SetPoint(iBin, iBin, fake*100.);
    }
  }
  */
}


/////////////////////////////////////////////////////////////////////////////////

void doPLots_Integral(std::string fileName){
  TFile* inF = TFile::Open(fileName.c_str(), "read");
  
  int colorsExt[9] = {kRed+1, kBlue, kPink, kGreen+2, kAzure+7, kOrange+7, kCyan, kWhite, kBlack};
  std::vector<int> colors;
  for(int posVec = 0; posVec<6; ++posVec){
    colors.push_back(colorsExt[posVec]);
  }
  
  std::vector<std::string> nameMCP;
  nameMCP.push_back("MiB1");
  nameMCP.push_back("MiB2");
  nameMCP.push_back("ScB");
  nameMCP.push_back("Planacon");
  nameMCP.push_back("MiB3");
  nameMCP.push_back("RomaX");

  TGraph* gCont_Bs[6];
  TGraph* gCont_SoverB[6];

  TH1F* hIntegralAll_Bs[9];
  TH1F* hIntegralAll_Ch[9];

  char h0[10];
  for (int iw=0;iw<6;iw++){
    sprintf (h0,"IntABs_Ch%d",iw);
    hIntegralAll_Bs[iw] = (TH1F*)inF->Get(h0);
    hIntegralAll_Bs[iw]->SetLineColor(colors.at(iw));
    hIntegralAll_Bs[iw]->SetLineWidth(2);

    sprintf (h0,"IntACh_Ch%d",iw);
    hIntegralAll_Ch[iw] = (TH1F*)inF->Get(h0);
    hIntegralAll_Ch[iw]->SetLineColor(colors.at(iw));
    hIntegralAll_Ch[iw]->SetLineWidth(2);

    gCont_Bs[iw] = new TGraph();
    gCont_Bs[iw]->SetMarkerColor(colors.at(iw)); 
    gCont_Bs[iw]->SetLineColor(colors.at(iw)); 
    gCont_Bs[iw]->SetMarkerStyle(21); 
    gCont_Bs[iw]->SetMarkerSize(0.6); 

    gCont_SoverB[iw] = new TGraph();
    gCont_SoverB[iw]->SetMarkerColor(colors.at(iw)); 
    gCont_SoverB[iw]->SetLineColor(colors.at(iw)); 
    gCont_SoverB[iw]->SetMarkerStyle(21); 
    gCont_SoverB[iw]->SetMarkerSize(0.6); 
 }

  TLegend* leg = new TLegend(0.88,0.65,0.98,0.85);
  leg->SetFillColor(0);
  leg->SetTextFont(41);
  for (int iw=0;iw<6;iw++) {
    if(iw == 2) continue;
    leg->AddEntry(hIntegralAll_Bs[iw],nameMCP.at(iw).c_str(),"l");
  }

  TCanvas** c = new TCanvas*[6]; 
  c[0] = new TCanvas();
  hIntegralAll_Bs[0]->Draw();
  hIntegralAll_Ch[0]->Draw("same");
  for (int iw=1;iw<6;iw++){
    c[iw] = new TCanvas();
    c[iw]->cd();
    hIntegralAll_Bs[iw]->Draw();
    hIntegralAll_Ch[iw]->Draw("same");
  }
  leg->Draw("same");


  for(int iBin=0; iBin<2500; ++iBin){
    for (int iw=0;iw<6;iw++){
      if(iw == 2) continue;
      double fake = -1;
      double num = hIntegralAll_Bs[iw]->Integral(1, hIntegralAll_Bs[iw]->FindBin(iBin - 2000.));
      //      if(iw == 0) std::cout << " num = " << num <<  std::endl;
      double den = hIntegralAll_Bs[iw]->Integral(1, hIntegralAll_Bs[iw]->GetNbinsX());
      if(den != 0) fake = num/den;
      
//        if(iw == 0) {
// 	 std::cout << " iBin = " << iBin <<  std::endl;
// 	 std::cout << " num = " << num <<  std::endl;
// 	 std::cout << " den = " << den <<  std::endl;
// 	 std::cout << " fake = " << fake << " thresh = " << iBin <<  std::endl;
//        }
      gCont_Bs[iw]->SetPoint(iBin, iBin-2000., fake*100.);

      double SoverB = -1;
      int numS = hIntegralAll_Ch[iw]->Integral(1, hIntegralAll_Ch[iw]->FindBin(iBin - 2000.));
      if(numS != 0) SoverB = num/numS;
      gCont_SoverB[iw]->SetPoint(iBin, iBin-2000., SoverB*100.);
    }
  }

  TCanvas* c2 = new TCanvas();
  gPad->SetGrid();
//   gCont[0]->GetXaxis()->SetRangeUser(-500, 60);
//   gCont[0]->GetYaxis()->SetRangeUser(0, 105);
  gCont_Bs[0]->GetXaxis()->SetTitle("Integral threshold (ADC)");
  gCont_Bs[0]->GetYaxis()->SetTitle("noise contamination (%)");
  gCont_Bs[0]->Draw("apl");
  for (int iw=1;iw<6;iw++){
    if(iw == 2) continue;
    gCont_Bs[iw]->Draw("pl,same");
  }
  leg->Draw("same");

  TCanvas* c3 = new TCanvas();
  gPad->SetGrid();
//   gCont[0]->GetXaxis()->SetRangeUser(-500, 60);
//   gCont[0]->GetYaxis()->SetRangeUser(0, 105);
  gCont_SoverB[0]->GetXaxis()->SetTitle("Integral threshold (ADC)");
  gCont_SoverB[0]->GetYaxis()->SetTitle("S/B (%)");
  gCont_SoverB[0]->Draw("apl");
  for (int iw=1;iw<6;iw++){
    if(iw == 2) continue;
    gCont_SoverB[iw]->Draw("pl,same");
  }
  leg->Draw("same");

}


void doAltro(){



  int nScan = 16;

  int goodCh[3] = {1, 3, 4}; 

  std::vector<std::string> nameMCP;
  nameMCP.push_back("MiB1");
  nameMCP.push_back("MiB2");
  nameMCP.push_back("ScB");
  nameMCP.push_back("Planacon");
  nameMCP.push_back("MiB3");
  nameMCP.push_back("RomaX");


  TFile* inF = TFile::Open("outHistos.root");
  TGraph* hRatio[9];

  int colorsExt[9] = {kRed+1, kBlue, kPink, kGreen+1, kCyan+2, kOrange+7, kAzure+7, kWhite, kBlack};
  std::vector<int> colors;
  for(int posVec = 0; posVec<6; ++posVec){
    colors.push_back(colorsExt[posVec]);
  }

  TLegend* leg = new TLegend(0.88,0.65,0.98,0.85);
  leg->SetFillColor(0);
  leg->SetTextFont(41);
  for (int iw=0;iw<6;iw++) {
    leg->AddEntry(hRatio[iw],nameMCP.at(iw).c_str(),"l");
  }

  TCanvas* c = new TCanvas();

  TH1F* hAmpMax[9];
  TH1F* hBaseLine[9];

  float thres = 20;

  for(int iSC=1; iSC<=nScan+1; ++iSC){
  //  for(int iSC=1; iSC<=2; ++iSC){
    for(int iCh=0; iCh<6; ++iCh){
      if(iCh == 2) continue;

      char h1[10];
      sprintf (h1,"AM_Ch%d_Scan%d",iCh, iSC);
      hAmpMax[iCh] = (TH1F*)inF->Get("h1");


      char h2[10];
      sprintf (h2,"Bs_Ch%d_Scan%d",iCh, iSC);
      hBaseLine[iCh] = (TH1F*)inF->Get("h2");
    
      std::cout << "preso " << hAmpMax[iCh]->FindBin(thres) << std::endl;
    float ratio = (hAmpMax[iCh]->Integral( hAmpMax[iCh]->FindBin(thres), hAmpMax[iCh]->MaximumBin()) /
		   hBaseLine[iCh]->Integral( hBaseLine[iCh]->FindBin(thres), hBaseLine[iCh]->MaximumBin()) );
    hRatio->SetPoint(iSc, iSc, ratio);
    hRatio[iw]->SetMarkerColor(colors.at(iw));
    
    c->cd();
    if(iScan == 1 && iCh == 0)    hRatio[iw]->Draw("apl");
    else hRatio[iw]->Draw("p,same");
    }
  }
 
  c->cd();
  leg->Draw("same");
}
