//this macro is the last analysis step, where the identified pixel positions and the pixel borders are used to re-analyse the single events and assign
//each events to a single pixel. After that, the energy spectrum of each pixel can be produced and analysed to equalized the block energy response on a
//single pixel basis.

#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <limits>
#include <math.h>
#include <vector>
#include <Riostream.h>

// Root
#include "TCanvas.h"
#include "TLegend.h"
#include "TFile.h"
#include "TTree.h"
#include "TLeaf.h"
#include "TH2D.h"
#include "TH1F.h"
#include "TSystem.h"
#include "TF2.h"
#include "TF1.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TString.h"
#include "TGraph.h"
#include "TMath.h"
#include "TSpectrum.h"
#include "functions.hh"


void pixel_Ecal(Int_t blockNum, Int_t runNumber){
  //string confPIX_filename = "../data/conf_files/confPIX_3166.txt";
  //<-------------------------GRAPHICAL STYLE-------------------------------->//
  gROOT->SetStyle("Plain");
  gROOT->ForceStyle();
  gStyle -> SetStatW(0.28);
  gStyle -> SetStatH(0.13);
  gStyle -> SetStatColor(0);
  gStyle -> SetStatX(0.87);
  gStyle -> SetStatY(0.85);
  gStyle -> SetStatFont(0);
  gStyle -> SetOptStat(111);
  gStyle->SetPalette(55);
  gErrorIgnoreLevel=kError;
  //<-------------------------------------------------------------------------->//
  std::ostringstream ss;
  ss << std::setfill('0') << std::setw(5) << runNumber;
  string s = ss.str();
  std::ostringstream ff;
  ff << blockNum;
  string f = ff.str();
  string inpixels_name = "../data/block_"+f+"/conf_files/confPIX.txt";
  string inData_name = "../data/block_"+f+"/RUN"+s+"/CalibratedData_RUN" + s +".txt";
  string incolxylinxy_name = "../data/block_"+f+"/conf_files/colxy_linxy.txt";
  //Histos for 1D integrated profiles
  //<--------------------------------PIXEL SETTINGS---------------------------->//
  Double_t ref_en = 1275.;
  TH1F *single_pix_histo[8][8];
  vector <vector <Double_t>> pixel_data;
  for(Int_t fl = 0; fl <64; fl ++){
    pixel_data.push_back({0});
  }
  Int_t pix_entries[8][8];
  Double_t PM_sum = 0.;
  Double_t sum_THR = 0.;
  Double_t distance = 0.;
  Double_t tmp_max_dist = 2.;
  Int_t chosen_pix[2];
  chosen_pix[0]=0;
  chosen_pix[1]=0;
  Int_t pourcentage=-1;
  Int_t pnum;
  Double_t pcxcoor;
  Double_t pcycoor;
  Double_t plxcoor;
  Double_t plycoor;
  Double_t Cx[56];//average values between pixel positions along X and Y
  Double_t Cy[56];
  Double_t Lx[56];
  Double_t Ly[56];
  Double_t Xpos = 0., Ypos = 0.;
  string line;

  ifstream confcxylxy;
  confcxylxy.open(incolxylinxy_name);//reads the position of the average points between pixels from the text file produced by pixel_id.C
  //getline(confcxylxy, line);
  while(confcxylxy>>pnum>>pcxcoor>>pcycoor>>plxcoor>>plycoor){
    Cx[pnum]=pcxcoor;
    Cy[pnum]=pcycoor;
    Lx[pnum]=plxcoor;
    Ly[pnum]=plycoor;
  }
  confcxylxy.close();

  Double_t temp[4];
  vector <Double_t> PM0_data, PM1_data, PM2_data, PM3_data;
  //reads the calibrated data event per event for the 4 PMs
  ifstream indata;
  indata.open(inData_name);
  if(indata.is_open()){
    getline(indata, line);//skip the first line
    while(indata>>temp[0]>>temp[1]>>temp[2]>>temp[3]){
      PM0_data.push_back(temp[0]);
      PM1_data.push_back(temp[1]);
      PM2_data.push_back(temp[2]);
      PM3_data.push_back(temp[3]);
    }
  }else{
    cout<<"FILE NOT FOUND"<<endl;
    return;
  }
  Int_t temp_index = 0;
  for(Int_t xx = 0; xx<8; xx++){
    for(Int_t yy = 0; yy<8; yy++){
      single_pix_histo[xx][yy] = new TH1F(" ", " ", 200, -0.5, 1999.5);
      single_pix_histo[xx][yy]->SetLineColor(((xx+yy)%8)+1);
    }
  }

Int_t prctg = floor(PM0_data.size()/100);
  for(UInt_t i = 0; i < PM0_data.size(); i++ ){
    //scan of all the events

    if(i%prctg==0){//used to create an analysis advancement bar
      pourcentage++;
      if (pourcentage<10){cout << pourcentage << "%" << endl;}
      if (pourcentage>=10 && pourcentage%10==0){cout << pourcentage << "%" << endl;}
    }
//starts the analysis by reconstructing the event interaction position and defining the points of minimal distance (along X and Y) and the relative position (right/left/up/down)
    distance = 0.;
    tmp_max_dist = 2.;
    PM_sum = 0.;
    PM_sum = PM0_data.at(i)+PM1_data.at(i)+PM2_data.at(i)+PM3_data.at(i);
    Xpos = ((PM3_data.at(i)+PM1_data.at(i))-(PM2_data.at(i)+PM0_data.at(i)))/PM_sum;
    Ypos = ((PM2_data.at(i)+PM3_data.at(i))-(PM1_data.at(i)+PM0_data.at(i)))/PM_sum;
    for(Int_t o=0; o<56; o++){
       distance = TMath::Sqrt(TMath::Power(Cx[o]-Xpos,2)+TMath::Power(Cy[o]-Ypos,2));
       if(distance < tmp_max_dist){
          tmp_max_dist = distance;
          if (Cx[o]-Xpos<0){//on est à droite
             chosen_pix[0] = (o%8);
             chosen_pix[1] = floor(o/8)+1;
          }
          if (Cx[o]-Xpos>0){//on est à gauche
             chosen_pix[0] = (o%8);
             chosen_pix[1] = floor(o/8);
          }
       }
       distance = TMath::Sqrt(TMath::Power(Lx[o]-Xpos,2)+TMath::Power(Ly[o]-Ypos,2));
       if(distance < tmp_max_dist){
          tmp_max_dist = distance;
          if (Ly[o]-Ypos<0){//on est au dessus
             chosen_pix[0] = floor(o/8)+1;
             chosen_pix[1] = (o%8);
          }
          if (Ly[o]-Ypos>0){//on est en dessous
             chosen_pix[0] = floor(o/8);
             chosen_pix[1] = (o%8);
          }

       }
    }
    single_pix_histo[chosen_pix[0]][chosen_pix[1]]->Fill(PM_sum);
    temp_index = chosen_pix[0]*8+chosen_pix[1];
    pixel_data[temp_index].push_back(PM_sum);
    //the pixel where to assign the event is univocally defined
  }

  Int_t max_entries = 0;
  Int_t max_pix[2];
  max_pix[0] = 0;
  max_pix[1] = 0;
  //traces the histo of the entries of each pixel
  for(Int_t cox = 0; cox <8; cox++){
    for(Int_t coy = 0; coy<8; coy++){
      pix_entries[cox][coy] = single_pix_histo[cox][coy]->GetEntries();
      if(pix_entries[cox][coy]>max_entries){max_entries = pix_entries[cox][coy];
        max_pix[0] = cox;
        max_pix[1] = coy;
    }
      single_pix_histo[cox][coy]->Scale(1/single_pix_histo[cox][coy]->GetMaximum());
      single_pix_histo[cox][coy]->GetXaxis()->SetTitle("Energy [keV]");
      single_pix_histo[cox][coy]->GetYaxis()->SetTitle("Normalized counts");
      single_pix_histo[cox][coy]->SetStats(kFALSE);
    }
  }

  //once the events are assigned, the energy calibration can start
  //ENERGY CALIBRATION ON SIGNLE PIXEL BASIS
  TCanvas *single_line_pixels[8];
  TGraph *temporary_peak[8][8];
  Double_t peaks_pos[8][8][2];

  Double_t temp_x, temp_y;
  for(Int_t slpx = 0; slpx<8; slpx++){
  single_line_pixels[slpx] = new TCanvas(Form("line %d", slpx), Form("line %d", slpx), 1400, 1400);
  single_line_pixels[slpx]->SetFillColor(0); //
  single_line_pixels[slpx]->SetBorderMode(0);
  single_line_pixels[slpx]->SetLeftMargin(0.1409396); //
  single_line_pixels[slpx]->SetRightMargin(0.14865772); //
  single_line_pixels[slpx]->Divide(4,2);

    for(Int_t slpy = 0; slpy<8; slpy++){
      temporary_peak[slpx][slpy] = new TGraph(2);
      temporary_peak[slpx][slpy]->SetMarkerStyle(23);
      temporary_peak[slpx][slpy]->SetMarkerSize(1);
      temporary_peak[slpx][slpy]->SetMarkerColor(2);
      single_line_pixels[slpx]->cd(slpy+1);
      single_pix_histo[slpx][slpy]->Smooth();
      single_pix_histo[slpx][slpy]->GetXaxis()->SetRangeUser(350, 650);
      temp_x = single_pix_histo[slpx][slpy]->GetBinCenter(single_pix_histo[slpx][slpy]->GetMaximumBin());
      temp_y = single_pix_histo[slpx][slpy]->GetBinContent(single_pix_histo[slpx][slpy]->GetMaximumBin());
      peaks_pos[slpx][slpy][0] = temp_x;
      temporary_peak[slpx][slpy]->SetPoint(0, temp_x, temp_y);
      single_pix_histo[slpx][slpy]->GetXaxis()->SetRangeUser(800, 1800);
      temp_x = single_pix_histo[slpx][slpy]->GetBinCenter(single_pix_histo[slpx][slpy]->GetMaximumBin());
      temp_y = single_pix_histo[slpx][slpy]->GetBinContent(single_pix_histo[slpx][slpy]->GetMaximumBin());
      peaks_pos[slpx][slpy][1] = temp_x;
      single_pix_histo[slpx][slpy]->GetXaxis()->SetRange(0, 800);
      temporary_peak[slpx][slpy]->SetPoint(1, temp_x, temp_y);
      single_pix_histo[slpx][slpy]->Draw("histo");
      temporary_peak[slpx][slpy]->Draw("P same");
    }
  }

  TH1F *single_pix_histo_cal[8][8];
  for(Int_t pxs = 0; pxs < 8; pxs++){
    for(Int_t pys = 0; pys < 8; pys++){
      single_pix_histo_cal[pxs][pys] = new TH1F(" ", " ", 200, -0.5, 1999.5);
      single_pix_histo_cal[pxs][pys]->SetLineColor(((pxs+pys)%8)+1);
    for(UInt_t entr = 1; entr<pixel_data[pxs*8+pys].size(); entr++){
      single_pix_histo_cal[pxs][pys]->Fill(pixel_data[pxs*8+pys].at(entr)*(ref_en/peaks_pos[pxs][pys][1]));
    }
    single_pix_histo_cal[pxs][pys]->Scale(1/single_pix_histo_cal[pxs][pys]->GetMaximum());
    single_pix_histo_cal[pxs][pys]->GetXaxis()->SetTitle("Energy [keV]");
    single_pix_histo_cal[pxs][pys]->GetYaxis()->SetTitle("Normalized counts");
    single_pix_histo_cal[pxs][pys]->SetStats(kFALSE);
  }
}
TCanvas *singlePix_canvas_simple = new TCanvas(" "," ", 1200, 800);
  singlePix_canvas_simple->SetFillColor(0); //
  singlePix_canvas_simple->SetBorderMode(0);
  singlePix_canvas_simple->SetLeftMargin(0.1409396); //
  singlePix_canvas_simple->SetRightMargin(0.14865772); //
  singlePix_canvas_simple->SetTitle(NULL);

  //singlePix_canvas_simple->Divide(2,1);
  //singlePix_canvas_simple->cd(1);//single pixel not calibrated
  //gPad->SetLogy();
  //gPad->SetTitle(NULL);
  for(Int_t hxn = 1; hxn<7; hxn++){
    for(Int_t hyn = 1; hyn<7; hyn++){
      single_pix_histo[hxn][hyn]->GetYaxis()->SetRangeUser(0, 1.1);
      if(hxn==0&&hyn==0){
        single_pix_histo[hxn][hyn]->Draw("C");}
      else{single_pix_histo[hxn][hyn]->Draw("C same");}
    }
  }
singlePix_canvas_simple->Update();

    TCanvas *singlePix_canvas = new TCanvas("singlePix_canvas","singlePix_canvas", 1000, 800);
      singlePix_canvas->SetFillColor(0); //
      singlePix_canvas->SetBorderMode(0);
      singlePix_canvas->SetLeftMargin(0.1409396); //
      singlePix_canvas->SetRightMargin(0.14865772); //
      singlePix_canvas->SetTitle(NULL);
      singlePix_canvas->Divide(2,1);
      singlePix_canvas->cd(1);//single pixel not calibrated
//      gPad->SetLogy();
      gPad->SetTitle(NULL);
      for(Int_t hxn = 1; hxn<7; hxn++){
        for(Int_t hyn = 1; hyn<7; hyn++){
          if(hxn==0&&hyn==0){single_pix_histo[hxn][hyn]->Draw("C");}
          else{single_pix_histo[hxn][hyn]->Draw("C same");}
        }
      }
      singlePix_canvas->cd(2);//single pixel calibrated
//      gPad->SetLogy();
      gPad->SetTitle(NULL);
      for(Int_t hx = 1; hx<7; hx++){
        for(Int_t hy = 1; hy<7; hy++){
          if(hx==0&&hy==0){single_pix_histo_cal[hx][hy]->Draw("C");}
          else{single_pix_histo_cal[hx][hy]->Draw("C same");}
        }
      }


//CREATION OF ENERGY SPECTRUM PLOTS FOR NON CALIBRATED AND CALIBRATED PIXELS
TH1F *Espect_noCal = new TH1F(" ","Energy spectrum before calibration", 200, -0.5, 1999.5);
TH1F *Espect_Cal = new TH1F(" ","Energy spectrum after calibration", 200, -0.5, 1999.5);
Int_t N_bins = single_pix_histo[0][0]->GetNbinsX();
  for(UInt_t lx = 0; lx<8; lx++){
    for(UInt_t ly = 0; ly<8; ly++){
      if(lx == 0 && ly == 0)continue;//the corners are excluded because their behaviour is different from the pthers due to the light loss
      if(lx == 0 && ly == 7)continue;
      if(lx == 7 && ly == 0)continue;
      if(lx == 7 && ly == 7)continue;
      for(Int_t nb = 0; nb<N_bins; nb++){
        Espect_noCal->Fill(single_pix_histo[lx][ly]->GetBinCenter(nb),single_pix_histo[lx][ly]->GetBinContent(nb));
        Espect_Cal->Fill(single_pix_histo_cal[lx][ly]->GetBinCenter(nb),single_pix_histo_cal[lx][ly]->GetBinContent(nb));
      }
    }
  }
  //fitting functions -> to extract the energy resolution of the block
  TF1 *g1 = new TF1("g1","gaus",300,650);
  TF1 *g2 = new TF1("g2","gaus",1000,1500);
  TF1 *g3 = new TF1("g3","gaus",400,600);
  TF1 *g4 = new TF1("g4","gaus",1100,1450);
  g1->SetLineColor(kRed);
  g2->SetLineColor(kRed);
  g3->SetLineColor(kRed);
  g4->SetLineColor(kRed);
  Double_t fit_par[12];
  TCanvas *Espectrum = new TCanvas("E spectra", "E spectra", 1200, 800);
  Espectrum->SetFillColor(0); //
  Espectrum->SetBorderMode(0);
  Espectrum->SetLeftMargin(0.1409396); //
  Espectrum->SetRightMargin(0.14865772); //
  Espectrum->SetTitle(NULL);
  Espectrum->Divide(2,1);
  Espectrum->cd(1);//Espectrum not calibrated
//  gPad->SetLogy();
  gPad->SetTitle(NULL);
  Espect_noCal->SetLineColor(kBlue);
  Espect_noCal->GetXaxis()->SetTitle("Energy [keV]");
  Espect_noCal->GetYaxis()->SetTitle("Normalized counts");
  Espect_noCal->Scale(1/Espect_noCal->GetMaximum());
  Espect_noCal->SetStats(kFALSE);
  Espect_noCal->Draw("hist C");
  Espect_noCal->Fit(g1, "EMR");
  Espect_noCal->Fit(g2, "EMR+");
  g1->Draw("same");
  g2->Draw("same");
  Espectrum->cd(2);
//  gPad->SetLogy();
  gPad->SetTitle(NULL);
  Espect_Cal->SetLineColor(kBlue);
  Espect_Cal->GetXaxis()->SetTitle("Energy [keV]");
  Espect_Cal->GetYaxis()->SetTitle("Normalized counts");
  Espect_Cal->Scale(1/Espect_Cal->GetMaximum());
  Espect_Cal->SetStats(kFALSE);
  Espect_Cal->Draw("hist C");
  Espect_Cal->Fit(g3, "EMR");
  Espect_Cal->Fit(g4, "EMR+");
  g3->Draw("same");
  g4->Draw("same");
  g1->GetParameters(&fit_par[0]);
  g2->GetParameters(&fit_par[3]);
  g3->GetParameters(&fit_par[6]);
  g4->GetParameters(&fit_par[9]);

    TCanvas *entries_c = new TCanvas("entries histo canvas","entries histo canvas", 600, 500 );
    TH1F *histo_entries[8];
    for(Int_t en = 0; en <8; en ++){
      histo_entries[en] = new TH1F(Form("row %d", en), Form("row %d", en), 8, 0, 8);
      histo_entries[en]->SetLineColor(en + 1);
      histo_entries[en]->GetXaxis()->SetTitle("Pixel #");
      histo_entries[en]->GetYaxis()->SetTitle("Counts [#]");
      for(Int_t px_en = 0; px_en <8; px_en++){
        histo_entries[en]->Fill(px_en, pix_entries[en][px_en]);
      }
    }
    entries_c->cd();
    histo_entries[max_pix[0]]->Draw("histo");
    for(Int_t xx = 0; xx<8; xx++){
        entries_c->cd();
        if( xx != max_pix[0]){histo_entries[xx]->Draw(" same histo");}
    }

    TCanvas *entries_single = new TCanvas("entries histo canvas single","entries histo canvas single", 1400, 1400);
    entries_single->SetFillColor(0); //
    entries_single->SetBorderMode(0);
    entries_single->SetLeftMargin(0.1409396); //
    entries_single->SetRightMargin(0.14865772); //
    entries_single->Divide(4,2);
    for(Int_t yy = 0; yy<8; yy++){
      entries_single->cd(yy+1);
      gPad->SetTitle(NULL);
      histo_entries[yy]->GetYaxis()->SetRangeUser(0, max_entries+1000);
      histo_entries[yy]->SetStats(kFALSE);
      histo_entries[yy]->Draw("histo");
    }

    cout<<"PARAMETERS FOR NON CALIBRATED DATA "<<endl;
    cout<<"@ 511 keV : Mean value = "<<fit_par[1]<<" --- FWHM = "<<fit_par[2]*2.355<<" -> E Resolution is "<<((fit_par[2]*2.355)/fit_par[1])*100<<" % FWHM"<<endl;
    cout<<"@ 1275 keV : Mean value = "<<fit_par[4]<<" --- FWHM = "<<fit_par[5]*2.355<<" -> E Resolution is "<<((fit_par[5]*2.355)/fit_par[4])*100<<" % FWHM"<<endl;
    cout<<"PARAMETERS FOR CALIBRATED DATA "<<endl;
    cout<<"@ 511 keV : Mean value = "<<fit_par[7]<<" --- FWHM = "<<fit_par[8]*2.355<<" -> E Resolution is "<<((fit_par[8]*2.355)/fit_par[7])*100<<" % FWHM"<<endl;
    cout<<"@ 1275 keV : Mean value = "<<fit_par[10]<<" --- FWHM = "<<fit_par[11]*2.355<<" -> E Resolution is "<<((fit_par[11]*2.355)/fit_par[10])*100<<" % FWHM"<<endl;

    TString output_root = "../data/block_"+f+"/RUN"+s+"/pixels_CalibratedOUT_RUN" + s +".root";
    TFile *outROOT = new TFile(output_root, "RECREATE");
    if(outROOT->IsOpen())cout<<"OUTPUT FILE OPENED"<<endl;
    //outROOT->cd();
    for(Int_t sa = 0; sa<8; sa++){
      single_line_pixels[sa]->Write(Form("Single line pixel %d", sa));
    }
    entries_single->Write("Entries single pixels per row");
    entries_c->Write("Entries single pixels per row overlap");
    Espectrum->Write("E spectra log");
    singlePix_canvas->Write("Single pixel spectra overlap");

    outROOT->Close();

    return;

}
//<-------------------------------------------------------------------------->//
