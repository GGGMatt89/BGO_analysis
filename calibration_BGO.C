//This macro makes use of the calibration parameters extracted with the first raw data analyse to equalize the PM gains and produce a new file with the same
//format as the raw data but with the PM gain equalized. This means that each data value is multiplied by the calibration factor.
//The code is very similar to the one in analyse_BGO.C
//It produces the same plots after the calibration: ADC value spectra, flood map, monodimensional distributions, energy spectrum

//IT MUST BE APPLIED ON BOTH THE NA22 DATA and CO60 DATA for CALIBRATION

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


void calibration_BGO(Int_t blockNum, Int_t runNumber, Bool_t px_an = kFALSE){
  Double_t sat_limit = 3484.; //upper ADC limit - > reject saturated events
  Double_t ref_en = 1275.; //keV - reference energy for Na22 source - > second peak
  std::ostringstream ff;
  ff << blockNum;
  string f = ff.str();
  //Modify this two paths according to the path of the configurations files
  string confGAIN_filename = "../data/conf_files/confGAIN_7627.txt";//simple file with the 4 calibration factors -> it is created from scratch, it contains one line with 4 values separated by tabs
  string confPIX_filename = "../data/block_"+f+"/conf_files/confPIX.txt";

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
  //<-------------------------------------------------------------------------->//
  //Defining data file
  std::ostringstream ss;
  ss << std::setfill('0') << std::setw(5) << runNumber;
  string s = ss.str();
  cout <<"Calibration data RUN"<< s << endl;
  //string filename = "../data/Scan_Point_1mm_40steps/RUN"+ s +"/DataMeas_RUN"+s+".txt";
  //string calib_filename = "../data/Scan_Point_1mm_40steps/RUN"+ s +"/CalibratedData_RUN"+s+".txt";
  //Modify the file path according to the location of the files
  string filename = "../data/Scan_Point_1mm_40steps/RUN"+ s +"/DataMeas_RUN"+s+".txt";//Input file -> raw data
  string calib_filename = "../data/Scan_Point_1mm_40steps/RUN"+ s +"/CalibratedData_RUN"+s+".txt";//Output file -> calibrated data

  ofstream output_cal_txt(calib_filename);
  cout<<filename<<endl;
  //Variables for analysis and calibration
  string line;
  Double_t temp[8];
  vector <Double_t> PM0_data, PM1_data, PM2_data, PM3_data;//vectors for PM raw data
  vector <Double_t> PM0_data_cal, PM1_data_cal, PM2_data_cal, PM3_data_cal;//vectors for PM data after calibration and selection

  //<--------------------------------PIXEL SETTINGS---------------------------->//
  //useful if the pixel analysis is switched on -> variable px_an
  Double_t temp_pos[2];
  temp_pos[0] = 0.;
  temp_pos[1] = 0.;
  Int_t PixLine = 0;
  Int_t PixColumn = 0;
  Double_t pixel_pos[8][8][2];
  TH1F *single_pix_histo[8][8];
  ifstream confPIX;
  //it opens the file with the pixel positions and stores them in the array pixel_pos
  if(px_an){
    confPIX.open(confPIX_filename);
    getline(confPIX, line);
    while(confPIX>>PixLine>>PixColumn>>temp_pos[0]>>temp_pos[1]){
      pixel_pos[PixLine][PixColumn][0] = temp_pos[0];
      pixel_pos[PixLine][PixColumn][1] = temp_pos[1];
    }
    confPIX.close();

    for(Int_t xx = 0; xx<8; xx++){
      for(Int_t yy = 0; yy<8; yy++){
        single_pix_histo[xx][yy] = new TH1F(Form("pix %d-%d", xx, yy), Form("pix %d-%d", xx, yy), 800, -0.5, 1999.5);
        single_pix_histo[xx][yy]->SetLineColor(((xx+yy)%8)+1);
      }
    }
  }
  //<-------------------------------------------------------------------------->//

  //<-----------------------------PLOT DEFINITION------------------------------>//
  TH1F *PM0_histo = new TH1F("PM 0", "PM 0", 800, -0.5, 1999.5);
  TH1F *PM1_histo = new TH1F("PM 1", "PM 1", 800, -0.5, 1999.5);
  TH1F *PM2_histo = new TH1F("PM 2", "PM 2", 800, -0.5, 1999.5);
  TH1F *PM3_histo = new TH1F("PM 3", "PM 3", 800, -0.5, 1999.5);//single PM amplitude profiles - raw data

  TH1F *PM0_histo_cal = new TH1F("PM 0 cal", "PM 0 cal", 800, -0.5, 1999.5);
  TH1F *PM1_histo_cal = new TH1F("PM 1 cal", "PM 1 cal", 800, -0.5, 1999.5);
  TH1F *PM2_histo_cal = new TH1F("PM 2 cal", "PM 2 cal", 800, -0.5, 1999.5);
  TH1F *PM3_histo_cal = new TH1F("PM 3 cal", "PM 3 cal", 800, -0.5, 1999.5);//single PM amplitude profiles - calibrated and selected data
  TH1F *PM_sum_histo = new TH1F("PM sum cal", "PM sum cal", 800, -0.5, 1999.5);
  //<-------------------------------------------------------------------------->//
  //Opening data FILE
  ifstream input_data(filename);
  if(input_data.is_open()){
    getline(input_data, line);//skip the first line
    output_cal_txt<<line;
    //read and store the raw data in the vectors and plots
    while(input_data>>temp[0]>>temp[1]>>temp[2]>>temp[3]>>temp[4]>>temp[5]>>temp[6]>>temp[7]){
      PM0_data.push_back(temp[0]);
      PM0_histo->Fill(temp[0]);
      PM1_data.push_back(temp[1]);
      PM1_histo->Fill(temp[1]);
      PM2_data.push_back(temp[2]);
      PM2_histo->Fill(temp[2]);
      PM3_data.push_back(temp[3]);
      PM3_histo->Fill(temp[3]);
    }
  }else{
    cout<<"FILE NOT FOUND"<<endl;
  }

  //<------------------------------CALIBRATION--------------------------------->//
  cout<<"Start calibration ..."<<endl;
  //It starts with the correction of the offset. Maximum at low ADC values translated to 0
  //Offset calculation
  Double_t max_PM0, max_PM1, max_PM2, max_PM3;
  max_PM0 = PM0_histo->GetBinCenter(PM0_histo->GetMaximumBin())-12;
  max_PM1 = PM1_histo->GetBinCenter(PM1_histo->GetMaximumBin())-12;
  max_PM2 = PM2_histo->GetBinCenter(PM2_histo->GetMaximumBin())-12;
  max_PM3 = PM3_histo->GetBinCenter(PM3_histo->GetMaximumBin())-12;//everything is translated towards 0 according to the spectrum maximum value
  cout<<"Maximum bin PM 0 = "<<max_PM0<<endl;
  cout<<"Maximum bin PM 1 = "<<max_PM1<<endl;
  cout<<"Maximum bin PM 2 = "<<max_PM2<<endl;
  cout<<"Maximum bin PM 3 = "<<max_PM3<<endl;
  //Gain equalization
  Double_t PM_0_corr, PM_1_corr, PM_2_corr, PM_3_corr;//correction value extracted (manually) from raw amplitude spectra
  ifstream confGAIN(confGAIN_filename);//read the configuration file with the calibration factors
  getline(confGAIN, line);
  confGAIN>>PM_0_corr>>PM_1_corr>>PM_2_corr>>PM_3_corr;
  confGAIN.close();
  Double_t sum_4PM = 0.;
  Double_t temp_ADC_0 = 0.;
  Double_t temp_ADC_1 = 0.;
  Double_t temp_ADC_2 = 0.;
  Double_t temp_ADC_3 = 0.;
  Int_t rej_cal = 0;
  for(UInt_t i = 0; i<PM0_data.size(); i++){
    //Offset subtraction + gain correction
    if(PM0_data.at(i)>sat_limit||PM1_data.at(i)>sat_limit||PM2_data.at(i)>sat_limit||PM3_data.at(i)>sat_limit){
      //rejection of each events with at least on PM saturated
      rej_cal++;
      continue;
    }
    //calibration + offset
    temp_ADC_0 = PM0_data.at(i)*(ref_en/PM_0_corr) - max_PM0;
    temp_ADC_1 = PM1_data.at(i)*(ref_en/PM_1_corr) - max_PM1;
    temp_ADC_2 = PM2_data.at(i)*(ref_en/PM_2_corr) - max_PM2;
    temp_ADC_3 = PM3_data.at(i)*(ref_en/PM_3_corr) - max_PM3;

    //store the calibrated data in vectors and histos
    if(temp_ADC_0>0. && temp_ADC_1>0. && temp_ADC_2>0. &&temp_ADC_3>0.){
      PM0_data_cal.push_back(temp_ADC_0);
      PM0_histo_cal->Fill(temp_ADC_0);
      PM1_data_cal.push_back(temp_ADC_1);
      PM1_histo_cal->Fill(temp_ADC_1);
      PM2_data_cal.push_back(temp_ADC_2);
      PM2_histo_cal->Fill(temp_ADC_2);
      PM3_data_cal.push_back(temp_ADC_3);
      PM3_histo_cal->Fill(temp_ADC_3);
      sum_4PM = temp_ADC_0+temp_ADC_1+temp_ADC_2+temp_ADC_3;
      if(sum_4PM>100.){
        PM_sum_histo->Fill(sum_4PM);
      }
    }
    else{
      rej_cal ++;
      continue;
    }
  }
cout<<"N of rejected events for pedestal sub or saturation "<<rej_cal<<endl;
//<---------------------------------------------------------------------------->//

//Histograms styling
PM0_histo_cal->SetLineColor(kBlack);
PM0_histo_cal->GetXaxis()->SetTitle("Energy [keV]");
PM0_histo_cal->GetYaxis()->SetTitle("Entries");
PM0_histo_cal->SetStats(0);
PM0_histo_cal->SetTitle(NULL);
PM1_histo_cal->SetLineColor(kRed);
PM1_histo_cal->SetStats(0);
PM1_histo_cal->SetTitle(NULL);
PM2_histo_cal->SetLineColor(kGreen);
PM2_histo_cal->SetStats(0);
PM2_histo_cal->SetTitle(NULL);
PM3_histo_cal->SetLineColor(kBlue);
PM3_histo_cal->SetStats(0);
PM3_histo_cal->SetTitle(NULL);
PM_sum_histo->SetTitle(NULL);
PM_sum_histo->SetLineColor(kOrange);
PM_sum_histo->SetLineStyle(5);
PM_sum_histo->SetLineWidth(4);
TLegend *leg_PM = new TLegend(0.5,0.7,0.65,0.9);
leg_PM->SetHeader("PMT"); // option "C" allows to center the header
leg_PM->AddEntry(PM0_histo_cal,"0","l");
leg_PM->AddEntry(PM1_histo_cal,"1","l");
leg_PM->AddEntry(PM2_histo_cal,"2","l");
leg_PM->AddEntry(PM3_histo_cal,"3","l");
leg_PM->AddEntry(PM_sum_histo, "Sum", "l");


//PM calibrated amplitude
TCanvas *PM_canvas = new TCanvas("PM_canvas", "PM_canvas", 600, 400);
PM_canvas->SetFillColor(0); //
PM_canvas->SetBorderMode(0);	//
PM_canvas->SetLeftMargin(0.1409396); //
PM_canvas->SetRightMargin(0.14865772); //
gStyle->SetOptStat(000); //
PM_canvas->cd(2);
PM0_histo_cal->Draw();
PM1_histo_cal->Draw("same");
PM2_histo_cal->Draw("same");
PM3_histo_cal->Draw("same");
PM_sum_histo->Draw("same");
leg_PM->Draw();
PM_canvas->SetLogy();

//FLOOD MAP & SPATIAL DISTRIBUTIONS
Double_t PM_sum = 0.;
Double_t Xpos, Ypos;
TH2D *floodMap = new TH2D("floodMap", "floodMap", 200,-1, 1,200,-1,1);
TH1F *profileX_int = new TH1F("profile X int cal", "profile X int cal", 200, -1, 1);
TH1F *profileY_int = new TH1F("profile Y int cal", "profile Y int cal", 200, -1, 1);
TH1F *energy_spectrum = new TH1F("energy_spectrum_cal", "energy_spectrum_cal", 600, -0.5, 5999.5);
Int_t ADC_sel_up = 3000;
Int_t ADC_sel_down = 0;
Int_t en_sel_down = 200;
Int_t en_sel_up = 3500;
UInt_t N_PM_sel = 0;
UInt_t PM_on=0;
UInt_t NPM_evts_count[5];
for(Int_t init = 0; init < 5; init++){NPM_evts_count[init] = 0;}
UInt_t sat_evts = 0;
Double_t sum_THR = 0.;
Double_t distance = 0.;
Double_t tmp_max_dist = 2.;
Int_t chosen_pix[2];
chosen_pix[0]=0;
chosen_pix[1]=0;


for(UInt_t i = 0; i < PM0_data_cal.size(); i++ ){
  //for(UInt_t i = 0; i < 10; i++ ){
    PM_on = 0;
    distance = 0.;
    tmp_max_dist = 2.;
    //Check N of PM over THR
    if(PM0_data_cal.at(i)>ADC_sel_down)PM_on++;
    if(PM1_data_cal.at(i)>ADC_sel_down)PM_on++;
    if(PM2_data_cal.at(i)>ADC_sel_down)PM_on++;
    if(PM3_data_cal.at(i)>ADC_sel_down)PM_on++;
    NPM_evts_count[PM_on]++;
    //Check N of PM under saturation THR
    if((PM3_data_cal.at(i)>ADC_sel_up)|(PM1_data_cal.at(i)>ADC_sel_up)|(PM2_data_cal.at(i)>ADC_sel_up)|(PM0_data_cal.at(i)>ADC_sel_up)){
      sat_evts++;
      continue;
    }else{
    if(PM_on>=N_PM_sel){//Selection on number of PM over thr -> selection of events with at least N PMs with a signal over threshold
      PM_sum = PM0_data_cal.at(i)+PM1_data_cal.at(i)+PM2_data_cal.at(i)+PM3_data_cal.at(i);
      Xpos = ((PM3_data_cal.at(i)+PM1_data_cal.at(i))-(PM2_data_cal.at(i)+PM0_data_cal.at(i)))/PM_sum;
      Ypos = ((PM2_data_cal.at(i)+PM3_data_cal.at(i))-(PM1_data_cal.at(i)+PM0_data_cal.at(i)))/PM_sum;
      if(px_an){//done only if px_an is TRUE - not the default
        //it finds the pixel at the minimum distance from the reconstructed event and assigns the event to such a pixel - pixel position map in the file confPIX
        for(Int_t px = 0; px<8; px++){
          for(Int_t py = 0; py<8; py++){
            //cout<<"Pixel "<<px<<"-"<<py<<endl;
            distance = TMath::Sqrt(TMath::Power(pixel_pos[px][py][0]-Xpos,2)+TMath::Power(pixel_pos[px][py][1]-Ypos,2));
            //cout<<"calculated distance "<<distance<<endl;
            if(distance < tmp_max_dist){
              chosen_pix[0] = px;
              chosen_pix[1] = py;
              tmp_max_dist = distance;
            }
          }
        }
        //cout<<"Chosen pixel "<<chosen_pix[0]<<"-"<<chosen_pix[1]<<endl;
        //cout<<"pixel filled with "<<PM_sum<<endl;
        //according to the assignement of the events to the pixels, an energy spectrum per pixel is created
          if(PM_sum>sum_THR)single_pix_histo[chosen_pix[0]][chosen_pix[1]]->Fill(PM_sum);
      }
          if(PM_sum< en_sel_up && PM_sum>en_sel_down ){floodMap->Fill(Xpos, Ypos); //energy selection
          profileX_int->Fill(Xpos);
          profileY_int->Fill(Ypos);
          energy_spectrum->Fill(PM_sum);
          //fill the text file with the calibrated data -> 4 values of the calibrated PMs
          output_cal_txt<<PM0_data_cal.at(i)<<"\t"<<PM1_data_cal.at(i)<<"\t"<<PM2_data_cal.at(i)<<"\t"<<PM3_data_cal.at(i)<<"\t"<<endl;
        }
      }
    }
  }
  output_cal_txt.close();
  cout<<"Total number of detected events : "<<PM0_data_cal.size()<<endl;
  cout<<"Total number of selected events : "<<floodMap->GetEntries()<<endl;
  cout<<"Event classification : "<<endl;
  cout<<"0 PM on : "<<NPM_evts_count[0]<<endl;
  cout<<"1 PM on : "<<NPM_evts_count[1]<<endl;
  cout<<"2 PM on : "<<NPM_evts_count[2]<<endl;
  cout<<"3 PM on : "<<NPM_evts_count[3]<<endl;
  cout<<"4 PM on : "<<NPM_evts_count[4]<<endl;


  if(px_an){
    for(Int_t cox = 0; cox <8; cox++){
      for(Int_t coy = 0; coy<8; coy++){
        cout<<"Entries in pixel "<<"\t"<<cox<<"-"<<coy<<" = "<<single_pix_histo[cox][coy]->GetEntries()<<endl;
        single_pix_histo[cox][coy]->Scale(1/single_pix_histo[cox][coy]->GetMaximum());
      }
    }
  }
  floodMap->GetXaxis()->SetTitle("Relative X position");
  floodMap->GetYaxis()->SetTitle("Relative Y position");
  floodMap->GetZaxis()->SetTitle("Entries");

  Int_t count = 0;
  TGraph *pixels = new TGraph(64);

  if(px_an){
    for(Int_t col = 0; col<8; col++){
      for(Int_t lin = 0; lin < 8; lin++){
        pixels->SetPoint(count, pixel_pos[col][lin][0], pixel_pos[col][lin][1]);
        count++;
      }
    }
    pixels->SetMarkerStyle(3);
    pixels->SetMarkerColor(kGreen);
    pixels->SetMarkerSize(1);
  }

  TCanvas *flood_canvas = new TCanvas("flood_canvas","flood_canvas ", 600, 500);
  flood_canvas->SetFillColor(0); //
  flood_canvas->SetBorderMode(0);
  flood_canvas->SetLeftMargin(0.1409396); //
  flood_canvas->SetRightMargin(0.14865772); //
  gStyle->SetOptStat(000); //
  floodMap->SetTitle(NULL);
  gStyle->SetPalette(55);
  floodMap->Draw("colz");
  if(px_an)pixels->Draw("P same");
  flood_canvas->Update();

  profileX_int->GetXaxis()->SetTitle("Relative position");
  profileX_int->GetYaxis()->SetTitle("Entries");
  profileX_int->SetLineColor(kBlue);

  profileY_int->GetXaxis()->SetTitle("Relative position");
  profileY_int->GetYaxis()->SetTitle("Entries");
  profileY_int->SetLineColor(kRed);


  TCanvas *profile_canvas = new TCanvas("profile_canvas_cal","profile_canvas_cal", 1000, 400);
  profile_canvas->SetFillColor(0); //
  profile_canvas->SetBorderMode(0);
  profile_canvas->SetLeftMargin(0.1409396); //
  profile_canvas->SetRightMargin(0.14865772); //
  profile_canvas->Divide(2,1);
  profile_canvas->cd(1);
  profileX_int->Draw();
  profile_canvas->cd(2);
  profileY_int->Draw();


  energy_spectrum->GetXaxis()->SetTitle("ADC channel");
  energy_spectrum->GetYaxis()->SetTitle("Entries");

  TCanvas *Espectrum_canvas = new TCanvas("Espectrum_canvas_cal","Espectrum_canvas_cal", 600, 400);
  Espectrum_canvas->SetFillColor(0); //
  Espectrum_canvas->SetBorderMode(0);
  Espectrum_canvas->SetLeftMargin(0.1409396); //
  Espectrum_canvas->SetRightMargin(0.14865772); //
  energy_spectrum->SetTitle(NULL);
  energy_spectrum->Draw();

  TCanvas *singlePix_canvas = new TCanvas("singlePix_canvas","singlePix_canvas", 600, 400);;
  if(px_an){
    singlePix_canvas->SetFillColor(0); //
    singlePix_canvas->SetBorderMode(0);
    singlePix_canvas->SetLeftMargin(0.1409396); //
    singlePix_canvas->SetRightMargin(0.14865772); //
    singlePix_canvas->SetTitle(NULL);
    for(Int_t hx = 1; hx<7; hx++){
      for(Int_t hy = 1; hy<7; hy++){
        if(hx==0&&hy==0){single_pix_histo[hx][hy]->Draw();}
        else{single_pix_histo[hx][hy]->Draw("same");}
      }
    }
  }

  TCanvas *singlePix_HC[8][8];
  //OUTPUT ROOT FILE
  //TString rootfilename = "../data/Scan_Point_1mm_40steps/RUN"+s+"/outputCalibration_RUN" + s +".root";
  TString rootfilename = "../data/Scan_Point_1mm_40steps/RUN"+s+"/outputCalibration_RUN" + s +".root";
  TFile root_out(rootfilename, "RECREATE");
  root_out.cd();


  if(px_an){
    for(Int_t hhx = 0; hhx<8; hhx++){
      for(Int_t hhy = 0; hhy<8; hhy++){
        singlePix_HC[hhx][hhy] = new TCanvas(Form("Pix %d-%d HC",hhx, hhy), Form("Pix %d-%d HC",hhx, hhy), 600, 400);
        single_pix_histo[hhx][hhy]->Draw();
        singlePix_HC[hhx][hhy]->SetLogy();
        singlePix_HC[hhx][hhy]->Write(Form("single pix %d-%d", hhx, hhy));
      }
    }
  }
  PM_canvas->Write("Cal_PMAmp");
  flood_canvas->Write("Cal_FLOODMAP");
  profile_canvas->Write("Cal_profileXY_int");
  Espectrum_canvas->Write("Cal_energy_spectrum");
  profileX_int->Write("X_profile");
  profileY_int->Write("Y_profile");
  if(px_an)singlePix_canvas->Write("Single_pix_ampliDistr");


  //deleting pointers
  delete PM0_histo;
  delete PM1_histo;
  delete PM2_histo;
  delete PM3_histo;
  delete PM0_histo_cal;
  delete PM1_histo_cal;
  delete PM2_histo_cal;
  delete PM3_histo_cal;
  delete energy_spectrum;
  delete profileY_int;
  delete profileX_int;
  //delete floodMap;
  delete PM_canvas;
  //delete flood_canvas;
  delete profile_canvas;
  delete Espectrum_canvas;
  delete leg_PM;
  delete singlePix_canvas;

  return;

  }
