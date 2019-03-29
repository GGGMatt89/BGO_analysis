//BGO data format with PXIe acquisition
//Raw data file name -> DataMeas_RUNXXXXX.txt   (NB: run number always with 5 digits)
//Header -> 1 line with channel names (0 to 15)
//data with one event per line, signal max amplitude per channel with tab separation
//EXAMPLE: with 8 channels (standard format - in general only one block connected, but 1 acquisition card switched on so 8 channels minimum - only 4 used)
//CH0	CH1	CH2	CH3	CH4	CH5	CH6	CH7	CH8	CH9	CH10	CH11	CH12	CH13	CH14	CH15
//56.200	223.880	10.941	10.307	261.921	10.416	11.704	9.773
//47.850	12.217	444.840	37.111	520.938	12.175	12.281	10.940
//.......

//analyse_BGO.C
//This macro reads the raw data file and draw and save the following plots:
// - Histo of the 4/8 channels ADC energy_spectrum
// - Flood map (interaction position 2D map)
// - 1D spatial distributions of events on X and Y
// - Raw energy spectrum
//Thanks to these plots, the 4 calibration parameters can be extracted and used in the following analysis step -> calibration_BGO.C

//IT MUST BE APPLIED ON THE NA22 DATA

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


void analyse_BGO(Int_t blockNumber, Int_t runNumber){ //studied block and related run number
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
gStyle->SetPalette(kRainBow);

//settings to create the file name related to the block number and run number
std::ostringstream ss;
ss << std::setfill('0') << std::setw(5) << runNumber; //ss is filled with the run number formatted with 5 digits. Example: run 53 -> 00053
string s = ss.str();
std::ostringstream ff;
ff << blockNumber;
string f = ff.str();
cout <<"Analysis RUN"<< s << endl;
string filename = "../data/block_"+f+"/RUN"+ s +"/DataMeas_RUN"+s+".txt"; //first part to be changed according to the data path
//string filename = "../data/Scan_Line_22steps_2mm_200kEv/RUN"+ s +"/DataMeas_RUN"+s+".txt";
cout<<filename<<endl;


string line;
Double_t temp[8];//array used to temporary storage of the 8 values of each data line
vector <Double_t> PM0_data, PM1_data, PM2_data, PM3_data; //vector for the raw data of the 4 PMs

//PLOT DEFINITION - ADC spectra of the 4 channels (PMs) of a single block
TH1F *PM0_histo = new TH1F("PM 0", "PM 0", 800, -0.5, 3999.5);
TH1F *PM1_histo = new TH1F("PM 1", "PM 1", 800, -0.5, 3999.5);
TH1F *PM2_histo = new TH1F("PM 2", "PM 2", 800, -0.5, 3999.5);
TH1F *PM3_histo = new TH1F("PM 3", "PM 3", 800, -0.5, 3999.5);

//data FILE
ifstream input_data(filename);

if(input_data.is_open()){
  getline(input_data, line);//skip the first line of the file (header) -> I get tge line and it is not used

//for each data line, the values are temporary stored in temp[i] and then used to fill the histograms and the vectors - only first 4 columns of the file used, where the block was connected
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

//histograms styling
PM0_histo->SetLineColor(kBlack);
PM0_histo->GetXaxis()->SetTitle("ADC channel");
PM0_histo->GetYaxis()->SetTitle("Entries");
PM0_histo->SetStats(0);
PM0_histo->SetTitle(NULL);
PM1_histo->SetLineColor(kRed);
PM1_histo->SetStats(0);
PM1_histo->SetTitle(NULL);
PM2_histo->SetLineColor(kGreen);
PM2_histo->SetStats(0);
PM2_histo->SetTitle(NULL);
PM3_histo->SetLineColor(kBlue);
PM3_histo->SetStats(0);
PM3_histo->SetTitle(NULL);

//legend for the canvas with the superposition of the 4 channel spectra
TLegend *leg_PM = new TLegend(0.5,0.7,0.65,0.9);
leg_PM->SetHeader("PMT"); // option "C" allows to center the header
leg_PM->AddEntry(PM0_histo,"0","l");
leg_PM->AddEntry(PM1_histo,"1","l");
leg_PM->AddEntry(PM2_histo,"2","l");
leg_PM->AddEntry(PM3_histo,"3","l");


//PM raw amplitude
TCanvas *PM_canvas = new TCanvas("PM_canvas", "PM_canvas", 600, 400);
PM_canvas->SetFillColor(0); //
PM_canvas->SetBorderMode(0);	//
PM_canvas->SetLeftMargin(0.1409396); //
PM_canvas->SetRightMargin(0.14865772); //
gStyle->SetOptStat(000); //

PM0_histo->Draw();
PM1_histo->Draw("same");
PM2_histo->Draw("same");
PM3_histo->Draw("same");
leg_PM->Draw();
PM_canvas->SetLogy();//log scale on y axis

//RAW FLOOD MAP & SPATIAL DISTRIBUTIONS
Double_t Xpos, Ypos;
TH2D *floodMap = new TH2D(" ", " ", 200,-1, 1,200,-1,1);//2D position map
TH1F *profileX_int = new TH1F("profile X int", "profile X int", 200, -1, 1);//1D distribution on X
TH1F *profileY_int = new TH1F("profile Y int", "profile Y int", 200, -1, 1);//1D distribution on Y
TH1F *energy_spectrum = new TH1F("energy_spectrum", "energy_spectrum", 600, -0.5, 5999.5);//energy spectrum -> deposited energy proportional to the sum of the 4 signal amplitudes (because the signals result from a shaping amplifier step)
Double_t PM_sum = 0.;

for(UInt_t i = 0; i < PM0_data.size(); i++ ){//For each event, the interaction position is calculated with a center of gravity of the 4 signal amplitudes
  PM_sum = PM0_data.at(i)+PM1_data.at(i)+PM2_data.at(i)+PM3_data.at(i);
  Xpos = ((PM3_data.at(i)+PM1_data.at(i))-(PM2_data.at(i)+PM0_data.at(i)))/PM_sum;
  Ypos = ((PM2_data.at(i)+PM3_data.at(i))-(PM1_data.at(i)+PM0_data.at(i)))/PM_sum;

  floodMap->Fill(Xpos, Ypos);
  profileX_int->Fill(Xpos);
  profileY_int->Fill(Ypos);
  energy_spectrum->Fill(PM_sum);
}

//plot styling and arrangement in dedicated canvases
floodMap->GetXaxis()->SetTitle("Relative X position");
floodMap->GetYaxis()->SetTitle("Relative Y position");
floodMap->GetZaxis()->SetTitle("Entries");

TCanvas *flood_canvas = new TCanvas("flood_canvas","flood_canvas ", 600, 500);
flood_canvas->SetFillColor(0); //
flood_canvas->SetBorderMode(0);
flood_canvas->SetLeftMargin(0.1409396); //
flood_canvas->SetRightMargin(0.14865772); //
gStyle->SetOptStat(000); //
floodMap->SetTitle(NULL);
floodMap->Draw("COLZ");


profileX_int->GetXaxis()->SetTitle("Relative position");
profileX_int->GetYaxis()->SetTitle("Entries");
profileX_int->SetLineColor(kBlue);

profileY_int->GetXaxis()->SetTitle("Relative position");
profileY_int->GetYaxis()->SetTitle("Entries");
profileY_int->SetLineColor(kRed);

TCanvas *profile_canvas = new TCanvas("profile_canvas","profile_canvas ", 1000, 400);
profile_canvas->SetFillColor(0); //
profile_canvas->SetBorderMode(0);
profile_canvas->SetLeftMargin(0.1409396); //
profile_canvas->SetRightMargin(0.14865772); //
gStyle->SetOptStat(000); //
profile_canvas->Divide(2,1);
profile_canvas->cd(1);
profileX_int->Draw();
profile_canvas->cd(2);
profileY_int->Draw();

energy_spectrum->GetXaxis()->SetTitle("ADC channel");
energy_spectrum->GetYaxis()->SetTitle("Entries");

TCanvas *Espectrum_canvas = new TCanvas("Espectrum_canvas","Espectrum_canvas ", 600, 400);
Espectrum_canvas->SetFillColor(0); //
Espectrum_canvas->SetBorderMode(0);
Espectrum_canvas->SetLeftMargin(0.1409396); //
Espectrum_canvas->SetRightMargin(0.14865772); //
gStyle->SetOptStat(000); //
energy_spectrum->SetTitle(NULL);
energy_spectrum->Draw();

//OUTPUT ROOT FILE -> save plot to file, stored in the same folder as the data
TString rootfilename = "../data/block_"+f+"/RUN"+s+"/output_RUN" + s +".root";//modify according to the data path or to the desired output path
//TString rootfilename = "../data/Scan_Line_22steps_2mm_200kEv/RUN"+s+"/output_RUN" + s +".root";
TFile root_out(rootfilename, "RECREATE");
root_out.cd();

PM_canvas->Write("Raw_PMAmp");
flood_canvas->Write("Raw_FLOODMAP");
profile_canvas->Write("Raw_profileXY_int");
Espectrum_canvas->Write("Raw_energy_spectrum");

/*delete PM0_histo;
delete PM1_histo;
delete PM2_histo;
delete PM3_histo;
delete energy_spectrum;
delete profileY_int;
delete profileX_int;
delete floodMap;*/
//delete PM_canvas;
//delete flood_canvas;
//delete profile_canvas;
//delete Espectrum_canvas;
//delete leg_PM;


return;

}
