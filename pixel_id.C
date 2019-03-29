//This macro analyses the calibrated BGO data (CO60 data after PM gain equalization), and automatically extracts the position of each pseudo-pixels
//It produces the map which will be used for the pixel-based energy calibration, final analysis step

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
#include "TSpectrum2.h"
#include "TSpectrum.h"
#include "TLine.h"
#include "functions.hh"

void pixel_id(Int_t blockNum, Int_t runNumber){

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
  //Definition of the data file paths and names
  std::ostringstream ss;
  ss << std::setfill('0') << std::setw(5) << runNumber;
  string s = ss.str();
  std::ostringstream ff;
  ff << blockNum;
  string f = ff.str();
  //TString rootfilename = "../data/pixelID_tests/RUN"+s+"/outputCalibration_RUN" + s +".root";
  //string outpixels_name = "../data/pixelID_tests/RUN"+s+"/pixels_RUN" + s +".txt";
  TString rootfilename = "../data/block_"+f+"/RUN"+s+"/outputCalibration_RUN" + s +".root";//name of the root file where the calibrated 2D map is stored
  string outpixels_name = "../data/block_"+f+"/conf_files/confPIX.txt";//text file where the position of each pseudopixel is stored
  string outcolxy_linxy_name = "../data/block_"+f+"/conf_files/colxy_linxy.txt";//text file where the positions of the avarage points between pixels are stored -> used to create the map for the event assignment to pixels
  TH2D *map_histo = nullptr; //nullptr is used to initialize a pointer to a null variable
  TCanvas *flood_canvas_old = nullptr;
  TFile calib_data(rootfilename, "READ");//opening calibrated data file
  if(calib_data.IsZombie()){cout<<"Problem accessing data file"<<endl;}
  else{
    //calib_data.ls();
    //calibrated data file opened, the flood map with the 2D event positions is recovered
    if(calib_data.GetListOfKeys()->Contains("Cal_FLOODMAP")){
      flood_canvas_old = (TCanvas*)calib_data.Get("Cal_FLOODMAP");//recover the canvas where the TH2D is stored
      cout<<"Calibrated FLOOD MAP found!"<<endl;
    }else{
      cout<<"Calibrated FLOOD MAP cannot be found in file "<<rootfilename<<endl;
      return;
    }
    map_histo = (TH2D*)flood_canvas_old->GetPrimitive("floodMap");//from the canvas, the flood map tH2D is recovered
  }
  calib_data.Close();
  cout<<"Stored data retrieved! "<<endl;
  cout<<"----------------ooooooooooooooooooo---------------------------"<<endl;

  //Histos for 1D integrated profiles
  Int_t NbinsX = 0;
  Int_t NbinsY = 0;
  NbinsX = map_histo->GetXaxis()->GetNbins();
  NbinsY = map_histo->GetYaxis()->GetNbins();
  Double_t x_temp = 0., y_temp = 0.;
  TH1F *profileX_int = new TH1F("profile X int cal", "profile X int cal", NbinsX, -1, 1);
  TH1F *profileY_int = new TH1F("profile Y int cal", "profile Y int cal", NbinsY, -1, 1);
  TGraph *pixels = new TGraph(64);//maximum number of points -> 2D graph with the pseudo-pixel positions

  //variables for the pixel-search algorithm
  Int_t nfoundX = 0;
  Int_t nfoundY = 0;
  vector <Double_t> x_pos_peaks, y_pos_peaks;
  vector <Double_t> pix_Xborders;
  vector <Double_t> pix_Yborders;
  vector <vector <Double_t>> x_pos_multiPeaks, y_pos_multiPeaks;
  vector <Int_t> nfoundX_single, nfoundY_single;
  Int_t check_X = 0;
  Int_t check_Y = 0;
  vector <vector <Int_t>> missedX, missedY;
  Int_t countPNT=0;
  Int_t count_row = 0;

  //functions included in the file functions.h
  integrated_peaks(profileX_int, profileY_int, map_histo, nfoundX, nfoundY, x_pos_peaks, y_pos_peaks);
  //Now we have two vectors of peak positions in X and Y
  //The flood map can be drawn with the identified peaks/valleys
  draw_map_defBorders(map_histo, nfoundX, nfoundY, x_pos_peaks, y_pos_peaks, pix_Xborders, pix_Yborders);
  //Now the valleys have been defined, so we can move to the real pixel identification algorithm.
  //Pixel lines and columns are now treated independently
  //PIXEL IDENTIFICATION
  cout<<"----------------ooooooooooooooooooo---------------------------"<<endl;
  cout<<"Start pixel identification process!"<<endl;
  cout<<"----------------ooooooooooooooooooo---------------------------"<<endl;
  cout<<"Working on "<<nfoundY<<" rows and "<<nfoundX<<" columns "<<endl;
  //Till here the code has identified the number of rows and columns (hopefully 8 and 8), and the valleys between rows and columns
  //histograms for single pixel rows and columns
  single_rowCol_def(map_histo, nfoundX, nfoundY, pix_Xborders, pix_Yborders, nfoundX_single, nfoundY_single, x_pos_multiPeaks, y_pos_multiPeaks);
  //after this function we have the position of the pixels (peaks) for each row and each column, so we have one vector per line and on vector per column with a number of elements
  //defined in nfoundY_single and nfoundX_single
  cout<<"Check single rows!"<<endl;
  check_X = check_borders(nfoundX, nfoundY_single, pix_Yborders, y_pos_multiPeaks);
  cout<<"Check single columns!"<<endl;
  check_Y = check_borders(nfoundY, nfoundX_single, pix_Xborders, x_pos_multiPeaks);
  //after the application of this function for the two axis, the number of row and columns is defined according to the single rwo/column analysis and checked
  if(check_Y != 0) nfoundX = check_Y;//the number of rows and columns could have to be changed after the check
  if(check_X != 0) nfoundY = check_X;
  cout<<"Defined number of rows "<<nfoundY<<endl;
  cout<<"Defined number of columns "<<nfoundX<<endl;

  if(check_X || check_Y){//in case a modification is necessary, also the 2D map borders and the definition of single rows and columns must be adapted, so that the same fuctions as before (or slightly modified) are applied with the new values
    cout<<"Reevaluation of rows and columns limits "<<endl;
    map_new_borders(map_histo, nfoundX, nfoundY, pix_Xborders, pix_Yborders);
    x_pos_multiPeaks.clear();
    y_pos_multiPeaks.clear();
    single_rowCol_def(map_histo, nfoundX, nfoundY, pix_Xborders, pix_Yborders, nfoundX_single, nfoundY_single, x_pos_multiPeaks, y_pos_multiPeaks);
  }
  //NOW ALL PEAKS ARE DEFINED
  //All peaks for single rows and columns have been identified
  //The algorithm requires a squared matrix of pixels, with size as the maximum number of peaks found in one of the two directions - the peak vectors are then filled with non-sense values (25) in the position of missed pixels
  for(UInt_t Nrr = 0; Nrr<x_pos_multiPeaks.size(); Nrr++){
    if(x_pos_multiPeaks[Nrr].size()>0 && x_pos_multiPeaks[Nrr].size()<(UInt_t)nfoundX){
      cout<<"Filling row "<<Nrr<<endl;
      missedX.push_back({(Int_t)Nrr});
      complete_peaks(nfoundX, x_pos_multiPeaks[Nrr], pix_Xborders, missedX);//complete and print vector for X
      print_vector(x_pos_multiPeaks[Nrr]);
    }
  }
  for(UInt_t Ncc = 0; Ncc<y_pos_multiPeaks.size(); Ncc++){
    if(y_pos_multiPeaks[Ncc].size()>0 && y_pos_multiPeaks[Ncc].size()<(UInt_t)nfoundY){
      cout<<"Filling column "<<Ncc<<endl;
      missedY.push_back({(Int_t)Ncc});
      complete_peaks(nfoundY, y_pos_multiPeaks[Ncc], pix_Yborders, missedY);//complete and print vector for Y
      print_vector(y_pos_multiPeaks[Ncc]);
    }
  }
//now I have a complete vector for each line and column with all the pseudo-pixel positions
  cout<<"PRINT VECTORS FOR VERIFICATION!"<<endl;
  for(UInt_t priX = 0; priX<x_pos_multiPeaks.size(); priX++){
    print_vector(x_pos_multiPeaks[priX]);
  }
  for(UInt_t priY = 0; priY<y_pos_multiPeaks.size(); priY++){
    print_vector(y_pos_multiPeaks[priY]);
  }
  //Till here, all the analysis is monodimensional: we found the position of pixels in rows and columns
  //now the X and Y positions must be associated to create 2D positions in the 2D map
  cout<<"----------------ooooooooooooooooooo---------------------------"<<endl;
  cout<<"Start search for peak points in 2D!"<<endl;
  cout<<"----------------ooooooooooooooooooo---------------------------"<<endl;
  //Rows and columns are here defined and the peaks positions are found in 2D
  //Rows and columns peaks must be coupled to find peak points in 2 dimensions
  for(UInt_t Xpix = 0; Xpix<x_pos_multiPeaks.size(); Xpix++){
    fill_pixels(nfoundY, nfoundX, Xpix, x_pos_multiPeaks[Xpix], y_pos_multiPeaks, pixels, countPNT, pix_Xborders, pix_Yborders);
  }
  //remove points not filled, and by default put in (0;0)
  Double_t pntX = 0., pntY = 0.;
  Int_t count_missed = 0;
  //the missing points are counted and removed to avoid graphical artifacts
  for(Int_t pix_pnt = 0; pix_pnt<64; pix_pnt++){
    pixels->GetPoint(pix_pnt, pntX, pntY);
    if(pntX==0 && pntY == 0){
      pixels->RemovePoint(pix_pnt);
      count_missed++;
    }
  }
  cout<<"Total missed pixels from an 8x8 matrix: "<<count_missed<<endl;
  //start filling output text file with the identified pixel positions
  ofstream out_pixels(outpixels_name);
  out_pixels<<"PIXEL L"<<"\t"<<"PIXEL C"<<"\t"<<"X"<<"\t"<<"Y"<<endl;
  Int_t temp_miss_pix = 0;
  Int_t temp_miss_piy = 0;
  vector <vector<Int_t>> missed_pixels;
  //check and identify missing pixels along X
  Int_t knt = 0;
  for(UInt_t mspx = 1; mspx<missedX.size(); mspx++){
    if(missedX[mspx].size()!=0){
      cout<<"vector X "<<missedX[knt][0]<<" has "<<missedX[mspx].size()<<" entries "<<endl;
      for(UInt_t mpix = 0; mpix<missedX[mspx].size(); mpix++){
        temp_miss_pix = missedX[mspx].at(mpix);
        cout<<" identified missing pixel "<<missedX[knt][0]<<"\t"<<temp_miss_pix<<endl;
        missed_pixels.push_back({missedX[knt][0], temp_miss_pix});
      }
    }
    knt+=2;
    mspx++;
  }
  knt = 0;
  //check and identify missing pixels along Y
  for(UInt_t mspy = 1; mspy<missedY.size(); mspy++){
    if(missedY[mspy].size()!=0){
      cout<<"vector Y "<<missedY[knt][0]<<" has "<<missedY[mspy].size()<<" entries "<<endl;
      for(UInt_t mpiy = 0; mpiy<missedY[mspy].size(); mpiy++){
        temp_miss_piy = missedY[mspy].at(mpiy);
        cout<<" identified missing pixel "<<temp_miss_piy<<"\t"<<missedY[knt][0]<<endl;
        missed_pixels.push_back({temp_miss_piy, missedY[knt][0]});
      }
    }
    knt+=2;
    mspy++;
  }
  //also the missing pixels have been identified, so that the map can be manually completed in case it is needed
  if(missed_pixels.size()>0)cout<<"The missing pixels are: "<<endl;
  for(UInt_t mss = 0; mss<missed_pixels.size(); mss++){
    print_vector_int(missed_pixels[mss]);
  }
  //Check if missing pixels have been identified more than one time
  vector <Int_t> help;
  for(UInt_t ck = 0; ck<missed_pixels.size(); ck++)
    {
      help = missed_pixels.at(ck);
      for(UInt_t ck2 = ck+1; ck2<missed_pixels.size(); ck2++)
        if(help == missed_pixels.at(ck2)){
          cout<<"erased element "<<missed_pixels[ck2][0]<<"\t"<<missed_pixels[ck2][1]<<endl;
          missed_pixels.erase(missed_pixels.begin() + (Int_t)ck2);
        }
    }
  Int_t count_x=0;
  Int_t count_y=0;
  TGraph *pixels_mod = new TGraph(64-count_missed);
  for(Int_t pix_pnt_new = 0; pix_pnt_new<64-count_missed; pix_pnt_new++){
    pixels->GetPoint(pix_pnt_new, pntX, pntY);
    pixels_mod->SetPoint(pix_pnt_new, pntX, pntY);

    for(UInt_t cmp = 0; cmp<missed_pixels.size();cmp++){
      if(missed_pixels[cmp][0] == count_x && missed_pixels[cmp][1]==count_y){
        out_pixels<<count_x<<"\t"<<count_y<<"\t"<<0.<<"\t"<<0.<<endl;
        if(count_y<7){count_y++;}
        else{count_x++; count_y = 0;}
      }
    }
        out_pixels<<count_x<<"\t"<<count_y<<"\t"<<pntX<<"\t"<<pntY<<endl;
        if(count_y<7){count_y++;}
        else{count_x++; count_y = 0;}
    }
  pixels_mod->SetMarkerSize(2);
  pixels_mod->SetMarkerColor(kGreen);
  pixels_mod->SetMarkerStyle(5);
  TLine *line[nfoundX+nfoundY - 2];

  //graphical arrangement of the results
  TCanvas *map_withPixels = new TCanvas("map with pixels", "map with pixels", 600, 500);
  map_withPixels->SetFillColor(0); //
  map_withPixels->SetBorderMode(0);
  map_withPixels->SetLeftMargin(0.1409396); //
  map_withPixels->SetRightMargin(0.14865772); //
  gStyle->SetOptStat(000); //
  map_withPixels->SetTitle(NULL);
  gStyle->SetPalette(55);
  map_histo->Draw("colz");
  pixels_mod->Draw("P same");

  //once all the pixels have been defined, the map is optimized with the following function to follow the geometrical distortions given by the block light response
  new_grid(map_withPixels, pixels, outcolxy_linxy_name, nfoundX, nfoundY, count_missed, pntX, pntY);

  if (sqrt(64-count_missed)-floor(sqrt(64-count_missed))!=0 || nfoundX<6){
     for(Int_t lx = 0; lx<nfoundX-1; lx++){//vertical lines
        line[lx] = new TLine((x_pos_peaks[lx]+x_pos_peaks[lx+1])/2, -1., (x_pos_peaks[lx]+x_pos_peaks[lx+1])/2,1.);
        line[lx]->SetLineColor(kYellow);
        line[lx]->SetLineWidth(2);
        map_withPixels->cd();
        line[lx]->Draw();
     }
     for(Int_t ly = 0; ly<nfoundY-1; ly++){//horizontal lines
        line[nfoundX-1+ly] =new  TLine(-1., (y_pos_peaks[ly]+y_pos_peaks[ly+1])/2, 1., (y_pos_peaks[ly]+y_pos_peaks[ly+1])/2);
        line[nfoundX-1+ly]->SetLineColor(kYellow);
        line[nfoundX-1+ly]->SetLineWidth(2);
        map_withPixels->cd();
        line[nfoundX-1+ly]->Draw();
     }
  }

  map_withPixels->Update();//the 2D map with the pixel definition is updated with the last results and the polylines determined in the function new_grid


  if(count_missed>0)cout<<count_missed<<" pixels cannot be found automatically. Look at the map to identify their position!"<<endl;
  if(count_missed>(Int_t)missed_pixels.size()&&missed_pixels.size()>0)cout<<"Even in the identified pixel matrix "<<missed_pixels.size()<<" pixels are not identified"<<endl;
  for(UInt_t chk = 0; chk<missed_pixels.size();chk++){
    cout<<"Pixel " <<chk<<"\t"<<missed_pixels[chk][0]<<"\t"<<missed_pixels[chk][1]<<endl;
  }
  if(count_missed>0)cout<<"Check the plot for the positions of the other pixels!"<<endl;
  cout<<"All position saved in "<<outpixels_name<<endl;
  cout<<"DONE! Bye!"<<endl;
  out_pixels.close();

}
