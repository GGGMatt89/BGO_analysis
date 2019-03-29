#include "TSpectrum2.h"
#include "TSpectrum.h"
#include "TLine.h"

void print_vector(vector <Double_t> peaks){//prints the content of a vector of double in a line
  if(peaks.size()!=0){
    for(UInt_t pr = 0; pr<peaks.size()-1; pr++){
      cout<<peaks[pr]<<"\t";
    }
    cout<<peaks[peaks.size()-1]<<endl;
  }
}
//---------------------------------------000000000000000000----------------------------------------//
void print_vector_int(vector <Int_t> peaks){//prints the content of a vector of int in a line
  if(peaks.size()!=0){
    for(UInt_t pr = 0; pr<peaks.size()-1; pr++){
      cout<<peaks[pr]<<"\t";
    }
    cout<<peaks[peaks.size()-1]<<endl;
  }
}
//---------------------------------------000000000000000000----------------------------------------//
void plot_to_vec(TH2D *map_histo, vector <Double_t> &x_pos, vector <Double_t> &y_pos){

  //this function extracts the mono-dimensional profiles (integrated) from the 2-dim flood map, and store the values in two vectors for the x and y axes

  //Definition of the total number of bins in the 2D flood map
  Int_t NbinsX = 0;
  Int_t NbinsY = 0;
  NbinsX = map_histo->GetXaxis()->GetNbins();
  NbinsY = map_histo->GetYaxis()->GetNbins();

  //Fill the 1D X and Y profile from the flood map data
  Double_t x_temp = 0., y_temp = 0.;
  Int_t count = 0;
  for(Int_t xx = 1; xx<=NbinsX; xx++){
    for(Int_t yy = 1; yy<=NbinsY; yy++){
      x_temp = map_histo->GetXaxis()->GetBinCenter(xx);
      y_temp = map_histo->GetYaxis()->GetBinCenter(yy);
      count = map_histo->GetBinContent(xx, yy);
      for(Int_t pp = 0; pp<count; pp++){
        x_pos.push_back(x_temp);
        y_pos.push_back(y_temp);
      }
    }
  }
}
//---------------------------------------000000000000000000----------------------------------------//
void integrated_peaks(TH1F *profileX_int, TH1F *profileY_int, TH2D *map_histo,Int_t &nfoundX, Int_t &nfoundY, vector <Double_t> &x_pos_peaks, vector <Double_t> &y_pos_peaks){

  //this function uses the 2D flood map to create the 1D integrated spatial distributions along X and Y
  //once the two 1D distributions are obtained, it searches the peaks (corresponding to the average center of the pixel lines and columns), and records the number of peaks found
  //per axis (nfoundX and nfoundY) and their positions (x_pos_peaks and y_pos_peaks) -> all these are passed to the function by reference, so that they can be modified by the fuction

    //Definition of the total number of bins in the 2D flood map
    Int_t NbinsX = 0;
    Int_t NbinsY = 0;
    NbinsX = map_histo->GetXaxis()->GetNbins();
    NbinsY = map_histo->GetYaxis()->GetNbins();

    //Fill the 1D X and Y profile from the flood map data
    Double_t x_temp = 0., y_temp = 0.;
    for(Int_t xx = 1; xx<=NbinsX; xx++){
      for(Int_t yy = 1; yy<=NbinsY; yy++){
        x_temp = map_histo->GetXaxis()->GetBinCenter(xx);
        y_temp = map_histo->GetYaxis()->GetBinCenter(yy);
        profileX_int->Fill(x_temp,map_histo->GetBinContent(xx, yy));
        profileY_int->Fill(y_temp,map_histo->GetBinContent(xx, yy));
      }
    }
    //1D profiles plot styling
    profileX_int->GetXaxis()->SetTitle("Relative position");
    profileX_int->GetYaxis()->SetTitle("Entries");
    profileX_int->SetLineColor(kBlue);
    profileX_int->SetStats(kFALSE);
    profileY_int->GetXaxis()->SetTitle("Relative position");
    profileY_int->GetYaxis()->SetTitle("Entries");
    profileY_int->SetLineColor(kRed);
    profileY_int->SetStats(kFALSE);
    //Canvas for the 1D X and Y profiles before the peak identification
    TCanvas *profile_canvas = new TCanvas("profile_canvas_cal","profile_canvas_cal", 1000, 400);
    profile_canvas->SetFillColor(0); //
    profile_canvas->SetBorderMode(0);
    profile_canvas->SetLeftMargin(0.1409396); //
    profile_canvas->SetRightMargin(0.14865772); //
    profile_canvas->Divide(2,1);
    profile_canvas->cd(1);
    profileX_int->Draw("histo");
    profile_canvas->cd(2);
    profileY_int->Draw("histo");

    //searching for peaks in 1D profiles
    //the object TSpectrum is used to analyse a profile and extract the desired information -> for example the peaks
    TSpectrum *spectX = new TSpectrum(8); //8 is the maximum number of peaks in the distribution -> 8 pixels per line/column
    TSpectrum *spectY = new TSpectrum(8);
    //Canvas for the 1D X and Y profiles after the peak identification -> chech for peak position
    TCanvas *profile_canvas_peaks = new TCanvas("profile_canvas_withpeaks","profile_canvas_withpeaks", 1000, 400);
    profile_canvas_peaks->SetFillColor(0); //
    profile_canvas_peaks->SetBorderMode(0);
    profile_canvas_peaks->SetLeftMargin(0.1409396); //
    profile_canvas_peaks->SetRightMargin(0.14865772); //
    profile_canvas_peaks->Divide(2,1);
    profile_canvas_peaks->cd(1);
    nfoundX = spectX->Search(profileX_int);//function of an object TSpectrum which identifies the peaks
    profile_canvas_peaks->cd(2);
    nfoundY = spectY->Search(profileY_int);
    //Copying peak positions in two vectors
    //xpeaks is an array, resulting from the function GetPositionX applied on spectX
    Double_t *xpeaks = spectX->GetPositionX();
    for (Int_t p=0;p<nfoundX;p++) {
      Double_t xp = xpeaks[p];
      x_pos_peaks.push_back(xp);
    }
    Double_t *ypeaks = spectY->GetPositionX();
      for (Int_t pp=0;pp<nfoundY;pp++) {
        Double_t yp = ypeaks[pp];
        y_pos_peaks.push_back(yp);
      }
    //the TSpectrum object does not list the peaks in the correct order, but in increasing weight order -> for us it is necessary to have their position in increasing order
    //The two vectors have to be sorted in increasing order
    std::sort (x_pos_peaks.begin(), x_pos_peaks.end());
    std::sort (y_pos_peaks.begin(), y_pos_peaks.end());
    cout<<nfoundX<<" peaks in X distribution : ";
    print_vector(x_pos_peaks);//function defined in this file (functions.h - top of the text file), just used to print the obtained peak positions
    cout<<nfoundY<<" peaks of Y distribution : ";
    print_vector(y_pos_peaks);
}
//---------------------------------------000000000000000000----------------------------------------//
void draw_map_defBorders(TH2D *map_histo, Int_t nfoundX, Int_t nfoundY, vector <Double_t> x_pos_peaks, vector <Double_t> y_pos_peaks, vector <Double_t> &pix_Xborders, vector <Double_t> &pix_Yborders){

  //this function takes the 2D flood map and draws it with the addition of the horizontal and vertical lines corresponding to the valleys between pixel xpeaks_single
  //the valley positions are, at this stage, calculated as the simple average between two following peak positions
  //the position of the valleys are also stored in two vectors, pix_Xborders and pix_Yborders

  TCanvas *flood_canvas = new TCanvas("flood_canvas","flood_canvas ", 600, 500);
  flood_canvas->SetFillColor(0); //
  flood_canvas->SetBorderMode(0);
  flood_canvas->SetLeftMargin(0.1409396); //
  flood_canvas->SetRightMargin(0.14865772); //
  gStyle->SetOptStat(000); //
  map_histo->SetTitle(NULL);
  gStyle->SetPalette(55);
  map_histo->Draw("colz");

  //Definition of valley lines for the flood map plot
  TLine *line[nfoundX+nfoundY - 2];
  for(Int_t lx = 0; lx<nfoundX-1; lx++){//vertical lines
    pix_Xborders.push_back((x_pos_peaks[lx]+x_pos_peaks[lx+1])/2);//average position between two next peaks
    line[lx] = new TLine(pix_Xborders[lx], -1., pix_Xborders[lx],1.);//definition of the line with (x initial pos, y initial pos, x final pos, y final pos) -> if x initial and final positions are the same, the line is vertical
    line[lx]->SetLineColor(kYellow);
    line[lx]->SetLineWidth(2);
    flood_canvas->cd();
    line[lx]->Draw();
  }
  for(Int_t ly = 0; ly<nfoundY-1; ly++){//horizontal lines
    pix_Yborders.push_back((y_pos_peaks[ly]+y_pos_peaks[ly+1])/2);
    line[nfoundX-1+ly] =new  TLine(-1., pix_Yborders[ly], 1., pix_Yborders[ly]);
    line[nfoundX-1+ly]->SetLineColor(kYellow);
    line[nfoundX-1+ly]->SetLineWidth(2);
    flood_canvas->cd();
    line[nfoundX-1+ly]->Draw();
  }
  flood_canvas->Update();
  cout<<"Identified column limits : ";
  print_vector(pix_Xborders);
  cout<<"Identified row limits : ";
  print_vector(pix_Yborders);
}
//---------------------------------------000000000000000000----------------------------------------//

//---------------------------------------000000000000000000----------------------------------------//
//void single_rowCol_def(TH2D *map_histo, Int_t nfoundX, Int_t nfoundY, vector <Double_t> pix_Xborders, vector <Double_t> pix_Yborders, vector <Int_t> &nfoundX_single, vector <Int_t> &nfoundY_single, vector <Double_t> &x_pos_peaks_0, vector <Double_t> &x_pos_peaks_1, vector <Double_t> &x_pos_peaks_2, vector <Double_t> &x_pos_peaks_3, vector <Double_t> &x_pos_peaks_4, vector <Double_t> &x_pos_peaks_5, vector <Double_t> &x_pos_peaks_6, vector <Double_t> &x_pos_peaks_7, vector <Double_t> &y_pos_peaks_0, vector <Double_t> &y_pos_peaks_1, vector <Double_t> &y_pos_peaks_2, vector <Double_t> &y_pos_peaks_3, vector <Double_t> &y_pos_peaks_4, vector <Double_t> &y_pos_peaks_5, vector <Double_t> &y_pos_peaks_6, vector <Double_t> &y_pos_peaks_7){

void single_rowCol_def(TH2D *map_histo, Int_t nfoundX, Int_t nfoundY, vector <Double_t> pix_Xborders, vector <Double_t> pix_Yborders, vector <Int_t> &nfoundX_single, vector <Int_t> &nfoundY_single, vector <vector <Double_t>> &x_pos_multiPeaks, vector <vector <Double_t>> &y_pos_multiPeaks){

  //this function uses the number of rows and columns and the borders between rows and columns to define the single row and column mono-dimensional distributions,
  //creates the histograms for each row and column and finds the peaks (pixel positions) with the TSpectrum objects for each row and column

  //Definition of histos for single line and column distributions
  TH1F *rows[nfoundX];
  for(Int_t defR = 0; defR<nfoundY; defR++){
    rows[defR] = new TH1F(Form("profile X single line %d", defR), Form("profile X single line %d", defR), 200, -1, 1);
  }
  TH1F *columns[nfoundY];
  for(Int_t defC = 0; defC<nfoundX; defC++){
    columns[defC] = new TH1F(Form("profile Y single column %d", defC), Form("profile Y single column %d", defC), 200, -1, 1);
  }

  Double_t x_temp = 0., y_temp = 0.;
  Int_t row_count = 0;
  Int_t column_count = 0;
  Int_t histo_N = 0;
  //Filling single row histograms
  Int_t NbinsX = 0;
  Int_t NbinsY = 0;
  NbinsX = map_histo->GetXaxis()->GetNbins();
  NbinsY = map_histo->GetYaxis()->GetNbins();
  for(Int_t yy = 1; yy<=NbinsY; yy++){
    if(map_histo->GetYaxis()->GetBinCenter(yy)>pix_Yborders[row_count]&&row_count<(nfoundY-1)){
      row_count++;
    }
    for(Int_t xx = 1; xx<=NbinsX; xx++){
      x_temp = map_histo->GetXaxis()->GetBinCenter(xx);
      rows[row_count]->Fill(x_temp,map_histo->GetBinContent(xx, yy));
    }
  }
  //Filling single column histograms
  for(Int_t xxc = 1; xxc<=NbinsX; xxc++){
    if(map_histo->GetXaxis()->GetBinCenter(xxc)>pix_Xborders[column_count]&&column_count<(nfoundX-1)){
      column_count++;
    }
    for(Int_t yyc = 1; yyc<=NbinsY; yyc++){
      y_temp = map_histo->GetYaxis()->GetBinCenter(yyc);
      columns[column_count]->Fill(y_temp,map_histo->GetBinContent(xxc, yyc));
    }
  }
  cout<<"Single pixel rows and columns defined!"<<endl;
  //Canvas for single rows and columns
  TCanvas *single_lines = new TCanvas("single lines canvas", "single lines canvas", 1400, 1400);
  single_lines->SetFillColor(0); //
  single_lines->SetBorderMode(0);
  single_lines->SetLeftMargin(0.1409396); //
  single_lines->SetRightMargin(0.14865772); //
  single_lines->Divide(4,2);
  TCanvas *single_columns = new TCanvas("single columns canvas", "single columns canvas", 1400, 1400);
  single_columns->SetFillColor(0); //
  single_columns->SetBorderMode(0);
  single_columns->SetLeftMargin(0.1409396); //
  single_columns->SetRightMargin(0.14865772); //
  single_columns->Divide(4,2);
  //Fill canvas for single rows and columns
  for(Int_t pltR = 0; pltR<nfoundY; pltR++){
    single_lines->cd(pltR+1);
    if(!rows[pltR]->IsZombie()){
      rows[pltR]->SetLineColor(kBlue);
      rows[pltR]->GetXaxis()->SetTitle("Relative position");
      rows[pltR]->GetYaxis()->SetTitle("Entries");
      if(rows[pltR]->GetEntries()!=0)rows[pltR]->Draw("histo");
    }
  }
  for(Int_t pltC = 0; pltC<nfoundX; pltC++){
    single_columns->cd(pltC+1);
    if(!columns[pltC]->IsZombie()){
      columns[pltC]->SetLineColor(kBlue);
      columns[pltC]->GetXaxis()->SetTitle("Relative position");
      columns[pltC]->GetYaxis()->SetTitle("Entries");
    if(columns[pltC]->GetEntries()!=0)columns[pltC]->Draw("histo");
      }
    }
    //Single rows and columns distributions have been defined and drawn
    cout<<"----------------ooooooooooooooooooo---------------------------"<<endl;
    cout<<"Start finding peaks on single rows and column distributions!"<<endl;
    cout<<"----------------ooooooooooooooooooo---------------------------"<<endl;
    TSpectrum *spectX_singleLine[nfoundY];//objects to define the pixel positions in each row and column
    TSpectrum *spectY_singleLine[nfoundX];
    //Canvas for single rows and columns with the identified peaks
    TCanvas *single_lines_withPeak = new TCanvas("single lines canvas with peaks", "single lines canvas with peaks", 1400, 1400);
    single_lines_withPeak->SetFillColor(0); //
    single_lines_withPeak->SetBorderMode(0);
    single_lines_withPeak->SetLeftMargin(0.1409396); //
    single_lines_withPeak->SetRightMargin(0.14865772); //
    single_lines_withPeak->Divide(4,2);
    for(Int_t slpx = 0; slpx<nfoundY; slpx++){
      single_lines_withPeak->cd(slpx+1);
      spectX_singleLine[slpx] =  new TSpectrum(8);
      nfoundX_single.push_back(spectX_singleLine[slpx]->Search(rows[slpx], 2));//vector of arrays, one per column
    }
    TCanvas *single_columns_withPeak = new TCanvas("single columns canvas with peaks", "single columns canvas with peaks", 1400, 1400);
    single_columns_withPeak->SetFillColor(0); //
    single_columns_withPeak->SetBorderMode(0);
    single_columns_withPeak->SetLeftMargin(0.1409396); //
    single_columns_withPeak->SetRightMargin(0.14865772); //
    single_columns_withPeak->Divide(4,2);
    for(Int_t scpx = 0; scpx<nfoundX; scpx++){
      single_columns_withPeak->cd(scpx+1);
      spectY_singleLine[scpx] =  new TSpectrum(8);
      nfoundY_single.push_back(spectY_singleLine[scpx]->Search(columns[scpx], 2));//vector of arrays, one per row
    }
    //Assign peaks values to a vector per row and per column in a vector of vectors
    Double_t *xpeaks_single[nfoundY];
    Double_t *ypeaks_single[nfoundX];
    Double_t xpll, ypll;
    for(Int_t llx = 0; llx<nfoundY; llx++){
      xpeaks_single[llx]= spectX_singleLine[llx]->GetPositionX();
    }
    for(Int_t lly = 0; lly<nfoundX; lly++){
      ypeaks_single[lly]= spectY_singleLine[lly]->GetPositionX();
    }
    vector <Double_t> temp_vec;
    for(Int_t llccR = 0; llccR<nfoundY; llccR++){
      temp_vec.clear();
      for (Int_t pls=0;pls<nfoundX_single[llccR];pls++) {
        Double_t xpll = xpeaks_single[llccR][pls];
          temp_vec.push_back(xpll);
        /*if(llccR == 0 && llccR<(nfoundY-1))x_pos_peaks_0.push_back(xpll);
        if(llccR == 1 && llccR<(nfoundY-1))x_pos_peaks_1.push_back(xpll);
        if(llccR == 2 && llccR<(nfoundY-1))x_pos_peaks_2.push_back(xpll);
        if(llccR == 3 && llccR<(nfoundY-1))x_pos_peaks_3.push_back(xpll);
        if(llccR == 4 && llccR<(nfoundY-1))x_pos_peaks_4.push_back(xpll);
        if(llccR == 5 && llccR<(nfoundY-1))x_pos_peaks_5.push_back(xpll);
        if(llccR == 6 && llccR<(nfoundY-1))x_pos_peaks_6.push_back(xpll);
        if(llccR == (nfoundY-1))x_pos_peaks_7.push_back(xpll);
        */
      }
      x_pos_multiPeaks.push_back(temp_vec);
    }
    for(Int_t llccC = 0; llccC<nfoundX; llccC++){
      temp_vec.clear();
      for (Int_t pcs=0;pcs<nfoundY_single[llccC];pcs++) {
        Double_t ypll = ypeaks_single[llccC][pcs];
        temp_vec.push_back(ypll);
        /*if(llccC == 0 && llccC<(nfoundX-1))y_pos_peaks_0.push_back(ypll);
        if(llccC == 1 && llccC<(nfoundX-1))y_pos_peaks_1.push_back(ypll);
        if(llccC == 2 && llccC<(nfoundX-1))y_pos_peaks_2.push_back(ypll);
        if(llccC == 3 && llccC<(nfoundX-1))y_pos_peaks_3.push_back(ypll);
        if(llccC == 4 && llccC<(nfoundX-1))y_pos_peaks_4.push_back(ypll);
        if(llccC == 5 && llccC<(nfoundX-1))y_pos_peaks_5.push_back(ypll);
        if(llccC == 6 && llccC<(nfoundX-1))y_pos_peaks_6.push_back(ypll);
        if(llccC == (nfoundX-1))y_pos_peaks_7.push_back(ypll);*/
      }
      y_pos_multiPeaks.push_back(temp_vec);
    }
    //Peak position vector must be sorted in increasing order
    cout<<"Found x positions : "<<endl;
    for(Int_t srtX = 0; srtX<nfoundY; srtX++){
      cout<<"Line "<<srtX<<" ";
      std::sort (x_pos_multiPeaks[srtX].begin(), x_pos_multiPeaks[srtX].end());
      print_vector(x_pos_multiPeaks[srtX]);
    }
    cout<<"Found y positions : "<<endl;
    for(Int_t srtY = 0; srtY<nfoundX; srtY++){
      cout<<"Column "<<srtY<<" ";
      std::sort (y_pos_multiPeaks[srtY].begin(), y_pos_multiPeaks[srtY].end());
      print_vector(y_pos_multiPeaks[srtY]);
    }
    /*
    if(x_pos_peaks_0.size()>0)std::sort (x_pos_peaks_0.begin(), x_pos_peaks_0.end());
    if(x_pos_peaks_1.size()>0)std::sort (x_pos_peaks_1.begin(), x_pos_peaks_1.end());
    if(x_pos_peaks_2.size()>0)std::sort (x_pos_peaks_2.begin(), x_pos_peaks_2.end());
    if(x_pos_peaks_3.size()>0)std::sort (x_pos_peaks_3.begin(), x_pos_peaks_3.end());
    if(x_pos_peaks_4.size()>0)std::sort (x_pos_peaks_4.begin(), x_pos_peaks_4.end());
    if(x_pos_peaks_5.size()>0)std::sort (x_pos_peaks_5.begin(), x_pos_peaks_5.end());
    if(x_pos_peaks_6.size()>0)std::sort (x_pos_peaks_6.begin(), x_pos_peaks_6.end());
    if(x_pos_peaks_7.size()>0)std::sort (x_pos_peaks_7.begin(), x_pos_peaks_7.end());
    if(y_pos_peaks_0.size()>0)std::sort (y_pos_peaks_0.begin(), y_pos_peaks_0.end());
    if(y_pos_peaks_1.size()>0)std::sort (y_pos_peaks_1.begin(), y_pos_peaks_1.end());
    if(y_pos_peaks_2.size()>0)std::sort (y_pos_peaks_2.begin(), y_pos_peaks_2.end());
    if(y_pos_peaks_3.size()>0)std::sort (y_pos_peaks_3.begin(), y_pos_peaks_3.end());
    if(y_pos_peaks_4.size()>0)std::sort (y_pos_peaks_4.begin(), y_pos_peaks_4.end());
    if(y_pos_peaks_5.size()>0)std::sort (y_pos_peaks_5.begin(), y_pos_peaks_5.end());
    if(y_pos_peaks_6.size()>0)std::sort (y_pos_peaks_6.begin(), y_pos_peaks_6.end());
    if(y_pos_peaks_7.size()>0)std::sort (y_pos_peaks_7.begin(), y_pos_peaks_7.end());
   */
}
//---------------------------------------000000000000000000----------------------------------------//
//Int_t check_borders(Int_t nfoundX, vector <Int_t> nfoundY_single, vector <Double_t> &pix_Yborders, vector <Double_t > y_pos_peaks_0, vector <Double_t > y_pos_peaks_1, vector <Double_t > y_pos_peaks_2, vector <Double_t > y_pos_peaks_3, vector <Double_t > y_pos_peaks_4, vector <Double_t > y_pos_peaks_5, vector <Double_t > y_pos_peaks_6, vector <Double_t > y_pos_peaks_7){
Int_t check_borders(Int_t nfoundX, vector <Int_t> nfoundY_single, vector <Double_t> &pix_Yborders, vector <vector <Double_t>> y_pos_multiPeaks){

  //this function verifies if the number of pixels identified with the single row/column analysis is coherent with the one found with the integrated distributions
  //if so, everything's ok
  //if not, the number of rows/columns is redefined for the further analysis -> pix_Xborders and pix_Yborders can be redefined after this function

  Int_t flag = 1;
  for(Int_t ckX = 1; ckX<nfoundX; ckX++){
    if((UInt_t)nfoundY_single[ckX]>(pix_Yborders.size()+1)&&nfoundY_single[ckX]==nfoundY_single[ckX-1]){
      flag++;
    }
  }//count how many times the number of peaks found in a row/columns exceeds the N of peaks found in the integrated distribution->
  //if this number is greater from the initial one and always the same in all rows and columns, the actual number of peaks is changed and the borders vectors are redefined
  //because this means that with the integrated distribution a pixel row/column was probably missed
  Double_t pYtemp = 0.;
  Int_t count_ck = 0;
  vector <Double_t> temp_peaks;
  if(flag == nfoundX){
    cout<<"Detected not-coherent number of peaks -> smoothing ongoing ... "<<endl;
    pix_Yborders.clear();

    for(Int_t ydef = 0; ydef<nfoundY_single[0]; ydef++){
      pYtemp=0.;
      count_ck = 0;
      for(UInt_t multi = 0; multi < y_pos_multiPeaks.size(); multi++){
        if(y_pos_multiPeaks[multi].size()>0){
          pYtemp += y_pos_multiPeaks[multi][ydef];
          count_ck++;
        }
    }
    /*  if(y_pos_peaks_1.size()>0){
        pYtemp += y_pos_peaks_1[ydef];
        count_ck++;
      }
      if(y_pos_peaks_2.size()>0){
        pYtemp += y_pos_peaks_2[ydef];
        count_ck++;
      }
      if(y_pos_peaks_3.size()>0){
        pYtemp += y_pos_peaks_3[ydef];
        count_ck++;
      }
      if(y_pos_peaks_4.size()>0){
        pYtemp += y_pos_peaks_4[ydef];
        count_ck++;
      }
      if(y_pos_peaks_5.size()>0){
        pYtemp += y_pos_peaks_5[ydef];
        count_ck++;
      }
      if(y_pos_peaks_6.size()>0){
        pYtemp += y_pos_peaks_6[ydef];
        count_ck++;
      }
      if(y_pos_peaks_7.size()>0){
        pYtemp += y_pos_peaks_7[ydef];
        count_ck++;
      }
    }*/
    pYtemp = pYtemp/count_ck;
    temp_peaks.push_back(pYtemp);
    }
    for(UInt_t pr1 = 1; pr1<temp_peaks.size(); pr1++){
      pix_Yborders.push_back((temp_peaks[pr1] + temp_peaks[pr1-1])/2);
    }
  return (Int_t)temp_peaks.size();
  }else{
    return 0;
  }
}
//---------------------------------------000000000000000000----------------------------------------//
void map_new_borders(TH2D *map_histo, Int_t nfoundX, Int_t nfoundY, vector <Double_t> pix_Xborders, vector <Double_t> pix_Yborders){
  //this function redefines the rwo and columns borders after the check performed with the single row and columns analysis
  //the 2D map with the pixel limits is updated
  TCanvas *new_flood_canvas = new TCanvas("new_flood_canvas","new_flood_canvas ", 600, 500);
  new_flood_canvas->SetFillColor(0); //
  new_flood_canvas->SetBorderMode(0);
  new_flood_canvas->SetLeftMargin(0.1409396); //
  new_flood_canvas->SetRightMargin(0.14865772); //
  gStyle->SetOptStat(000); //
  map_histo->SetTitle(NULL);
  gStyle->SetPalette(55);
  map_histo->Draw("colz");
  //Definition of valley lines for the flood map plot
  TLine *line[nfoundX+nfoundY - 2];
  for(Int_t lx = 0; lx<nfoundX-1; lx++){//vertical lines
    line[lx] = new TLine(pix_Xborders[lx], -1., pix_Xborders[lx],1.);
    line[lx]->SetLineColor(kYellow);
    line[lx]->SetLineWidth(2);
    new_flood_canvas->cd();
    line[lx]->Draw();
  }
  for(Int_t ly = 0; ly<nfoundY-1; ly++){//horizontal lines
    line[nfoundX-1+ly] =new  TLine(-1., pix_Yborders[ly], 1., pix_Yborders[ly]);
    line[nfoundX-1+ly]->SetLineColor(kYellow);
    line[nfoundX-1+ly]->SetLineWidth(2);
    new_flood_canvas->cd();
    line[nfoundX-1+ly]->Draw();
  }
  new_flood_canvas->Update();
}
//---------------------------------------000000000000000000----------------------------------------//
void complete_peaks(Int_t N, vector <Double_t> &peaks, vector <Double_t> borders, vector <vector <Int_t>> &missed){

  //this function fills a vector of pixel positions with the retrieverd positions: in case the previous analysis found a missing pixel, it fills the missing positions
  //with a nonsense value (25)
  vector <Int_t> temp;
  if(peaks.size()<(UInt_t)N){
    for(Int_t scan = 0; scan < N; scan++){
      if(scan == 0){//side point
        if(peaks[scan]<borders[scan]&&peaks[scan]>=-1.){
          //the 0 position is already fill, do nothing
        }else{
          peaks.insert(peaks.begin()+scan, 25);//if position missing, fill with value 25 and add the index to the missed points vector
          temp.push_back(scan);
        }
      }
      if(scan>0 && scan < (N-1)){//central points
        if(peaks[scan]<borders[scan]&&peaks[scan]>borders[scan-1]){
        }else{
          peaks.insert(peaks.begin()+scan, 25);
          temp.push_back(scan);
        }
      }
      if(scan==(N-1)){//side point
        if(peaks[scan]>borders[N-2]&&peaks[scan]<=1.){
        }else{
          peaks.insert(peaks.begin()+scan, 25);
          temp.push_back(scan);
        }
      }
    }
  }
  missed.push_back(temp);
  temp.clear();
}
//---------------------------------------000000000000000000----------------------------------------//
void check_and_fill(Int_t Nrow, Int_t Ncol, Int_t row, Int_t &countX, vector <Double_t> xpeaks, vector <Double_t> ypeaks, TGraph *pixels, Int_t &countPNT, vector <Double_t> borders_x, vector <Double_t> borders_y){

  //according to the number of pixel positions in 1D, this function associate the positions to create 2D points and set them in the pixels graph
  Bool_t flag = kFALSE;
  if(countX<(Ncol-1)&&row < (Nrow-1)&&xpeaks[countX]<borders_x[countX] && ypeaks[row]<borders_y[row]&&xpeaks[countX]!=25&&ypeaks[row]!=25){
    pixels->SetPoint(countPNT,xpeaks[countX], ypeaks[row]);
    cout<<"Point "<<countPNT<<" X = "<< xpeaks[countX]<<" Y = "<<ypeaks[row]<<endl;
    countX++;
    countPNT++;
    flag = !flag;
  }//standard case, pixel in a central position
  else if(countX==(Ncol-1) && row < (Nrow-1) && xpeaks[countX]>borders_x[countX-1] && ypeaks[row]<borders_y[row]&&xpeaks[countX]!=25&&ypeaks[row]!=25){
    pixels->SetPoint(countPNT,xpeaks[countX], ypeaks[row]);
    cout<<"Point "<<countPNT<<" X = "<< xpeaks[countX]<<" Y = "<<ypeaks[row]<<endl;
    countPNT++;
    flag = !flag;
  }//pixels in the last column
  else if(countX<(Ncol-1) && row == (Nrow-1) && xpeaks[countX]<borders_x[countX] && ypeaks[row]>borders_y[row-1]&&xpeaks[countX]!=25&&ypeaks[row]!=25){
    pixels->SetPoint(countPNT,xpeaks[countX], ypeaks[row]);
    cout<<"Point "<<countPNT<<" X = "<< xpeaks[countX]<<" Y = "<<ypeaks[row]<<endl;
    countX++;
    countPNT++;
    flag = !flag;
  }//pixel in the last row but not in the last column
  else if(countX==(Ncol-1) && row == (Nrow-1) && xpeaks[countX]>borders_x[countX-1] && ypeaks[row]>borders_y[row-1]&&xpeaks[countX]!=25&&ypeaks[row]!=25){
    pixels->SetPoint(countPNT,xpeaks[countX], ypeaks[row]);
    cout<<"Point "<<countPNT<<" X = "<< xpeaks[countX]<<" Y = "<<ypeaks[row]<<endl;
    countPNT++;
    flag = !flag;
  }//pixel in the last row and in the last column
  if(!flag){
    countX++;
  }//pixel not identified
}
//void fill_pixels(Int_t Nrow, Int_t Ncol, Int_t row, vector <Double_t> xpeaks, vector <Double_t> ypeaks_0, vector <Double_t> ypeaks_1,vector <Double_t> ypeaks_2,vector <Double_t> ypeaks_3,vector <Double_t> ypeaks_4,vector <Double_t> ypeaks_5,vector <Double_t> ypeaks_6,vector <Double_t> ypeaks_7, TGraph *pixels, Int_t &countPNT, vector <Double_t> borders_x, vector <Double_t> borders_y){
void fill_pixels(Int_t Nrow, Int_t Ncol, Int_t row, vector <Double_t> xpeaks, vector <vector <Double_t>> ypeaks, TGraph *pixels, Int_t &countPNT, vector <Double_t> borders_x, vector <Double_t> borders_y){
  //this function scans the pixel positions in rows and columns to couple the positions and create 2D points
  cout<<"Start defining pixels for line "<<row<<endl;
  Int_t countX = 0;
  for(UInt_t col = 0; col<ypeaks.size(); col++){
    cout<<"Filling row "<<row<<" column "<<col<<endl;
    check_and_fill(Nrow, Ncol, row, countX, xpeaks, ypeaks[col], pixels, countPNT, borders_x, borders_y);
  }
}
//---------------------------------------000000000000000000----------------------------------------//

void new_grid(TCanvas *map_withPixels, TGraph *pixels, string outcolxy_linxy_name, Int_t nfoundX, Int_t nfoundY, Int_t count_missed, Double_t pntX, Double_t pntY){//make the non-square grid and right the point needed for pixels_Ecal

//this function takes the results of the whole analysis and optimize the 2D map by creating polylines
//Polylines can be well adapted to represent the geometrical distortions due to the detector light response.
//With the polylines, the real separations between pixels can be reproduced in the 2D map

  TGraph *linpnt_mod = new TGraph(56);
  TGraph *colpnt_mod = new TGraph(56);
  Double_t Px[64-count_missed];
  Double_t Py[64-count_missed];

  for (Int_t numpix=0; numpix<64-count_missed; numpix++){
      pixels->GetPoint(numpix,pntX,pntY);
      Px[numpix]=pntX;
      Py[numpix]=pntY;
  }
  linpnt_mod->SetMarkerSize(2);
  linpnt_mod->SetMarkerColor(kRed);
  linpnt_mod->SetMarkerStyle(5);
  linpnt_mod->Draw("P same");
  colpnt_mod->SetMarkerSize(2);
  colpnt_mod->SetMarkerColor(0);
  colpnt_mod->SetMarkerStyle(5);
  colpnt_mod->Draw("P same");
  TPolyLine *lineplus[nfoundX+nfoundY-2];

  ofstream out_colxy_linxy(outcolxy_linxy_name);
 // out_colxy_linxy<< "NumPt" << "\t" << "CX" << "\t" << "CY" << "\t" << "LX" << "\t" << "LY" << endl;

  if (nfoundX==8 && sqrt(64-count_missed)-floor(sqrt(64-count_missed))==0){
     for (Int_t numli=0; numli<nfoundX-1; numli++){

         Double_t CX[10];
         Double_t CY[10];
         Double_t LX[10];
         Double_t LY[10];

         CX[0]=((((Px[numli]+Px[numli+1])/2)-((Px[numli+nfoundX]+Px[numli+nfoundX+1])/2))*(1+((Py[numli+nfoundX]+Py[numli+nfoundX+1])/2))/(((Py[numli+nfoundX]+Py[numli+nfoundX+1])/2)-(Py[numli]+Py[numli+1])/2))+((Px[numli+nfoundX]+Px[numli+nfoundX+1])/2);
         CY[0]=-1;
         LX[0]=-1;
         LY[0]=((Py[nfoundX*numli]+Py[(nfoundX*numli)+nfoundX])/2)+((((Py[(nfoundX*numli)+1]+Py[(nfoundX*numli)+nfoundX+1])/2)-((Py[nfoundX*numli]+Py[(nfoundX*numli)+nfoundX])/2))*(1+((Px[nfoundX*numli]+Px[(nfoundX*numli)+nfoundX])/2))/(((Px[nfoundX*numli]+Px[(nfoundX*numli)+nfoundX])/2)-((Px[(nfoundX*numli)+1]+Px[(nfoundX*numli)+nfoundX+1])/2)));

         for (Int_t m=1; m<9; m++){
             CX[m]=(Px[numli+((m-1)*nfoundX)]+Px[numli+((m-1)*nfoundX)+1])/2;
             CY[m]=(Py[numli+((m-1)*nfoundX)]+Py[numli+((m-1)*nfoundX)+1])/2;
             LX[m]=(Px[(nfoundX*numli)+m-1]+Px[(nfoundX*numli)+nfoundX+m-1])/2;
             LY[m]=(Py[(nfoundX*numli)+m-1]+Py[(nfoundX*numli)+nfoundX+m-1])/2;
             out_colxy_linxy<< (numli*8)+m-1 << "\t" << CX[m] << "\t" << CY[m] << "\t" << LX[m] << "\t" << LY[m] << endl;
             linpnt_mod->SetPoint((numli*8)+m-1,LX[m],LY[m]);
             colpnt_mod->SetPoint((numli*8)+m-1,CX[m],CY[m]);
         }

         CX[9]=((((Px[numli+(6*nfoundX)]+Px[numli+(6*nfoundX)+1])/2)-((Px[numli+(7*nfoundX)]+Px[numli+(7*nfoundX)+1])/2))*(1-((Py[numli+(6*nfoundX)]+Py[numli+(6*nfoundX)+1])/2))/(((Py[numli+(6*nfoundX)]+Py[numli+(6*nfoundX)+1])/2)-(Py[numli+(7*nfoundX)]+Py[numli+(7*nfoundX)+1])/2))+((Px[numli+(6*nfoundX)]+Px[numli+(6*nfoundX)+1])/2);
         CY[9]=1;
         LX[9]=1;
         LY[9]=((Py[(nfoundX*numli)+6]+Py[(nfoundX*numli)+nfoundX+6])/2)+((((Py[(nfoundX*numli)+6]+Py[(nfoundX*numli)+nfoundX+6])/2)-((Py[(nfoundX*numli)+7]+Py[(nfoundX*numli)+nfoundX+7])/2))*(1-((Px[(nfoundX*numli)+6]+Px[(nfoundX*numli)+nfoundX+6])/2))/(((Px[(nfoundX*numli)+6]+Px[(nfoundX*numli)+nfoundX+6])/2)-((Px[(nfoundX*numli)+7]+Px[(nfoundX*numli)+nfoundX+7])/2)));

         lineplus[numli] = new TPolyLine(nfoundX+2,CX,CY);
         lineplus[numli]->SetLineColor(kYellow);
         lineplus[numli]->SetLineWidth(1);
         map_withPixels->cd();
         lineplus[numli]->Draw();
         lineplus[numli+nfoundX-1] = new TPolyLine(nfoundX+2,LX,LY);
         lineplus[numli+nfoundX-1]->SetLineColor(kYellow);
         lineplus[numli+nfoundX-1]->SetLineWidth(1);
         map_withPixels->cd();
         lineplus[numli+nfoundX-1]->Draw();
     }
  }

  if (nfoundX==7 && sqrt(64-count_missed)-floor(sqrt(64-count_missed))==0){
     for (Int_t numli=0; numli<nfoundX-1; numli++){
         Double_t CX[9];
         Double_t CY[9];
         Double_t LX[9];
         Double_t LY[9];

         CX[0]=((((Px[numli]+Px[numli+1])/2)-((Px[numli+nfoundX]+Px[numli+nfoundX+1])/2))*(1+((Py[numli+nfoundX]+Py[numli+nfoundX+1])/2))/(((Py[numli+nfoundX]+Py[numli+nfoundX+1])/2)-(Py[numli]+Py[numli+1])/2))+((Px[numli+nfoundX]+Px[numli+nfoundX+1])/2);
         CY[0]=-1;
         LX[0]=-1;
         LY[0]=((Py[nfoundX*numli]+Py[(nfoundX*numli)+nfoundX])/2)+((((Py[(nfoundX*numli)+1]+Py[(nfoundX*numli)+nfoundX+1])/2)-((Py[nfoundX*numli]+Py[(nfoundX*numli)+nfoundX])/2))*(1+((Px[nfoundX*numli]+Px[(nfoundX*numli)+nfoundX])/2))/(((Px[nfoundX*numli]+Px[(nfoundX*numli)+nfoundX])/2)-((Px[(nfoundX*numli)+1]+Px[(nfoundX*numli)+nfoundX+1])/2)));

         for (Int_t m=1; m<8; m++){
             CX[m]=(Px[numli+((m-1)*nfoundX)]+Px[numli+((m-1)*nfoundX)+1])/2;
             CY[m]=(Py[numli+((m-1)*nfoundX)]+Py[numli+((m-1)*nfoundX)+1])/2;
             LX[m]=(Px[(nfoundX*numli)+m-1]+Px[(nfoundX*numli)+nfoundX+m-1])/2;
             LY[m]=(Py[(nfoundX*numli)+m-1]+Py[(nfoundX*numli)+nfoundX+m-1])/2;
         }

         CX[8]=((((Px[numli+(5*nfoundX)]+Px[numli+(5*nfoundX)+1])/2)-((Px[numli+(6*nfoundX)]+Px[numli+(6*nfoundX)+1])/2))*(1-((Py[numli+(5*nfoundX)]+Py[numli+(5*nfoundX)+1])/2))/(((Py[numli+(5*nfoundX)]+Py[numli+(5*nfoundX)+1])/2)-(Py[numli+(6*nfoundX)]+Py[numli+(6*nfoundX)+1])/2))+((Px[numli+(5*nfoundX)]+Px[numli+(5*nfoundX)+1])/2);
         CY[8]=1;
         LX[8]=1;
         LY[8]=((Py[(nfoundX*numli)+5]+Py[(nfoundX*numli)+nfoundX+5])/2)+((((Py[(nfoundX*numli)+5]+Py[(nfoundX*numli)+nfoundX+5])/2)-((Py[(nfoundX*numli)+6]+Py[(nfoundX*numli)+nfoundX+6])/2))*(1-((Px[(nfoundX*numli)+5]+Px[(nfoundX*numli)+nfoundX+5])/2))/(((Px[(nfoundX*numli)+5]+Px[(nfoundX*numli)+nfoundX+5])/2)-((Px[(nfoundX*numli)+6]+Px[(nfoundX*numli)+nfoundX+6])/2)));

         lineplus[numli] = new TPolyLine(nfoundX+2,CX,CY);
         lineplus[numli]->SetLineColor(kYellow);
         lineplus[numli]->SetLineWidth(2);
         map_withPixels->cd();
         lineplus[numli]->Draw();
         lineplus[numli+nfoundX-1] = new TPolyLine(nfoundX+2,LX,LY);
         lineplus[numli+nfoundX-1]->SetLineColor(kYellow);
         lineplus[numli+nfoundX-1]->SetLineWidth(2);
         map_withPixels->cd();
         lineplus[numli+nfoundX-1]->Draw();
     }
  }



  if (nfoundX==6 && sqrt(64-count_missed)-floor(sqrt(64-count_missed))==0){
     for (Int_t numli=0; numli<nfoundX-1; numli++){
         Double_t CX[8];
         Double_t CY[8];
         Double_t LX[8];
         Double_t LY[8];

         CX[0]=((((Px[numli]+Px[numli+1])/2)-((Px[numli+nfoundX]+Px[numli+nfoundX+1])/2))*(1+((Py[numli+nfoundX]+Py[numli+nfoundX+1])/2))/(((Py[numli+nfoundX]+Py[numli+nfoundX+1])/2)-(Py[numli]+Py[numli+1])/2))+((Px[numli+nfoundX]+Px[numli+nfoundX+1])/2);
         CY[0]=-1;
         LX[0]=-1;
         LY[0]=((Py[nfoundX*numli]+Py[(nfoundX*numli)+nfoundX])/2)+((((Py[(nfoundX*numli)+1]+Py[(nfoundX*numli)+nfoundX+1])/2)-((Py[nfoundX*numli]+Py[(nfoundX*numli)+nfoundX])/2))*(1+((Px[nfoundX*numli]+Px[(nfoundX*numli)+nfoundX])/2))/(((Px[nfoundX*numli]+Px[(nfoundX*numli)+nfoundX])/2)-((Px[(nfoundX*numli)+1]+Px[(nfoundX*numli)+nfoundX+1])/2)));

         for (Int_t m=1; m<7; m++){
             CX[m]=(Px[numli+((m-1)*nfoundX)]+Px[numli+((m-1)*nfoundX)+1])/2;
             CY[m]=(Py[numli+((m-1)*nfoundX)]+Py[numli+((m-1)*nfoundX)+1])/2;
             LX[m]=(Px[(nfoundX*numli)+m-1]+Px[(nfoundX*numli)+nfoundX+m-1])/2;
             LY[m]=(Py[(nfoundX*numli)+m-1]+Py[(nfoundX*numli)+nfoundX+m-1])/2;
         }

         CX[7]=((((Px[numli+(4*nfoundX)]+Px[numli+(4*nfoundX)+1])/2)-((Px[numli+(5*nfoundX)]+Px[numli+(5*nfoundX)+1])/2))*(1-((Py[numli+(4*nfoundX)]+Py[numli+(4*nfoundX)+1])/2))/(((Py[numli+(4*nfoundX)]+Py[numli+(5*nfoundX)+1])/2)-(Py[numli+(5*nfoundX)]+Py[numli+(5*nfoundX)+1])/2))+((Px[numli+(4*nfoundX)]+Px[numli+(4*nfoundX)+1])/2);
         CY[7]=1;
         LX[7]=1;
         LY[7]=((Py[(nfoundX*numli)+4]+Py[(nfoundX*numli)+nfoundX+4])/2)+((((Py[(nfoundX*numli)+4]+Py[(nfoundX*numli)+nfoundX+4])/2)-((Py[(nfoundX*numli)+5]+Py[(nfoundX*numli)+nfoundX+5])/2))*(1-((Px[(nfoundX*numli)+4]+Px[(nfoundX*numli)+nfoundX+4])/2))/(((Px[(nfoundX*numli)+4]+Px[(nfoundX*numli)+nfoundX+4])/2)-((Px[(nfoundX*numli)+5]+Px[(nfoundX*numli)+nfoundX+5])/2)));

         lineplus[numli] = new TPolyLine(nfoundX+2,CX,CY);
         lineplus[numli]->SetLineColor(kYellow);
         lineplus[numli]->SetLineWidth(2);
         map_withPixels->cd();
         lineplus[numli]->Draw();
         lineplus[numli+nfoundX-1] = new TPolyLine(nfoundX+2,LX,LY);
         lineplus[numli+nfoundX-1]->SetLineColor(kYellow);
         lineplus[numli+nfoundX-1]->SetLineWidth(2);
         map_withPixels->cd();
         lineplus[numli+nfoundX-1]->Draw();
     }
  }
}
//---------------------------------------000000000000000000----------------------------------------//
