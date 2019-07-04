//c++ classes
#include <stdlib.h>
#include <strings.h>
#include <stdio.h>
#include <arpa/inet.h>
#include <unistd.h>
#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <stdint.h>
#include <cstring>
#include <cstdio>
#include <cstddef>
#include <dirent.h>//class used to manipulate and traverse directories
#include <stdio.h>
#include <math.h>
#include <mutex>
#include <vector>

//ROOT classes
#include "TLatex.h"
#include "TROOT.h"
#include "TObject.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TF1.h"
#include "TFile.h"
#include "TApplication.h"
#include "TRint.h"
#include "TAxis.h"
#include "TAttLine.h"
#include "TTimer.h"
#include "TStopwatch.h"
#include "TSystem.h"
#include "TStyle.h"
#include "TGWindow.h"
#include "TGClient.h"
#include "TPaveText.h"
#include "TMath.h"

//custom classes
#include "functions_new.h"

extern int alphasort();

void decimal_to_binary (int* arr,int &x, int &y, int decimal)
{
  int i=0,r;
  while(decimal!=0)
  {
   r = decimal%2;
   arr[i] = r;
   decimal /= 2;
   i++;
  }
  /*for (int l=0;l<=7;l++) {
  arr[l]=arr[l]*pow(2,16+l);
  x=x+arr[l];
  }*/

  if (arr[14]==1) {y = 1;}
  else {y = 0;} 
    
}

int file_select(const struct dirent *entry)//selection of files to be included in the return of the folder scan
{//ignore all files with the following features
  if((strcmp(entry->d_name, ".directory") == 0) ||(strcmp(entry->d_name, ".DS_Store") == 0) || (strcmp(entry->d_name, ".") == 0) || (strcmp(entry->d_name, "..") == 0) || strstr(entry->d_name,"tmp")!=NULL||strstr(entry->d_name,"pulse")!=NULL||strstr(entry->d_name,"T")!=NULL){
    return 0;
  }  else {
    return 1;
  }
}

int main (int argc,char ** argv)
{
  int runN = 0;//to be inserted by the user
  std::cout<<"Type the run number : (return to confirm)"<<std::endl;
  std::cin>>runN;
  char path[100] = "/media/oreste/DATA/BGO/Aquisitions_BGO_ASM/";//Cosmic_DATA/";//"/home/daq/gamhadron/daqGh/data/";//folder where the data are saved
  char filemain[100];//main part of the file name - only numbering stuff after that
  int xx = sprintf(filemain, "%d-", runN);
  char format[30] =".dat";//data format -> the acquisition software creates .dat files

  //-----------------------------------------------------------------------------------------------------------------------------//
  //---------------------------------------------DON'T MODIFY SINCE NOW ON (if not needed :-) )----------------------------------//
  //-----------------------------------------------------------------------------------------------------------------------------//

  //ROOT STYLING
  gROOT->SetStyle("Plain");
  gROOT->ForceStyle();
  gStyle -> SetStatW(0.28);
  gStyle -> SetStatH(0.13);
  gStyle -> SetStatColor(0);
  gStyle -> SetStatX(0.87);
  gStyle -> SetStatY(0.85);
  gStyle -> SetStatFont(0);
  gStyle -> SetOptStat(111);
  gStyle -> SetPalette(1);
  gStyle -> SetOptTitle(1);

//-------------------------------//
  std::cout<<"STARTING analysis for run "<<runN<<" ... "<<std::endl;
  //Creating data structure to read the data from file
  FILEHEADER file_beg;//file header - it contains 0xF0F0 (16 bit), the run number (32 bit), the total number of events recorder in the file (32 bit)
  FILEEND file_end;//file footer - it contains 0xF1F1 (16 bit), the total number of events recorder in the file (32 bit), the total number of octets recorded in the file (32 bit)
  EVENTHEADER event_id;//event header (one per event) - it contains 0xABCD (16 bit), the event number (32 bit), the trigger number (24 bit), the number of hits recorded in this trigger (16 bit)
  DATAMAIN data_struct;//main structure of data (one per event) - it contains 0x00EB (16 bit), the front end number (8 bit), the trigger number (24 bit), the mode number (8 bit), the number of involved detector modules (8 bit - example number of touched blocks in the absorber)
  hit_block_receive hit_block_rec;//data for each hit mode 1 - test (one per involved detector module) - it contains the number of the hit module (8 bit), the recorded time for the interaction (32 bit) and the recorded charge (16 bit)
  sample_receive sample_rec;//data contained in samples
  //Create ID to check what I read
  uint16_t beginning_file_id = 0xF0F0;
  uint16_t end_file_id = 0xF1F1;
  uint16_t beginning_event = 0xABCD;
  uint8_t beginning_data_main = 0xEB;
  int N_ev_toRead = 0;//Set the number of events you want to read from each file for monitoring purpose - depending on the speed and accuracy you need
  int length;
  int file_to_an = 0;
  int N_event=0;
  int N_event_valid = 0;
  int N_event_valid_coinc = 0;
  double N_event_valid_X =0;
  double N_event_valid_Y =0;

  int count = 0;
  int count1 = 0;
  int files_read = 0;
  int events_read = 0;
  int wait_loop = 0;//counter for waiting loops
  int time_to_wait = 200;//how many waiting loops you accept before stopping the program

  //Initializing char buffers for reading from file
  char init_file[10];//file header
  char end_file[8];//file end - for debugging
  uint8_t event_header[11];//event header
  uint8_t data_main_structure[7];//data main structure
  uint8_t hit_structure_mode0[17];//mode 20
  uint8_t hit_structure_mode1[7];//mode 22
  uint8_t sample_structure[2];

  //Starting ROOT in the c++ program
  TApplication theApp("App", &argc, argv);
  if (gROOT->IsBatch()) {
    fprintf(stderr, "%s: cannot run in batch mode\n", argv[0]);
    return 1;
  }
//--------------------------------------------------------------------------------//
//---------------------------***PLOTS SETTINGS***---------------------------------//
//--------------------------------------------------------------------------------//



  //Variables for files reading and data extraction
  int cnt = 0;
  double sum_sample=0.;
  int threshold_high=4096; // modify if needed
  int threshold_low=0; // modify if needed
  double hist_sample[4][1024];
  for (int fill_col=0;fill_col<4;fill_col++){
    for (int fill_line=0; fill_line<1024; fill_line++){
    hist_sample[fill_col][fill_line]=0.;
    }
  }
  int invalid_sple_val =0;
  int invalid_sample[4];
  for (int fill=0; fill<4; fill++) {invalid_sample[fill]=0;}
  int PM_signal_valid_counter=0;


  //Variables for signal treatment (crossing high threshold, baseline, events removal)
  double baseline[4], mean_1[4], mean_2[4], integral[4], sum_1[4], sum_2[4]; // mean_1 and mean_2 are used to verify if the variation of values in the baseline calculation is low. 
  for (int fill_bsl=0;fill_bsl<4;fill_bsl++){
  baseline[fill_bsl]=0.;
  mean_1[fill_bsl]=0.;
  mean_2[fill_bsl]=0.;
  integral[fill_bsl]=0.;
  }
  int PM_1_binary =0, PM_2_binary = 0, PM_3_binary = 0, PM_4_binary = 0;
  int FourPM_event1=0, FourPM_event2=0, FourPM_event3=0, FourPM_event4=0, FourPM_event_tot=0; 
  int Bit_14[16]; // used to count the number of events where the four PM are hit

  //Creation file for writing raw data
  TString path_raw = "/media/oreste/DATA/BGO/Data_Extraction/";
  int blockNumber;

  //Creating string to check that the selected files refer to the present run - all other data files are ignored
  std::stringstream  s;
  s<<runN;
  s<<"-";
  //Creating structures and variables to search in the data folder for files
  struct dirent **files;
  std::stringstream  ss;
  struct dirent **files2;
  std::stringstream  sss;
  std::fstream fp;
  bool cond = true;
  bool flag = 0;

  count = scandir(path, &files, file_select, alphasort);//how many files in the folder?
  if(count>0){
  std::cout<<"Number of data file found "<<count<<std::endl;
  std::cout<<"Recovered file list : "<<std::endl;
    for(int x = 0; x<count; x++){
      if(strstr(files[x]->d_name, s.str().c_str())){
	      printf("-> %s", files[x]->d_name);
	      std::cout<<std::endl;
	      cnt++;
	    }else{
	      continue;
	    }
    }
  }


if(cnt>0){
    std::ofstream integral_calculated_data(path_raw + Form("Complete_Data_RUN_%i.txt", runN));
    std::ofstream charge_calculation_data(path_raw + Form("Bit_charge_calculation_content_RUN_%i.txt", runN));
    std::ofstream FourPM_event(path_raw + Form("4_PM_event_content_RUN_%i.txt", runN));
    std::cout<<"We collected "<<cnt<<" files for run "<<runN<<std::endl;
    std::cout<<"How many files you want to analyse ? (Type here and return (-1 for all))"<<std::endl;
    std::cin>>file_to_an;
    if(file_to_an==-1){file_to_an = cnt - 2;}
    for(int fileAn = 0; fileAn < file_to_an; fileAn++){
      if(fp.is_open()){fp.close();}
      sss.str("");
      sss<<path;
      sss<<filemain;
      sss<<fileAn;
      sss<<format;
      fp.open(sss.str().c_str(), std::fstream::in | std::fstream::binary);
      if(fp){
        std::cout<<"FILE open"<<std::endl;
        std::cout<<"Analysing file "<<sss.str().c_str()<<" ..... "<<std::endl;
        fp.seekg(0, std::ios::end);
        length = fp.tellg();
        std::cout<<"The file size is "<<length<<std::endl;
        fp.seekg(0, std::ios::beg);//put the read pointer at the beginning of the file+
        //Reading the file header
        fp.read((char *)init_file, sizeof(init_file));//reading the file header
        unpack_fileHead((unsigned char *)init_file, file_beg);
        if(file_beg.file_beg_id==beginning_file_id){
          std::cout<<"Beginning of new file read well "<<std::endl;
          std::cout<<"This is run "<<file_beg.run_number<<" and this file contains "<<file_beg.Ntot_events<<std::endl;
        }else{
          std::cout<<"PROBLEM IN READING FILE HEADER -> CHECK THE FILE STRUCTURE! Moving to the next file"<<std::endl;
          continue;
        }

        for(int Nevnts = 0; Nevnts<file_beg.Ntot_events; Nevnts++){
          fp.read((char *)event_header, sizeof(event_header));//reading the event header
          printf("event_header = x%02x %02x %02x %02x %02x %02x %02x %02x %02x %02x %02x\n",
        		  event_header[0],event_header[1],event_header[2],event_header[3],
        		  event_header[4],event_header[5],event_header[6],event_header[7],
				  event_header[8],event_header[9],event_header[10]);
          unpack_eventHead((unsigned char *)event_header, event_id);
          if(event_id.event_beg_id != beginning_event){
            std::cout<<"PROBLEM IN READING DATA -> CHECK THE DATA STRUCTURE! Moving to the next file"<<std::endl;
            printf("event_beg_id=x%x\t event_number=x%x\t trigger_number=x%x\t hit_in_trig=x%x\n",
            	        		  event_id.event_beg_id,
            					  event_id.event_number,
            					  event_id.trigger_number,
            					  event_id.hit_in_trig);
              std::cout<<beginning_event<<std::endl;
	          std::cout<<event_id.event_beg_id<<std::endl;
	          std::cout<<event_id.event_number<<std::endl;
	          std::cout<<event_id.trigger_number<<std::endl;
	          std::cout<<event_id.hit_in_trig<<std::endl;
	          break;
          }
          int N_PM = 0;
          for(N_PM = 0; N_PM<event_id.hit_in_trig; N_PM++){
          fp.read((char *)data_main_structure, sizeof(data_main_structure));//reading the data main part

          unpack_dataMain((unsigned char *)data_main_structure, data_struct);
        if(data_struct.data_beg_id==beginning_data_main){
          std::cout<<"Data header read well "<<std::endl;
          std::cout<<"N module debug :"<< data_struct.modules_num<<std::endl;
          if(data_struct.modules_num>0){
            N_event_valid++;
            std::cout<<N_event_valid<<std::endl;
            std::cout<<data_struct.mode_num<<std::endl;
            std::cout<<data_struct.modules_num<<std::endl;
            if(data_struct.mode_num==20){
		 std::cout<<"OPTIMIZED MODE NOT PROGRAMMED YET..."<<std::endl;
            }
            if(data_struct.mode_num==22){ //mode 22
		 for(int w = 0; w<data_struct.modules_num; w++){//reading the data from each involved block
	         fp.read((char *)hit_structure_mode1, sizeof(hit_structure_mode1));
	         unpack_data((unsigned char *)hit_structure_mode1, hit_block_rec);
                 blockNumber=trunc((hit_block_rec.N_PM_rec/4)+1);
                 //std::ofstream file_raw_data(path_raw + Form("Raw_Data_RUN_%i_block_%i_trigg_%i.txt", runN, blockNumber, event_id.trigger_number));  //--> used for debugg
                 std::cout<<hit_block_rec.N_PM_rec<<'\t'<<hit_block_rec.hit_time_rec<<'\t'<<hit_block_rec.Nb_samples_rec<<std::endl;
                   for (int sple = 0; sple<hit_block_rec.Nb_samples_rec; sple++){
                    fp.read((char *)sample_structure, sizeof(sample_structure));
                    unpack_sample((unsigned char *)sample_structure, sample_rec);
                    //std::cout<<sample_rec.sample<<std::endl;
//--------------------------------------------------------------------------------------------//
//----------------------***Handling raw data and Writing valid events***----------------------//
//--------------------------------------------------------------------------------------------//
                    hist_sample[hit_block_rec.N_PM_rec%4][sple] = sample_rec.sample;
                    //sample-> Fill(sample_rec.sample);
                     if (sample_rec.sample>=4095) { // Verification of saturation
                     invalid_sple_val = invalid_sple_val+1;
                     }
                   

                     if (hit_block_rec.N_PM_rec>=0){
                     invalid_sample[hit_block_rec.N_PM_rec%4]=invalid_sple_val;
                   } 
                   else{std::cout<<"None PM hit"<<std::endl;}

                   } // end loop for sple
                     std::cout<<"invalid_sample :"<<'\t';
                     for (int invalid=0; invalid<4; invalid++){std::cout<< invalid_sample[invalid]<<'\t';
                     }
                   std::cout<<""<<std::endl;
                   invalid_sple_val = 0.;
                   

               PM_signal_valid_counter=0;
               for (int verification=0; verification<4; verification++){
                 if (invalid_sample[verification]<=40) {
                 PM_signal_valid_counter=PM_signal_valid_counter+1;
                 }
               }
               sum_sample=0.;
               for(int sum_col=0; sum_col<4; sum_col++) { // Sum of sample contents to verify if matrix is null
                for(int sum_line=0; sum_line<1024; sum_line++) {
              sum_sample=sum_sample+hist_sample[sum_col][sum_line];
                }
              }
                        if(N_event_valid%4==0 && PM_signal_valid_counter==4) {
                       /* for (int sample_write=0 ; sample_write<1024 ; sample_write++) {
                        file_raw_data<<hist_sample[0][sample_write]<<'\t'<<hist_sample[1][sample_write]<<'\t'<<hist_sample[2][sample_write]<<'\t'<<hist_sample[3][sample_write]<<std::endl;
                        }*/

//--------------------------------------------------------------------------------------------//
//---------------------------------***Integral calculation***---------------------------------//
//--------------------------------------------------------------------------------------------//
                            //baseline calculation
			    for(int baseline_calc=0; baseline_calc<4; baseline_calc++){
			     for (int i=2;i<18;i++){ // the two first samples are dedicated to charge calculation and baseline calculation in data format
			     sum_1[baseline_calc]= sum_1[baseline_calc] + hist_sample[baseline_calc][i];
			     }
			     mean_1[baseline_calc]=sum_1[baseline_calc]/16;
			     for (int ii=18;ii<33;ii++){
			     sum_2[baseline_calc]= sum_2[baseline_calc] + hist_sample[baseline_calc][ii];
			     }
			     mean_2[baseline_calc]=sum_2[baseline_calc]/16;
			     if (mean_1[baseline_calc]-mean_2[baseline_calc]<50 || mean_1[baseline_calc]-mean_2[baseline_calc]<-50){
			     baseline[baseline_calc]=(mean_1[baseline_calc]+mean_2[baseline_calc])/2;}
			     else {baseline[baseline_calc]=mean_1[baseline_calc];}
			    }
                           std::cout<<"Baseline : "<<std::endl;
                           for (int test=0; test<4; test++){std::cout<<baseline[test]<<'\t';}
                           std::cout<<""<<std::endl;
                           //integral calculation
                            for(int integ_calc=0; integ_calc<4; integ_calc++){
			     for (int l=2;l<1023;l++){
                               if (hist_sample[integ_calc][l] - baseline[integ_calc] < 0){
                               integral[integ_calc]=integral[integ_calc]+0;
                               }
                               else {integral[integ_calc]=integral[integ_calc]+hist_sample[integ_calc][l]- baseline[integ_calc];}
                             }
                           
                            }
                            for(int charge_size_adjustement=0; charge_size_adjustement<4; charge_size_adjustement++){
                             integral[charge_size_adjustement] = integral[charge_size_adjustement]/379;
                            }
                             std::cout<<"Integrale : "<<std::endl;
                             for (int test1=0; test1<4; test1++){std::cout<<integral[test1]<<'\t';}
                             std::cout<<""<<std::endl;
                           //writting in output file
                           for(int fill_complete_file=0; fill_complete_file<data_struct.modules_num; fill_complete_file++){
			    integral_calculated_data<<integral[0]<<'\t'<<integral[1]<<'\t'<<integral[2]<<'\t'<<integral[3]<<std::endl;

                           //Charge calculation verification 
                            PM_1_binary =0, PM_2_binary = 0, PM_3_binary = 0, PM_4_binary = 0, FourPM_event1=0, FourPM_event2=0, FourPM_event2=0, FourPM_event4=0, FourPM_event_tot=0;
                            for (int i2=0; i2<16; i2++){Bit_14[i2]=0;}
                            decimal_to_binary(Bit_14, PM_1_binary, FourPM_event1, hist_sample[0][0]);
                            std::cout<<"Bit_14 : ";
                            for (int i1=0; i1<16; i1++){
                            std::cout<<Bit_14[i1]<<'\t';}
                            std::cout<<"Y : "<<FourPM_event1;
                            std::cout<<""<<std::endl;
                            for (int i2=0; i2<16; i2++){Bit_14[i2]=0;}
                            decimal_to_binary(Bit_14, PM_2_binary, FourPM_event2, hist_sample[1][0]);
                            std::cout<<"Bit_14 : ";
                            for (int i1=0; i1<16; i1++){
                            std::cout<<Bit_14[i1]<<'\t';}
                            std::cout<<"Y : "<<FourPM_event2;
                            std::cout<<""<<std::endl;
                            for (int i2=0; i2<16; i2++){Bit_14[i2]=0;}
                            decimal_to_binary(Bit_14, PM_3_binary, FourPM_event3, hist_sample[2][0]);
                            std::cout<<"Bit_14 : ";
                            for (int i1=0; i1<16; i1++){
                            std::cout<<Bit_14[i1]<<'\t';}
                            std::cout<<"Y : "<<FourPM_event3;
                            std::cout<<""<<std::endl;
                            for (int i2=0; i2<16; i2++){Bit_14[i2]=0;}
                            decimal_to_binary(Bit_14, PM_4_binary, FourPM_event4, hist_sample[3][0]);
                            std::cout<<"Bit_14 : ";
                            for (int i1=0; i1<16; i1++){
                            std::cout<<Bit_14[i1]<<'\t';}
                            std::cout<<"Y : "<<FourPM_event4;
                            std::cout<<""<<std::endl;

                            PM_1_binary= PM_1_binary + hist_sample[0][1];
                            PM_2_binary= PM_2_binary + hist_sample[1][1];
                            PM_3_binary= PM_3_binary + hist_sample[2][1];
                            PM_4_binary= PM_4_binary + hist_sample[3][1];

                            for (int i2=0; i2<16; i2++){Bit_14[i2]=0;}

                            FourPM_event_tot = FourPM_event1 + FourPM_event2 + FourPM_event3 + FourPM_event4;
                            if (FourPM_event_tot == 4) { FourPM_event<<integral[0]<<'\t'<<integral[1]<<'\t'<<integral[2]<<'\t'<<integral[3]<<std::endl;}
                            

			    charge_calculation_data<<"Charge_PM1 : "<<PM_1_binary<<'\t'<<"Charge_PM2 : "<<PM_2_binary<<'\t'<<"Charge_PM3 : "<<PM_3_binary<<'\t'<<"Charge_PM4 : "<<PM_4_binary<<std::endl;
//"PM1_Bit1 : "<<hist_sample[0][0]<<'\t'<<"PM2_Bit1 : "<<hist_sample[1][0]<<'\t'<<"PM3_Bit1 : "<<hist_sample[2][0]<<'\t'<<"PM4_Bit1 : "<<hist_sample[3][0]<<'\n'<<"PM1_Bit2 : "<<hist_sample[0][1]<<'\t'<<"PM2_Bit2 : "<<hist_sample[1][1]<<'\t'<<"PM3_Bit2 : "<<hist_sample[2][1]<<'\t'<<"PM4_Bit2 : "<<hist_sample[3][1]<<std::endl;
                           }

                        for (int fill_col=0;fill_col<4;fill_col++){ // Table content removal
                          integral[fill_col]=0.;
                          sum_1[fill_col]=0.;
                          sum_2[fill_col]=0.;
                          mean_1[fill_col]=0.;
                          mean_2[fill_col]=0.;
                          baseline[fill_col]=0.;
                          for (int fill_line=0; fill_line<1024; fill_line++){
                          hist_sample[fill_col][fill_line]=0.;
                          }
                        }
                           
                       } //end if N_event_valid




            }// end loop fo w
            bzero(hit_structure_mode1, 7);
          }//end mode num 22
            else {std::cout<<"PROBLEM MODE NUMBER"<<std::endl;}
          } // end module num
        } // end if data_beg_id
        else {printf("PROBLEM IN READING DATA HEADER : beg_id = x%x\n", data_struct.data_beg_id);}        
       } // end for N_PM
      for (int clear_invalid_sample=0; clear_invalid_sample<4; clear_invalid_sample++){
         invalid_sample[clear_invalid_sample]=0;
         }
      } // end Nevnt
      } //end if fp
    } // end if fileAn

} //end if cnt

  return 0;
}

