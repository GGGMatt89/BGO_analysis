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
#include <dirent.h>//class used to manipulate and traverse directories
#include <stdio.h>
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


  //Creation file for writing raw data
  TString path_raw = "/media/oreste/DATA/BGO/";

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
          for(N_PM = 0; N_PM<4; N_PM++){
          std::ofstream file_raw_data(path_raw + Form("Raw_Data_RUN_%i_%i_%i.txt", runN, fileAn, N_event_valid/4+1));
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
                     if (sample_rec.sample>4096) { // Verification of saturation
                     invalid_sple_val = invalid_sple_val+1;
                     }
                   } // end loop for sple

                   if (hit_block_rec.N_PM_rec>=0){
                   invalid_sample[hit_block_rec.N_PM_rec%4]=invalid_sple_val;
                   std::cout<<"invalid_sample :"<<'\t';
                   for (int test=0; test<4; test++){std::cout<< invalid_sample[test]<<'\t';
                   }
                   std::cout<<""<<std::endl;
                   invalid_sple_val = 0.;
                   }
                   else{std::cout<<"None PM hit"<<std::endl;}

               PM_signal_valid_counter=0;
               for (int verification=0; verification<4; verification++){
                 if (invalid_sample[verification]<=1024) {
                 PM_signal_valid_counter=PM_signal_valid_counter+1;
                 }
               }
               sum_sample=0.;
               for(int sum_col=0; sum_col<4; sum_col++) { // Sum of sample contents to verify if matrix is null
                for(int sum_line=0; sum_line<1024; sum_line++) {
              sum_sample=sum_sample+hist_sample[sum_col][sum_line];
                }
              }

            }// end loop fo w
            bzero(hit_structure_mode1, 7);
          }//end mode num 22
            else {std::cout<<"PROBLEM MODE NUMBER"<<std::endl;}
          } // end module num
        } // end if data_beg_id
        else {printf("PROBLEM IN READING DATA HEADER : beg_id = x%x\n", data_struct.data_beg_id);}
        if(N_event_valid%4==0) {
                        for (int sample_write=0 ; sample_write<1024 ; sample_write++) {
                        file_raw_data<<hist_sample[0][sample_write]<<'\t'<<hist_sample[1][sample_write]<<'\t'<<hist_sample[2][sample_write]<<'\t'<<hist_sample[3][sample_write]<<std::endl;
                        }

                        for (int fill_col=0;fill_col<4;fill_col++){ // Table content removal
                          for (int fill_line=0; fill_line<1024; fill_line++){
                          hist_sample[fill_col][fill_line]=0.;
                          }
                        }
                      } //end if sum_sample
                      else {continue;}
                      for (int clear_invalid_sample=0; clear_invalid_sample<4; clear_invalid_sample++){
                      invalid_sample[clear_invalid_sample]=0;
                      }
       } // end for N_PM
      } // end Nevnt
      } //end if fp
    TFile outroot(Form("/media/oreste/DATA/BGO/Output_Analysis/run%d_%i.root", runN, fileAn), "RECREATE");
    } // end if fileAn

} //end if cnt

  return 0;
}

