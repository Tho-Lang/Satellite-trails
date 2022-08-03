// STL
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <cstdlib>
#include <sys/stat.h>

// ROOT


// HESS
#include <sash/HESSArray.hh>
#include <sash/Telescope.hh>
#include <sash/TelescopeConfig.hh>
#include <sash/Pixel.hh>
#include <sash/ListIterator.hh>
#include <sash/List.hh>
#include <sash/PixelConfig.hh>
#include <sash/RunHeader.hh>
#include <sash/DataSet.hh>
#include <sash/Pointer.hh>
#include <sash/PointerSet.hh>
#include <sash/PointerSetIterator.hh>

#include <calibration/TelescopeNSB.hh>
#include <calibration/PixelNSB.hh>

#include <utilities/SimpleStats.hh>

#include <sashfile/HandlerC.hh>
#include <stash/Coordinate.hh>

stringstream ss;
stringstream ss_avg_pix;
stringstream ss_delta_pix;
stringstream ss_high_pix;
stringstream ss_roll_avg;


string home_path = "/home/hpc/caph/mppi117h/";



void PrintAverageNSB(int firstrun) {

  if (firstrun % 200 !=0) {
    std::cout << "Choose number divisible by 200" << std::endl;
    return;
  }
  
  stringstream run_numbers;
  run_numbers << "run"<<firstrun<<"-"<<firstrun+199<<"/";
  
  
  stringstream make_eval_folder;
  string eval_folder = "eval_data";
  make_eval_folder << "mkdir -p "<< eval_folder<< "/";
  std::cout <<  make_eval_folder.str() << std::endl;
  const int dir_make_eval_folder = system( make_eval_folder.str().c_str());
  if (dir_make_eval_folder == -1) {
    printf("Error in creating evaluation directory!");
    exit(1);
  }
  
  int chdir_eval_folder = chdir(eval_folder.c_str());
 
  
  //Loop on all runs from x to x+200 like in $HESSDATA
  for (int run_number = 0; run_number<200; ++run_number) {
  int nrun = firstrun +run_number;

  SashFile::HandlerC han("");
  
  
  if (!han.ConnectFile("NSB",nrun)) {
    std::cout << "Can't find the NSB file for run " << nrun << std::endl;
    continue;
    return;
  }
  
  Sash::DataSet *ds_run = han.GetDataSet("run");
  if (!ds_run) {
    std::cout << "Can't find the DataSet [run]" << std::endl;
    return;
  }
  
  ofstream avg_tmp;
  ofstream avg_pix;
  ofstream delta_pix;
  ofstream high_pix;
  ofstream roll_avg;
  
  ds_run->GetEntry(0); // Load the first entry in memory. This contain the information about the run and the list of telescope that where involved.

  Sash::HESSArray *fhess = &Sash::HESSArray::GetHESSArray(); // Retrieve the singleton class that is the access point to every other storage class
  
  const Sash::RunHeader *runh = fhess->Get<Sash::RunHeader>();

  //Get zenith of run:
  const Stash::Coordinate &coord =  runh->GetObservationPosition().GetCoordinate(*fhess->GetGroundSystem());
  Double_t alt = coord.GetBeta().GetDegrees();
  Double_t az  = coord.GetLambda().GetDegrees();
  

  

  //Create directories 


  stringstream nrun_avg_pix; 
  nrun_avg_pix << "mkdir -p avg_per_pix/"<< run_numbers.str()<< nrun;
  std::cout <<  nrun_avg_pix.str() << std::endl;
  stringstream nrun_avg_tmp;
  nrun_avg_tmp <<  "mkdir -p avg_per_tmp/"<< run_numbers.str()<< nrun;
  std::cout <<  nrun_avg_tmp.str() << std::endl;
  stringstream nrun_delta_pix;
  nrun_delta_pix << "mkdir -p delta_pixel/"<< run_numbers.str() << nrun;
  std::cout <<  nrun_delta_pix.str() << std::endl;
  stringstream nrun_high_pix;
  nrun_high_pix << "mkdir -p high_pixel/"<< run_numbers.str() << nrun;
  std::cout <<  nrun_high_pix.str() << std::endl;
  stringstream nrun_roll_avg;
  nrun_roll_avg << "mkdir -p roll_avg/"<< run_numbers.str() << nrun;
  std::cout <<  nrun_roll_avg.str() << std::endl;
    
  const int dir_avg_pix = system( nrun_avg_pix.str().c_str());
  const int dir_avg_tmp = system( nrun_avg_tmp.str().c_str());
  const int dir_delta_pix = system( nrun_delta_pix.str().c_str());
  const int dir_high_pix = system( nrun_high_pix.str().c_str());
  const int dir_roll_avg = system( nrun_roll_avg.str().c_str());
  
  if (dir_avg_pix == -1) {
    printf("Error in creating average pixel directory!");
    exit(1);
  }
  if (dir_avg_tmp == -1) {
    printf("Error in creating average timestamp directory!");
    exit(1);
  }
  if (dir_delta_pix == -1) {
    printf("Error in creating delta pixel directory!");
    exit(1);
  }
  if (dir_high_pix == -1) {
    printf("Error in creating high pixel directory!");
    exit(1);
  }
  if (dir_roll_avg == -1) {
    printf("Error in creating roll_avg directory!");
    exit(1);
  }

  
  // Now loop on the telescope in run
  const Sash::PointerSet<Sash::Telescope>& telsinrun = runh->GetTelsInRun();
  for( Sash::PointerSet<Sash::Telescope>::iterator tel = telsinrun.begin(); tel != telsinrun.end(); ++tel) {
    int telId = (*tel)->GetId();
    
    // At the moment, only the NSB derived from Pedestal width is valid for CT1-4 and HESS2_2048 camera
    // For FlashCam, this is the one derived on event that is valid
    std::string ds_suffix;
    std::string class_name;
    if ( (*tel)->GetConfig()->GetSetupID() == Sash::TelescopeConfig::FlashCam_1764 ) {
      ds_suffix = "EvtNSB";
      class_name = "NSBEvent";
    }
    else {
      ds_suffix = "PedNSB";
      class_name = "NSBPed";
    }
    
    // Retrieve the DataSet
    
    std::ostringstream oss_ds_name;
    oss_ds_name << "CT" << telId << "_" << ds_suffix;

    std::cout << "Running for CT" << telId << std::endl;

    Sash::DataSet *ds_nsb = han.GetDataSet(oss_ds_name.str());
    if (!ds_nsb) {
      std::cout << "Can't find the DataSet [" << oss_ds_name.str() << "]" << std::endl;
      continue;
    }

    //Open and name the files 

    ss <<           home_path << eval_folder << "/avg_per_tmp/" << run_numbers.str().c_str() << nrun << "/avg_per_tmp_"<< nrun <<"_CT_"<<telId<<".txt";
    ss_avg_pix <<   home_path << eval_folder  << "/avg_per_pix/" << run_numbers.str().c_str() << nrun << "/avg_per_pix_"<< nrun <<"_CT_"<<telId<<".txt";
    ss_delta_pix << home_path <<  eval_folder  << "/delta_pixel/" << run_numbers.str().c_str() << nrun << "/delta_pix_"<< nrun <<"_CT_"<<telId<<".txt";
    ss_high_pix << home_path << eval_folder  << "/high_pixel/" << run_numbers.str().c_str() << nrun << "/high_pix_"<< nrun <<"_CT_"<<telId<<".txt";
    ss_roll_avg << home_path << eval_folder  << "/roll_avg/" << run_numbers.str().c_str() << nrun << "/roll_avg_"<< nrun <<"_CT_"<<telId<<".txt";
    
    std::cout << ss.str()<<std::endl;
    std::cout << ss_avg_pix.str()<<std::endl;
    std::cout << ss_delta_pix.str()<<std::endl;
    std::cout << ss_high_pix.str()<<std::endl;
    std::cout << ss_roll_avg.str()<<std::endl;
    
    avg_tmp.open(ss.str().c_str());
    avg_pix.open(ss_avg_pix.str().c_str());
    delta_pix.open(ss_delta_pix.str().c_str());
    high_pix.open(ss_high_pix.str().c_str());
    roll_avg.open(ss_roll_avg.str().c_str());

    const int nevents = (int)ds_nsb->GetEntries();
    
    std::cout << nevents << std::endl;
    std::cout << ds_nsb->GetEntries() << std::endl;

  
   


    int pix_count = 0;
    for (Sash::Pointer<Sash::Pixel> pix = (*tel)->beginPixel(), end_pix = (*tel)->endPixel(); pix != end_pix; ++pix) {
      ++pix_count;
    }
    std::cout << pix_count << std::endl;
    const int size = pix_count;

    std::vector<float> pixel_avg_vector;
    float pixel_avg[size];    
    float pix_cur[size]; //last saved value of pixels 
    float pix_dif[size]; //difference of current to last value in pixel
    float threshold = 100;
    float deltathreshold = 30;
    
    if (telId == 5){
      threshold = 900;
      if (firstrun < 153001) {
	threshold = 20;
      }
    }
    
    
    
   
    
    
    //simple moving average: gives average of last N values. If later local avgerage is needed more memory will be needed to store last N*pix_count values
    std::vector<int> counter_vec;
    std::vector<float> sum_sma_vec;
    float sum_sma[size];
    
    int N = 40;
    if (telId == 5) {
      N = 120;
    }
    float x_val_minus_n_sma[N][size];
    std::vector< vector<float> > x_val_minus_n_sma_vec;
    
    for (int i = 0;i<N;++i) {
      std::vector<float> x_val_minus_i_sma_vec;
      for (int j = 0;j<size; ++j) {
	x_val_minus_i_sma_vec.push_back(0.);
	x_val_minus_n_sma[i][j] = 0.;
      }
      x_val_minus_n_sma_vec.push_back(x_val_minus_i_sma_vec);
    }
    
      // Now loop on all the dataset entry (one per pedestal evaluation)
    Utilities::SimpleStats uss_nsb;
    for ( int ie = 0 ; ie < ds_nsb->GetEntries(); ++ie ) { // For pedestal, the first entry is not necessarily a great estimate because the pedestal is not estimated on lot of events.... skip it ?
      ds_nsb->GetEntry( ie ); // Load the entry in memory
      
      
      Utilities::SimpleStats uss_nsb_evt; // new each event / time stamp

      // Get the Telescope NSB information
      const Calibration::TelescopeNSB *TelNSB = (*tel)->Get<Calibration::TelescopeNSB>(class_name.c_str()); // Get retrieve a const pointer (Handle will create the object if it does not exist in memory)
      if (!TelNSB) {
	std::cout << "Can't find the Calibration::TelescopeNSB class named " << class_name << " for CT" << telId << std::endl;
	continue;
      }
      TelNSB->LoadAllMembers(); // Not needed because TelescopeNSB contain simple data but nice to keep in mind that sometimes it is needed to force complex stored information like pointers) to be loaded in memory
      Sash::Time teltime_nsb = TelNSB->GetTimeStamp();

      //std::cout << "NSB recorded at: "  << teltime_nsb.GetUTC() << std::endl;
      avg_tmp<< teltime_nsb.GetTime()<< ";";

      pix_count = 0;
      float roll_avg_cur;
      float roll_avg_cur_vec;
      // Now loop on pixels
      for (Sash::Pointer<Sash::Pixel> pix = (*tel)->beginPixel(), end_pix = (*tel)->endPixel(); pix != end_pix; ++pix) {
	const Calibration::PixelNSB *PixNSB = TelNSB->GetMonitor(pix);
      
	if (!PixNSB) {
	  std::cout << "Can't find the PixelNSB class for CT" << telId << " and Pixel: "  << pix->GetConfig()->GetPixelID() << std::endl;
	  continue;
	}	
	if(ie==0) {
	  pixel_avg[pix_count]=0.;
	  pixel_avg_vector.push_back(0.);
	  counter_vec.push_back(0);
	  sum_sma_vec.push_back(0);
	  //std::cout <<"Value of sum_sma_vec in first timeslot: "<< sum_sma_vec[pix_count] <<std::endl;
	}
	
	if ( PixNSB->GetFlag() != 0 ) {	  
	  pix_cur[pix_count] = PixNSB->GetNSBValue();
	  if (pix_cur[pix_count]>0) {
	  uss_nsb.Fill( PixNSB->GetNSBValue() ); 
	  uss_nsb_evt.Fill( PixNSB->GetNSBValue() );
	  counter_vec[pix_count] = counter_vec[pix_count]+1;
	  
	  sum_sma[pix_count] = sum_sma[pix_count] + pix_cur[pix_count] - x_val_minus_n_sma[0][pix_count];
	  sum_sma_vec[pix_count] = sum_sma_vec[pix_count] + pix_cur[pix_count] - x_val_minus_n_sma_vec[0][pix_count];
	  //std::cout << sum_sma_vec[pix_count] <<std::endl; 
	  for (int i=0; i<N-2; ++i) {
	    x_val_minus_n_sma[i][pix_count] = x_val_minus_n_sma[i+1][pix_count];
	    x_val_minus_n_sma_vec[i][pix_count] = x_val_minus_n_sma_vec[i+1][pix_count];
	  }
	  x_val_minus_n_sma[N-1][pix_count] = pix_cur[pix_count]; //insert current values at end of array and use the first entry for x(n-N)
	  x_val_minus_n_sma_vec[N-1][pix_count] = pix_cur[pix_count];
	  if (counter_vec[pix_count]<N) {
	    roll_avg_cur = sum_sma[pix_count]/counter_vec[pix_count];
	    roll_avg_cur_vec = sum_sma_vec[pix_count]/counter_vec[pix_count];
	  }
	  else {
	    roll_avg_cur = sum_sma[pix_count]/N; 
	    roll_avg_cur_vec = sum_sma_vec[pix_count]/N; 
	  }
	  if(pix_cur[pix_count]> roll_avg_cur_vec + threshold) {
	    roll_avg << pix_count << "; "<< pix_cur[pix_count]<< ";"<< teltime_nsb.GetTime() << ";"
		     <<(*pix).GetConfig()->GetPos().GetX() << ";" <<(*pix).GetConfig()->GetPos().GetY()
		     << "; " << roll_avg_cur << "; " << sum_sma[pix_count]
		     <<std::endl;
	    std::cout << telId <<" counter: "<< counter_vec[pix_count] << ", "<< sum_sma[pix_count] << " vs. vector:  " << sum_sma_vec[pix_count] <<std::endl;
	  }
	  pixel_avg_vector[pix_count] = pixel_avg_vector[pix_count] + PixNSB->GetNSBValue();
	  pixel_avg[pix_count] = pixel_avg[pix_count] + PixNSB->GetNSBValue()/nevents;
	 	  
	  if (ie>0) {
	    pix_dif[pix_count] = PixNSB->GetNSBValue()- pix_cur[pix_count];
	    if (pix_dif[pix_count]>deltathreshold) {
	      //std::cout << "Pixel " << pix_count << " with dif> "<< pix_dif[pix_count] << " at"  << pix_cur[pix_count]<< " MHz at " <<teltime_nsb.GetUTC() << std::endl;
	      delta_pix << pix_count << "; "<< pix_dif[pix_count]<< ";"<< teltime_nsb.GetTime() << ";"
			<<(*pix).GetConfig()->GetPos().GetX() << ";" <<(*pix).GetConfig()->GetPos().GetY()
			<<std::endl;
	    }
	  }	  
	  if (pix_cur[pix_count]>threshold) {
	    //std::cout << "Pixel " << pix_count << " at " << pix_cur[pix_count]<< " MHz at " <<teltime_nsb.GetUTC() << std::endl;
	    high_pix << pix_count << ";" << pix_cur[pix_count]<< ";"<< teltime_nsb.GetTime() << ";"
		     <<(*pix).GetConfig()->GetPos().GetX() << ";" << (*pix).GetConfig()->GetPos().GetY()
		     << std::endl;
	  }
	  }
	}//remove one bracket, if condition pix_cur[pix_count]>0 is removed
	/*
	else {
	  std::cout << "CT" <<telId << ": pixel value of pixel " <<pix_count << " below zero at timestamp "<< teltime_nsb.GetTime() << std::endl;
	  }
*/

	/*
	if (ie == 0) {
	  
	  if (pix_count == 50) {
	    std::cout << "Pixel " << pix->GetConfig()->GetPixelID() << " located at "
		  << pix->GetConfig()->GetPos().GetX() << ", " << pix->GetConfig()->GetPos().GetY()
		  << std::endl;
	    std::cout <<"Neighbours:";
	    
	    const Sash::List<Sash::Pixel>& neighbours = pix->GetNeighbourList();
	    Utilities::SimpleStats uss_nsb_neighbours;
	    
	    for(Sash::List<Sash::Pixel>::const_iterator it = neighbours.begin(); it != neighbours.end(); ++it) {
	      int neighbourID = (*it)->GetConfig()->GetPixelID();
	      const Calibration::PixelNSB *NeighbourNSB = TelNSB->GetMonitor((*it));
	      uss_nsb_neighbours.Fill(NeighbourNSB->GetNSBValue());
	      
	      std::cout << neighbourID << ", ";
	    }
	    std::cout<< " with average NSB " << uss_nsb_neighbours.GetMean()<< " MHz"<< std::endl;
	  }
	}*/
	++pix_count;
      }     
      avg_tmp << uss_nsb_evt.GetMean() << std::endl;
      
      if ( ie == ds_nsb->GetEntries()) {
	  std::cout << "last measured time: "<< teltime_nsb.GetTime()<< std::endl;
	}
      
     }
  
    std::cout << ds_nsb->GetEntries() << std::endl;
    
    pix_count = 0;
    
    
    
    for (Sash::Pointer<Sash::Pixel> pix = (*tel)->beginPixel(), end_pix = (*tel)->endPixel(); pix != end_pix; ++pix) {
      if (pixel_avg_vector[pix_count]<=3) {
	const Sash::List<Sash::Pixel>& neighbours = pix->GetNeighbourList();
	int n_neighbours = 0;
	for(Sash::List<Sash::Pixel>::const_iterator it = neighbours.begin(); it != neighbours.end(); ++it) {
	  int neighbourID = (*it)->GetConfig()->GetPixelID();
	  pixel_avg_vector[pix_count] = pixel_avg_vector[pix_count] + pixel_avg_vector[neighbourID];
	  ++n_neighbours; 
	}
	pixel_avg_vector[pix_count] = pixel_avg_vector[pix_count]/n_neighbours;
	std::cout <<"Pixel " << pix_count +1 << " was averaged over " <<n_neighbours<< "to be "<<pixel_avg_vector[pix_count]/nevents<<  std::endl;
     
      
      }
      
      //std::cout << pix_count<<"bzw., " <<  pix->GetConfig()->GetPixelID() <<" " << ": "<< pixel_avg[pix_count] << " vs " << pixel_avg_vector[pix_count]/nevents << std::endl;
      avg_pix << pix_count << "; " << pixel_avg_vector[pix_count]/nevents<<std::endl;
      ++pix_count;
    }
    pixel_avg_vector.clear();
    // std::cout << "pixel_avg_vector: "<< pixel_avg_vector << std::endl;
    /*
	  if (PixNSB->GetNSBValue()<0) {
	    Utilities::SimpleStats uss_nsb_neighbours;
	    const Sash::List<Sash::Pixel>& neighbours = pix->GetNeighbourList();
	    R
	    for(Sash::List<Sash::Pixel>::const_iterator it = neighbours.begin(); it != neighbours.end(); ++it) {
	      int neighbourID = (*it)->GetConfig()->GetPixelID();
	      const Calibration::PixelNSB *NeighbourNSB = TelNSB->GetMonitor((*it));
	      uss_nsb_neighbours.Fill(NeighbourNSB->GetNSBValue());
	    }
	    pixel_avg[pix_count] = pixel_avg[pix_count] + uss_nsb_neighbours.GetMean()/ds_nsb->GetEntries();
	    if (ie==0) {
	    std::cout << pix_count <<" has -1, average of neighbours is "<< uss_nsb_neighbours.GetMean()<< std::endl;
	    }
	  }
     */

    
    
    
    avg_tmp.close();
    ss.str("");
    avg_pix.close();
    ss_avg_pix.str("");
    delta_pix.close();
    ss_delta_pix.str("");
    high_pix.close();
    ss_high_pix.str("");
    roll_avg.close();
    ss_roll_avg.str("");

    std::cout << "CT" << telId << " Average Camera NSB is " << uss_nsb.GetMean() << " MHz (RMS: " << uss_nsb.GetRMS() << ")" << std::endl;
    /* std::cout << "Average pixel values are ";
    for (int i = 0; i<size; i++) {
      std::cout << i << ", " << pixel_avg[i] << std::endl;
      }*/
    
  }
  }
}
