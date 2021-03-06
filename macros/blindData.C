#include <map>
#include <string>
#include <vector>
#include <iostream>

#include <TKey.h>
#include <TH1F.h>
#include <TFile.h>
#include <TROOT.h>
#include <TMath.h>
#include <TRandom3.h>
#include <TRegexp.h>
#include <TString.h>
#include <TCollection.h>

#include "HiggsAnalysis/HiggsToTauTau/macros/Utils.h"


/**
   \class   blindData blindData.C "HiggsAnalysis/HiggsToTauTau/macros/blindData.C"

   \brief   macro to replace data_obs by the exact expectation from all indicated background and potentially signal histograms.

   Macro to replace the data_obs histogram by the exact expectation from all indicated 
   background and potantially signal histograms. The list of background and signal 
   histograms and the input file have to be given as agruments. The macro will survery
   the input file for one directory layer according to the conventions for htt limit 
   combinations and in this directory structure search for the data_obs histogram and 
   the indicated background (and signal) historams. The data_obs histogram will be 
   reset and refilled by the sum of all indicated background (and signal) histograms.
   The bin-wise uncertainties will be set up the sqrt(N) of the corresponding bin 
   content to give for proper uncertainties, when the signal strength is calculated 
   based on the blinded data_obs histogram.

   The macro can be used to inject signal for a given mass into the data_obs. To obtain 
   the expected signal only from background just pass an empty string for the signal 
   samples as argument. In case you intend to inject signal keep in mind that you should 
   have scaled the signal of your input file by the corresponding cross section to obtain 
   a proper yield estimate for the signal. 

   The input od data_obs can be randomized according to a Poissonian distributrion. If 
   no signal nor background samples are given as function arguments the content of data_obs 
   will be randomized based on the existing content. This can be used to do toy studies 
   based on the actual observation in data. 

   Function arguments are:
   
   filename               : input file that is supposed to keep the input histograms
   background_patterns    : list of background samples (should be separated by comma 
                            or whitespaces).
   signal_patterns        : list of signal samples (should be separated by comma or 
                            whitespaces).
   directory_patterns     : give a list of directory\ies to which the blinding should be 
                            applied only or '*' to apply it to all directories. 
   armed                  : write modified data_obs histograms back to file.
   rnd                    : randomized the content of data_obs with given integer value 
                            as random seed. If -1 is given data_obs is not randomized. 
   signal_scale           : in case you want to add signal into the content of data_obs 
                            you can give an arbitary scale here. 1 corresponds to the 
			    SM case. 
   outputLabel            : if non empty and armed is true the output of data_obs (only!)
                            will be written to a new otuput file with postfix outputLabel. 
			    If armed and outputLabel is an empty string the inputfile 
			    will be updated with the new data_obs histogram(s).
   debug                  : invoke several debug output levels. 
*/

void adjustUncerts(TH1F* hist){
  for(int idx=0; idx<hist->GetNbinsX(); ++idx){
    hist->SetBinError(idx+1, TMath::Sqrt(hist->GetBinContent(idx+1)));    
  }
}

void randomize(TH1F* hist, unsigned int seed, unsigned int debug=0.)
{
  TRandom3* rnd = new TRandom3(seed); rnd->SetSeed();
  for(int idx=0; idx<hist->GetNbinsX(); ++idx){
    if(debug>0){
      std::cerr << "[" << idx+1 << "] : " << "mean=" << hist->GetBinContent(idx+1) << "  rnd=" << rnd->Poisson(hist->GetBinContent(idx+1)) << std::endl;  
    }
    float value = rnd->Poisson(hist->GetBinContent(idx+1));
    hist->SetBinContent(idx+1, value); hist->SetBinError(idx+1, TMath::Sqrt(value));
  }
  // Make sure there is no rounding error, and the total is truly an integer.
  hist->Scale(TMath::Nint(hist->Integral())/hist->Integral());
  delete rnd;
}

bool inPatterns(const std::string& test, const char* patterns)
{
  std::vector<std::string> samples;
  string2Vector(cleanupWhitespaces(patterns), samples);
  for(std::vector<std::string>::const_iterator sample = samples.begin(); sample!=samples.end(); ++sample){
    TRegexp matcher(sample->c_str(), true);  // true = make it "glob-style" i.e. "*" instead of ".*"
    if(TString(test).Index(matcher) > -1){
      return true;
    }
  }
  return false;
}

void blindData(const char* filename, const char* background_patterns="Fakes, EWK, ttbar, Ztt", const char* signal_patterns="ggH125, qqH125, VH125", const char* directory_patterns="*", bool armed=false, int rnd=-1, float signal_scale=1., const char* outputLabel="", unsigned int debug=1)
{
  /// prepare input parameters
  std::vector<std::string> signals;
  string2Vector(cleanupWhitespaces(signal_patterns), signals);
  std::vector<std::string> samples;
  string2Vector(cleanupWhitespaces(background_patterns), samples);
  samples.insert(samples.end(), signals.begin(), signals.end());

  // in case data_obs is supposed to be written to an extra output file
  // open the file, otherwise the data_obs in the input file will be
  // overwritten
  TFile* outputFile = 0; 
  if(!std::string(outputLabel).empty()){
    std::string out = std::string(filename); 
    outputFile = new TFile((out.substr(0, out.rfind("."))+"_"+outputLabel+".root").c_str(), "update"); 
  }

  TKey* idir;
  TH1F* buffer = 0;
  TH1F* blind_data_obs = 0;
  TFile* file = new TFile(filename, "update");
  TIter nextDirectory(file->GetListOfKeys());
  while((idir = (TKey*)nextDirectory())){
    if( idir->IsFolder() ){
      file->cd(); // make sure to start in directory head 
      if( debug>0 ){ std::cerr << "Found directory: " << idir->GetName() << std::endl; }
      // check if we want to muck w/ this directory. For the vhtt case, we 
      // have different background types in the same root file, so we have 
      // to run blindData twice.
      if (!inPatterns(std::string(idir->GetName()), directory_patterns)) {
        if( debug>0 ){ 
	  std::cerr << "WARNING: Skipping directory: " << idir->GetName() << std::endl;
	  std::cerr << "         No match found in pattern: " << directory_patterns << std::endl; 
	}
        continue;
      }
      if( file->GetDirectory(idir->GetName()) ){
	file->cd(idir->GetName()); // change to sub-directory
	buffer = (TH1F*)file->Get((std::string(idir->GetName())+"/data_obs").c_str());
	if(!buffer){
	  std::cout << "WARNING: Did not find histogram data_obs in directory: " << idir->GetName() << std::endl;
	  std::cout << "         Will skip directory: " << std::endl;
	  continue;
	}
        blind_data_obs = (TH1F*)buffer->Clone("data_obs"); 
	if(!samples.empty()){
	  blind_data_obs->Reset();
	  std::cout << "INFO  : Blinding datacads now." << std::endl;
	  for(std::vector<std::string>::const_iterator sample = samples.begin(); sample!=samples.end(); ++sample){
	    if( debug>0 ){ std::cerr << "Looking for histogram: " << (std::string(idir->GetName())+"/"+(*sample)) << std::endl; }
	    // add special treatment for et/mt ZLL,ZJ,ZL here. You can run the
	    // macro with ZLL, ZL, ZJ in the background samples. Those samples,
	    // which do not apply for one or the other event category are 
	    // skipped here.
	    if(std::string(idir->GetName()).find("vbf") != std::string::npos && (*sample == std::string("ZL") || *sample == std::string("ZJ"))){
	      continue;
	    }
	    else if(std::string(idir->GetName()).find("vbf") == std::string::npos && *sample == std::string("ZLL")){
	      continue;
	    }
	    buffer = (TH1F*)file->Get((std::string(idir->GetName())+"/"+(*sample)).c_str()); 
	    if (!buffer) {
	      std::cerr << "ERROR : Could not get histogram from: " << std::string(idir->GetName())+"/"+(*sample) << std::endl;
	      std::cerr << "        Histogram will be skipped   : " << std::string(idir->GetName())+"/"+(*sample) << std::endl;
	      continue;
	    }
	    if(inPatterns(*sample, signal_patterns)) {
	      if( debug>1 ){
		std::cerr << "INFO  : Scale signal sample " << *sample << " by scale " << signal_scale << std::endl;
	      }
	      buffer->Scale(signal_scale);
	    }
	    blind_data_obs->Add(buffer);
	    if (debug > 1){
	      std::cerr << "INFO  : Adding: " << buffer->GetName() << " -- " << buffer->Integral() << " --> New value: " << blind_data_obs->Integral() << std::endl;
	    }
	  }
	}
	else{
	  std::cout << "INFO  : Data are not blinded." << std::endl;
	}
	if(rnd>=0){
	  // randomize histogram; this will automatically have integer integral
	  std::cerr << "-- R A N D O M I Z I N G --" << std::endl;
	  randomize(blind_data_obs, rnd, debug);
	}
	else{
	  // use expected mean with signal injected
	  blind_data_obs->Scale(TMath::Nint(blind_data_obs->Integral())/blind_data_obs->Integral());
	  // adjust uncertaintie
	  adjustUncerts(blind_data_obs);
	}
	std::cout << "INFO  : New data_obs yield: " << idir->GetName() << "   " << TMath::Nint(blind_data_obs->Integral()) << std::endl;
	if(armed){
          if (debug > 1){
            std::cerr << "INFO  : Writing to file: " << blind_data_obs->GetName() << std::endl;
	  }
	  if(outputFile){
	    // write to a dedicated new file with name output in case output has been specified
	    outputFile->mkdir(idir->GetName()); outputFile->cd(idir->GetName());
	    blind_data_obs->Write("data_obs");
	  }
	  else{
	    // override old data_obs in the inputfile otherwise
	    file->cd(idir->GetName());
	    blind_data_obs->Write("data_obs", TObject::kOverwrite); 
	  }
	}
      }
    }
  }
  file->Close();
  if(outputFile){
    outputFile->Close();
  }
  return;
}

