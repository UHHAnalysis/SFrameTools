#include <cmath>
#include "include/BaseHists.h"

BaseHists::BaseHists(const char* name)
{
  // named constructor
  m_name = name;
  SetLogName("BaseHists_" + m_name);
}

BaseHists::~BaseHists()
{
  // default destructor, does nothing
}

void BaseHists::WriteObj( const TObject& obj ) throw( SError ) {

  SCycleBaseHist::WriteObj( obj, m_name.Data() );
  return;
}


TH1* BaseHists::Hist( const char* name ) {

  TH1* his = SCycleBaseHist::Hist( name, m_name.Data() );
  return his;
}



double* BaseHists::MakeLogBinning(int n_bins, double xmin, double xmax)
{
  // return pointer to an array of bins in logscale             
  // usage:  Double_t* ptbins = MakeLogBinning(50,100,2000.);  
  //        ... TH1F("Pt","Jet PT", 50, ptbins);       

   double *binning = new double[n_bins+1];
   
   double delta_x_log = (log(xmax)-log(xmin))/n_bins;
 
   binning[0]=xmin;
   // put equidistant binwidth on a logarithmic scale
   for(int i=1;i<=n_bins;++i){
     binning[i] = exp(log(binning[i-1])+ delta_x_log);
   }
   return binning; 
}
