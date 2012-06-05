// Dear emacs, this is -*- c++ -*-
#ifndef BaseHists_H
#define BaseHists_H

// STL include(s):
#include <map>
#include <string>

// ROOT include(s):
#include <TObject.h>
#include <TString.h>

// Local include(s):
//#include "include/ISCycleBaseHist.h"
#include "include/SCycleBaseHist.h"
#include "include/SCycleBaseBase.h"
#include "include/SError.h"

// Forward declaration(s):
class TDirectory;
class TH1;
class TList;

/**
 *   Base Class for the histogramming functionality
 *
 *   Class which books and fills histograms. This class is a base 
 *   class, all histogramming classes should inherit from it. 
 *   The booking and naming of the histograms is taken care of 
 *   through this class. Any derived class has to overwrite the 
 *   Init, Fill and Finish methods.
 *   The derived objects have to be initialised in each 
 *   SCycleBase::BeginInputData(...)
 *
 *   @version $Revision: 0 $
 */

class BaseHists : public SCycleBaseHist {


public:
   /// Named constructor
   BaseHists(const char* name);

   /// Default destructor
   ~BaseHists();

   virtual void Init() = 0;

   virtual void Fill() = 0;

   virtual void Finish() = 0;

   double* MakeLogBinning(int n_bins, double xmin, double xmax);

   TString GetName() {return m_name;}
   void SetName(TString name) {m_name = name;}

protected:
   /// Function placing a ROOT object in the output file
   template< class T > T* Book( const T& histo ) throw( SError );

   /// Function searching for a ROOT object in the output file
   template< class T > T* Retrieve( const char* name ) throw( SError );

   /// Function for persistifying a ROOT object to the output
   void WriteObj( const TObject& obj ) throw( SError );

   /// Function searching for 1-dimensional histograms in the output file
   TH1* Hist( const char* name );

private:
   // private constructor, use the named one
   BaseHists(){}

   TString m_name;


}; // class BaseHists

// Don't include the templated function(s) when we're generating
// a dictionary:
#ifndef __CINT__
#include "BaseHists.icc"
#endif



#endif // BaseHists_H
