// Dear emacs, this is -*- c++ -*-
#ifndef BaseHists_H
#define BaseHists_H

#include <map>
#include <string>

#include <TObject.h>
#include <TString.h>
#include <TList.h>

#include "core/include/SCycleBase.h"
#include "core/include/SError.h"
#include "SFrameTools/include/Objects.h"

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
 *   @version $Revision: 1.1 $
 */

class BaseHists : public SCycleBase {
public:
   explicit BaseHists(const char* name);

   virtual ~BaseHists();

   virtual void Init() = 0;

   virtual void Scale(double scale) = 0;

   virtual void Fill() = 0;
  virtual void Fill2(TopJet topjet, double mva_value) {}

   // This is called at the end of the histograms. IMPORTANT: it's best not to use this method, since
   // it's called on each proof worker node (in proof node), and thus things like taking ratios of histograms
   // might work in local mode, but do not make sense in proof mode, as the divided histograms are combined
   // by adding them ...
   virtual void Finish(){}

   double* MakeLogBinning(int n_bins, double xmin, double xmax);// DEPRECATED_MSG("use log_binning in Utils.h instead");

   TString GetName() {return m_name;}
   void SetName(const TString & name) {m_name = name;}

  // class has to inherit from SCycleBase to have access to histogram functionality in SFrame
  // implement virtual routines but throw exception in case they are called
   void BeginCycle()throw( SError ) {  m_logger << ERROR << "This should not happen: BeginCycle called for BaseHist class" << SLogger::endmsg; } ;
   void EndCycle()throw( SError ){  m_logger << ERROR << "This should not happen: EndCycle called for BaseHist class" << SLogger::endmsg; } ;
   void BeginInputData(const SInputData&)throw( SError ){  m_logger << ERROR << "This should not happen: BeginInputData called for BaseHist class" << SLogger::endmsg; } ;
   void EndInputData(const SInputData&)throw( SError ){  m_logger << ERROR << "This should not happen: EndInputData called for BaseHist class" << SLogger::endmsg; } ;
   void BeginInputFile(const SInputData&)throw( SError ){  m_logger << ERROR << "This should not happen: BeginInputFile called for BaseHist class" << SLogger::endmsg; } ;
   void ExecuteEvent(const SInputData&, Double_t)throw( SError ){  m_logger << ERROR << "This should not happen: ExecuteEvent called for BaseHist class" << SLogger::endmsg; } ;


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
