#include "include/PDFWeights.h"

PDFWeights::PDFWeights(E_SystShift syst_shift, TString pdfname, TString pdfweightdir) : m_logger( "PDFWeights" )
{
    m_syst_shift = syst_shift;

    if( ( gSystem->Load("libLHAPDF") )==-1){
      m_logger << ERROR << "libLHAPDF not found, no pdf weights will be applied. To apply pdf re-weighting, add path to libLHAPDF.so to LD_LIBRARY_PATH" << SLogger::endmsg;
      m_libvalid=false;
      return;
    }
    m_libvalid=true;


    LHAPDF::initPDFSet(1, (string)(pdfname+".LHgrid"));
    m_N_unc = LHAPDF::numberPDF();

    m_normalize_to_total_sum=false;
    if(pdfweightdir!=""){
      m_normalize_to_total_sum=true;
      
      TString filename = pdfweightdir;
      filename +="_";
      filename += pdfname;
      filename += "_weights.txt";

      m_logger << INFO << "Do pdf re-weighting with respect to weights in file: " << filename << SLogger::endmsg;

      ifstream infile (((string)filename).c_str());
      
      infile>>m_N_tot;

      do{
	double number;
	infile>>number;
	m_sumofweights.push_back(number);
      }while ( !infile.eof() );
      //erase last entry which is loaded twice
      m_sumofweights.erase(m_sumofweights.end()-1);

      m_logger << DEBUG << "total number of events before selection: " << m_N_tot <<  SLogger::endmsg;
      for(unsigned int i=0; i< m_sumofweights.size(); ++i){
	m_logger << DEBUG << "sum of weights for pdf set " << i+1 << ": " << m_sumofweights[i] << SLogger::endmsg;
      }

      infile.close();
    
      if(m_sumofweights.size()!=m_N_unc){
	m_logger << ERROR << "Number of event weights in input file ("<< m_sumofweights.size()<< ") != number of parameters of chosen pdf set ("<< m_N_unc<< ")"<< SLogger::endmsg;
      }
    }
}

std::vector<double> PDFWeights::GetWeightList(){
  std::vector<double> pdf_weights;

  if(!m_libvalid) return pdf_weights;

  //pdf weighting code taken from https://twiki.cern.ch/twiki/bin/view/CMS/TWikiTopRefSyst#PDF_uncertainties

  EventCalc* calc = EventCalc::Instance();

  double x1=calc->GetGenInfo()->pdf_x1();
  double x2=calc->GetGenInfo()->pdf_x2();
  
  int id1 = calc->GetGenInfo()->pdf_id1();
  int id2 = calc->GetGenInfo()->pdf_id2();

  double q = calc->GetGenInfo()->pdf_scalePDF();

  LHAPDF::usePDFMember(1,0);
  double xpdf1 = LHAPDF::xfx(1, x1, q, id1);
  double xpdf2 = LHAPDF::xfx(1, x2, q, id2);

  double w0 = xpdf1 * xpdf2;
  for(unsigned int i=1; i <=m_N_unc; ++i){
    LHAPDF::usePDFMember(1,i);
    double xpdf1_new = LHAPDF::xfx(1, x1, q, id1);
    double xpdf2_new = LHAPDF::xfx(1, x2, q, id2);
    double weight = xpdf1_new * xpdf2_new / w0;

    if(m_normalize_to_total_sum){
      pdf_weights.push_back(weight/m_sumofweights[i-1]*m_N_tot);
    }
    else{ 
      pdf_weights.push_back(weight);
    }
  }

  return pdf_weights;
  
}


double PDFWeights::GetWeight(unsigned int index){

  if(!m_libvalid) return 1.;

  std::vector<double> pdf_weights = GetWeightList();

  if(index>pdf_weights.size() || index<1){
    m_logger << ERROR << "PDF index "  << index << " out of range, should be >=1 and <= " << pdf_weights.size() << SLogger::endmsg;
    return 1.;
  }

  return pdf_weights.at(index-1);

}
