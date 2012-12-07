#ifndef MCDataScaleFactors_H
#define MCDataScaleFactors_H

#include "include/Utils.h"
#include "include/EventCalc.h"

/**
 *  @short module to apply data-MC lepton scale factors for trigger and ID
 *
 *  
 */
class LeptonScaleFactors{
 public:
 /**
  * constructor
  *
  * first argument: list of corrections to be applied together with weight factors for luminosity, example: "MuonRunA 1.5 MuonRunB 2.6 MuonRunC 7.8"
  *
  * second argument: systematic shift
  * @see E_SystShift
  */
  LeptonScaleFactors(std::vector<std::string> correctionlist, E_SystShift syst_shift=e_Default);
  ///Default destructor
  ~LeptonScaleFactors(){};

  ///return the weighted correction factor
  double GetWeight();

 private:
  E_SystShift m_syst_shift;
  std::vector<std::pair<std::string, double> > m_correctionlist;

};


class BtagFunction
{
 public:
  BtagFunction();
  
  virtual ~BtagFunction()
    {
    };
  
  virtual float value(const float &x) const = 0;
  virtual float value_plus(const float &x) const = 0;
  virtual float value_minus(const float &x) const = 0;
  
 protected:
  const unsigned int find_bin(const float &jet_pt) const;
  
 private:
  std::vector<float> _bins;
};

class BtagScale: public BtagFunction
{
  // CSVT operating point
 public:
  BtagScale();
  
  virtual float value(const float &jet_pt) const;
  virtual float value_plus(const float &jet_pt) const
  {
                return value(jet_pt) + error(jet_pt);
  }
  
  virtual float value_minus(const float &jet_pt) const
  {
    const float value_ = value(jet_pt) - error(jet_pt);
    
    return value_ > 0 ? value_ : 0;
  }
  
 protected:
  virtual float error(const float &jet_pt) const;
  
 private:
  std::vector<float> _errors;
};

class CtagScale: public BtagScale
{
 protected:
  virtual float error(const float &jet_pt) const;
};

class LightScale: public BtagFunction
{
 public:
  virtual float value(const float &jet_pt) const;
  virtual float value_plus(const float &jet_pt) const;
  virtual float value_minus(const float &jet_pt) const;
  
 private:
  float value_max(const float &jet_pt) const;
  float value_min(const float &jet_pt) const;
};

class BtagEfficiency: public BtagFunction
{
  // Errors are not provided ... yet
 public:
  BtagEfficiency();
  
  virtual float value(const float &jet_pt) const;
  virtual float value_plus(const float &jet_pt) const
  {
    return value(jet_pt);
  }
  
  virtual float value_minus(const float &jet_pt) const
  {
    return value(jet_pt);
  }
  
 private:
  std::vector<float> _values;
};

class CtagEfficiency: public BtagEfficiency
{
};

class LightEfficiency: public BtagEfficiency
{
  // Errors are not provided ... yet
 public:
  LightEfficiency();
  
  virtual float value(const float &jet_pt) const;
  
 private:
  std::vector<float> _values;
    };

class LightEfficiencyData: public BtagEfficiency
{
  // Errors are not provided ... yet
 public:
  virtual float value(const float &jet_pt) const;
};




/**
 *  @short module to apply data-MC scale factors for b tagging
 *
 *  
 */
class BTaggingScaleFactors{
 public:
 /**
  * constructor
  *
  * second argument: systematic shift
  * @see E_SystShift
  */
  BTaggingScaleFactors(E_SystShift syst_shift=e_Default);
  ///Default destructor
  ~BTaggingScaleFactors(){};

  ///return the weighted correction factor
  double GetWeight();

 private:

  E_SystShift m_syst_shift;

  float scale(const bool &is_tagged,
	      const float &jet_pt,
	      const BtagFunction* sf,
	      const BtagFunction* eff,
	      const E_SystShift &sytematic);
  
  float scale_data(const bool &is_tagged,
		   const float &jet_pt,
		   const BtagFunction* sf,
		   const BtagFunction* eff,
		   const E_SystShift &sytematic);


  BtagFunction* _scale_btag;
  BtagFunction* _eff_btag;
  
  BtagFunction* _scale_ctag;
  BtagFunction* _eff_ctag;
  
  BtagFunction* _scale_light;
  BtagFunction* _eff_light;
};




#endif
