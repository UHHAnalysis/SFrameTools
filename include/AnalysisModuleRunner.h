#ifndef ANALYSIS_MODULE_RUNNER_H
#define ANALYSIS_MODULE_RUNNER_H

#include <memory>

#include "core/include/SCycleBase.h"
#include "SFrameTools/include/AnalysisModule.h"

/** \brief The SFrame specific implementation for the Context
 * 
 * See the Context class for a description of the methods.
 */
class SFrameContext: public Context{
public:
    SFrameContext(SCycleBase & base, std::string event_treename);
    
    virtual std::string get_setting(const std::string & key) const;
    virtual void set_setting(const string & key, const string & value);
    
    virtual void put(const identifier & id, TH1 * t);
    
    ~SFrameContext();

private:
    
    // event tree i/o
    virtual void do_declare_event_input(const char * name, void * addr, const type_info & ti);
    virtual void do_declare_event_output(const char * name, const void * addr, const type_info & ti);
    
    // other tree output:
    virtual void do_declare_output(const identifier & tree_id, const char * name, const void * addr, const type_info & ti);
    virtual void write_output(const identifier & tree_id);
    
    TDirectory * get_tmpdir();
    
    SCycleBase & base;
    std::string event_treename;
    std::map<string, string> settings;
    std::map<identifier, TH1*> hists;
    
    std::list<void*> ptrs;
    
    std::map<identifier, TTree*> output_trees;
};


/** \brief The SFrame cycle used to run an AnalysisModule
 * 
 * It takes care of
 *  - preparing the input tree according to the configuration
 *  - constructing the AnalysisModule class
 *  - calling AnalysisModule::begin_dataset for each dataset and AnalysisModule::process for each event
 *  - writing the tree of selected events
 */
class AnalysisModuleRunner: public SCycleBase{
public:
    AnalysisModuleRunner();
    
    // called at the beginning of the cycle, but for proof only in one process!
    void BeginCycle() throw( SError ){}
    void EndCycle() throw( SError ){}

    // called at the beginning of the input data, on all proof nodes:
    void BeginInputData( const SInputData& ) throw( SError );
    void EndInputData  ( const SInputData& ) throw( SError ){}

    // called at the beginning of an input file, on the proof nodes:
    void BeginInputFile( const SInputData& ) throw( SError );
    void ExecuteEvent( const SInputData&, Double_t ) throw( SError );
    
    // called after processing the dataset, only on the proof master, not on the proof nodes:
    void EndMasterInputData(const SInputData &) throw (SError);
        
    ClassDef(AnalysisModuleRunner, 0);
  
private:
    void setup_output();
    
    std::string m_JetCollection, m_GenJetCollection, m_ElectronCollection, m_MuonCollection, 
      m_TauCollection, m_PhotonCollection, m_PrimaryVertexCollection, m_METName, m_TopJetCollection, m_TopTagJetCollection, m_HiggsTagJetCollection, m_TopJetCollectionGen,
       m_PrunedJetCollection, m_GenParticleCollection, m_PFParticleCollection;
    bool m_readTTbarReco, m_readCommonInfo, m_addGenInfo;
    
    BaseCycleContainer m_bcc;
    bool first_event_inputdata;
    
    std::auto_ptr<SFrameContext> context;
    
    std::string m_selection_output;
    identifier selid_output;
       
    // the actual analysis module to run:
    std::string module_classname;
    std::auto_ptr<AnalysisModule> analysis;
};

#endif
