#ifndef ANALYSIS_MODULE_RUNNER_H
#define ANALYSIS_MODULE_RUNNER_H

#include <memory>

#include "core/include/SCycleBase.h"
#include "SFrameTools/include/AnalysisModule.h"

/** \brief The SFrame specific implementation for the Context
 * 
 * See the Context class for a description of the methods.
 * 
 * This class is intended to be used from a SCycle; see AnalysisModuleRunner for
 * an example use:
 *  - create a new instance of SFrameContext in SCycle::BeginInputData and set the dataset.type and dataset.version
 *    according to the current dataset
 *  - call SFrameContext::begin_input_file from SCycle::BeginInputFile to make sure that input addresses are set up correctly.
 */
class SFrameContext: public Context{
public:
    SFrameContext(SCycleBase & base, const SInputData& sin);
    
    virtual void put(const identifier & id, TH1 * t);
    
    ~SFrameContext();
    
    // this should be called from SCycleBase::BeginInputFile. It takes care
    // of setting up the input addresses
    void begin_input_file();
    void begin_event();

private:
    
    // event tree i/o
    // note that the declare_input does nothing but saving the parameters for later; the actual setup of
    // input branches is done in begin_input_file.
    virtual void do_declare_event_input(const char * name, void * addr, const type_info & ti);
    virtual void do_declare_event_output(const char * name, const void * addr, const type_info & ti);
    
    // other tree output:
    virtual void do_declare_output(const identifier & tree_id, const char * name, const void * addr, const type_info & ti);
    virtual void write_output(const identifier & tree_id);
    
    TDirectory * get_tmpdir();
    
    SCycleBase & base;
    std::string event_treename;
    
    // input (event) tree stuff:
    TTree * input_tree;
    // input event branch name to branch info:
    struct branchinfo {
        TBranch * branch;
        const type_info * ti; // this is always a non-pointer type.
        void * addr; // address of an object of type ti.
    };
    std::map<std::string, branchinfo> bname2bi;
    
    
    // output tree stuff:
    // pointers to the objects of the output tree(s)
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
    
    virtual void Initialize( TXMLNode* node ) throw( SError );
    virtual void SetConfig(const SCycleConfig& config);
        
    ClassDef(AnalysisModuleRunner, 0);
  
private:
    void setup_output();
    
#ifndef __CINT__ // we don't have anything for CINT to serialize ...:
    
    std::string m_JetCollection, m_GenJetCollection, m_ElectronCollection, m_MuonCollection, 
      m_TauCollection, m_PhotonCollection, m_PrimaryVertexCollection, m_METName, m_TopJetCollection, m_TopTagJetCollection, m_HiggsTagJetCollection, m_TopJetCollectionGen,
       m_PrunedJetCollection, m_GenParticleCollection, m_PFParticleCollection;
    bool m_readTTbarReco, m_readCommonInfo, m_addGenInfo;
    
    BaseCycleContainer m_bcc;
    bool first_event_inputdata;
    
    std::auto_ptr<SFrameContext> context;
    
    identifier selid_output;
       
    // the actual analysis module to run:
    std::auto_ptr<AnalysisModule> analysis;
    
    std::map<std::string, std::string> dummyConfigVars;
    
#endif
};

#endif
