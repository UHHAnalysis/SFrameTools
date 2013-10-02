#include "SFrameTools/include/AnalysisModuleRunner.h"
#include "SFrameTools/include/identifier.h"
#include "SFrameTools/include/registry.h"

#include "TH1.h"
#include "TTree.h"
#include "TDataType.h"
#include "TEmulatedCollectionProxy.h"

#include <stdexcept>

using namespace std;

namespace{
    
    // see TTree.cxx:
    char DataTypeToChar(EDataType datatype){

    switch(datatype) {
        case kChar_t:     return 'B';
        case kUChar_t:    return 'b';
        case kBool_t:     return 'O';
        case kShort_t:    return 'S';
        case kUShort_t:   return 's';
        case kCounter:
        case kInt_t:      return 'I';
        case kUInt_t:     return 'i';
        case kDouble_t:
        case kDouble32_t: return 'D';
        case kFloat_t:
        case kFloat16_t:  return 'F';
        case kLong_t:     return 0; // unsupported
        case kULong_t:    return 0; // unsupported?
        case kchar:       return 0; // unsupported
        case kLong64_t:   return 'L';
        case kULong64_t:  return 'l';

        case kCharStar:   return 'C';
        case kBits:       return 0; //unsupported

        case kOther_t:
        case kNoType_t:
        default:
            return 0;
    }
    }
}


SFrameContext::SFrameContext(SCycleBase & base_, string event_treename_): base(base_), event_treename(std::move(event_treename_)){
    vector<pair<string, string> > props = base.GetConfig().GetProperties();
    for(vector<pair<string, string> >::const_iterator it = props.begin(); it != props.end(); ++it){
        settings[it->first] = it->second;
    }
    if(settings.size() != props.size()){
        throw invalid_argument("duplicate key in Sframe settings");
    }
}

std::string SFrameContext::get_setting(const std::string & key) const{
    map<string, string>::const_iterator it = settings.find(key);
    if(it==settings.end()) throw runtime_error("SFrameContext: did not find setting with key '" + key + "'");
    return it->second;
}
    
void SFrameContext::set_setting(const string & key, const string & value){
    settings[key] = value;
}
    
TDirectory * SFrameContext::get_tmpdir(){
    // see SCycleBaseHist::GetTempDir (which is private ...):
    static TDirectory* tempdir = 0;

    if(!tempdir){
        gROOT->cd();
        tempdir = gROOT->mkdir( "SFrameTempDir" );
        if(!tempdir){
            throw runtime_error( "SFrameTempDir could not be created" );
        }
    }
    return tempdir;
}

void SFrameContext::put(const identifier & id, TH1 * t){
    TDirectory * tempdir = get_tmpdir();
    
    // find out directory name; this is the part of the full name up to (excluding) the last '/'
    string fullname = id.name();
    
    // split full name into directory name (everything before last '/') and actual name (eveything beyond last '/'):
    string dirname, objname;
    size_t pos = fullname.rfind('/');
    if(pos != string::npos){
        dirname = fullname.substr(0, pos);
        objname = fullname.substr(pos+1);
    }else{
        objname = fullname;
    }
    
    // see implementation of template SCycleBaseHist::Book with inFile=false (which is the default).
    // instead of calling t->Clone, try to call SetDirectory. This is not a method of TObject,
    // so for now, make sure it works for the common cases of a histogram and tree:
    tempdir->cd();
    t->SetDirectory(tempdir);
    t->SetName(objname.c_str());
    SCycleOutput* out = new SCycleOutput(t, fullname.c_str(), dirname.c_str());
    TList * outlist = base.GetHistOutput();
    outlist->AddLast(out);
    gROOT->cd();
    hists[id] = t;
}


namespace{

// TTree::Branch is really strange and does not let me do in a simple
// way what I want to do: given a pointer to T, setup a branch. While this is possible
// with templates, it's not if the type is only known at runtime (like here), so
// for this case, we have to so some painful stuff by readin root code ...
void tree_branch(TTree * tree, const char * name, void * addr, void ** addraddr, const type_info & ti){
    //emulate
    // template<typename T>
    // outtree->Branch(name, (T*)addr);
    // which expands to
    // outtree->BranchImpRef(name, TBuffer::GetClass(ti), TDataType::GetType(ti), const_cast<void*>(addr), 32000, 99),
    // so look there to understand how this code works.
    TClass * class_ = TBuffer::GetClass(ti);
    if(class_){
        TClass * actualclass =  class_->GetActualClass(addr);
        assert(actualclass!=0);
        assert((class_ == actualclass) || actualclass->InheritsFrom(class_));
        assert(dynamic_cast<TEmulatedCollectionProxy*>(actualclass->GetCollectionProxy())==0); // if this fails, it means T is a STL without compiled STL dictionary
        // I'd like to do now:
        //  outtree->BronchExec(name, actualClass->GetName(), addr, *** kFALSE ****, bufsize, splitlevel);
        // but root doesn't let me, so work around by passing a pointer-to-pointer after all:
        tree->Bronch(name, actualclass->GetName(), addraddr);
    }
    else{
        EDataType dt = TDataType::GetType(ti);
        if(dt==kOther_t or dt==kNoType_t) throw invalid_argument("unknown class");
        char c = DataTypeToChar(dt);
        assert(c!=0);
        tree->Branch(name, addr, (string(name) + "/" + c).c_str());
    }
}
    
    
}


// event tree i/o
void SFrameContext::do_declare_event_input(const char * name, void * addr, const type_info & ti){
    
}

void SFrameContext::do_declare_event_output(const char * name, const void * caddr, const type_info & ti){
    TTree * outtree = base.GetOutputTree(event_treename.c_str());
    assert(outtree);
    void * addr = const_cast<void*>(caddr);
    ptrs.push_back(addr);
    tree_branch(outtree, name, addr, &ptrs.back(), ti);
}

SFrameContext::~SFrameContext(){
}
    
// other tree output:
void SFrameContext::do_declare_output(const identifier & tree_id, const char * name, const void * caddr, const type_info & ti){
    auto it = output_trees.find(tree_id);
    TTree * tree;
    if(it==output_trees.end()){
        tree = base.GetOutputMetadataTree(tree_id.name().c_str());
        assert(tree);
        output_trees[tree_id] = tree;
    }
    else{
        tree = it->second;
    }
    void * addr = const_cast<void*>(caddr);
    ptrs.push_back(addr);
    tree_branch(tree, name, addr, &ptrs.back(), ti);
}


void SFrameContext::write_output(const identifier & tree_id){
    auto it = output_trees.find(tree_id);
    if(it==output_trees.end()){
        throw runtime_error("called write_output for tree '" + tree_id.name() + "' which has not been declared previously");
    }
    it->second->Fill();
}


AnalysisModuleRunner::AnalysisModuleRunner(){
    DeclareProperty( "JetCollection", m_JetCollection );
    DeclareProperty( "GenJetCollection", m_GenJetCollection );
    DeclareProperty( "ElectronCollection", m_ElectronCollection );
    DeclareProperty( "MuonCollection", m_MuonCollection );
    DeclareProperty( "TauCollection", m_TauCollection );
    DeclareProperty( "PhotonCollection", m_PhotonCollection );
    DeclareProperty( "PrimaryVertexCollection", m_PrimaryVertexCollection );
    DeclareProperty( "METName", m_METName );
    DeclareProperty( "TopJetCollection", m_TopJetCollection );
    DeclareProperty( "TopJetCollectionGen", m_TopJetCollectionGen );
    DeclareProperty( "PrunedJetCollection", m_PrunedJetCollection );
    DeclareProperty( "GenParticleCollection", m_GenParticleCollection);
    DeclareProperty( "PFParticleCollection", m_PFParticleCollection);
    DeclareProperty( "readTTbarReco", m_readTTbarReco);
    DeclareProperty( "readCommonInfo", m_readCommonInfo);
    
    DeclareProperty("AnalysisModule", module_classname);
    DeclareProperty("OutputSelection", m_selection_output);
}

void AnalysisModuleRunner::BeginInputData( const SInputData& in ) throw( SError ){
    EventCalc* calc = EventCalc::Instance();
    calc->SetBaseCycleContainer(&m_bcc);
    analysis = Registry<AnalysisModule>::instance().build(module_classname);
    assert(analysis.get());
    selid_output = m_selection_output;
    context.reset(new SFrameContext(*this, "AnalysisTree"));
    assert(analysis.get());
    analysis->begin_dataset(*context);
    first_event_inputdata = true;
    // note: output tree setup is defered after processing the first event in the input dataset, see ExecuteEvent.
}

void AnalysisModuleRunner::BeginInputFile( const SInputData& ) throw( SError ){
    m_bcc.reset();

    if(m_ElectronCollection.size()>0) ConnectVariable( "AnalysisTree", m_ElectronCollection.c_str(), m_bcc.electrons);
    if(m_MuonCollection.size()>0) ConnectVariable( "AnalysisTree", m_MuonCollection.c_str(), m_bcc.muons);
    if(m_TauCollection.size()>0) ConnectVariable( "AnalysisTree", m_TauCollection.c_str(), m_bcc.taus);
    if(m_JetCollection.size()>0) ConnectVariable( "AnalysisTree", m_JetCollection.c_str(), m_bcc.jets);
    if(m_addGenInfo && m_GenJetCollection.size()>0) ConnectVariable( "AnalysisTree", m_GenJetCollection.c_str(), m_bcc.genjets);
    if(m_PhotonCollection.size()>0) ConnectVariable( "AnalysisTree", m_PhotonCollection.c_str(), m_bcc.photons);
    if(m_METName.size()>0) ConnectVariable( "AnalysisTree", m_METName.c_str(), m_bcc.met);
    if(m_PrimaryVertexCollection.size()>0) ConnectVariable( "AnalysisTree", m_PrimaryVertexCollection.c_str() , m_bcc.pvs);
    if(m_TopJetCollection.size()>0) ConnectVariable( "AnalysisTree", m_TopJetCollection.c_str(), m_bcc.topjets);
    if(m_addGenInfo && m_TopJetCollectionGen.size()>0) ConnectVariable( "AnalysisTree", m_TopJetCollectionGen.c_str() , m_bcc.topjetsgen);
    if(m_PrunedJetCollection.size()>0) ConnectVariable( "AnalysisTree", m_PrunedJetCollection.c_str() , m_bcc.prunedjets);
    if(m_addGenInfo && m_GenParticleCollection.size()>0) ConnectVariable( "AnalysisTree", m_GenParticleCollection.c_str() , m_bcc.genparticles);
    if(m_PFParticleCollection.size()>0) ConnectVariable( "AnalysisTree", m_PFParticleCollection.c_str() , m_bcc.pfparticles);
    if(m_addGenInfo && m_readCommonInfo) ConnectVariable( "AnalysisTree", "genInfo" , m_bcc.genInfo);
    if(m_readTTbarReco) ConnectVariable( "AnalysisTree", "recoHyps", m_bcc.recoHyps);

    ConnectVariable( "AnalysisTree", "run" , m_bcc.run);
    ConnectVariable( "AnalysisTree", "luminosityBlock" , m_bcc.luminosityBlock);
    ConnectVariable( "AnalysisTree" ,"event" ,m_bcc.event);
    
    ConnectVariable( "AnalysisTree", "rho" , m_bcc.rho);
    ConnectVariable( "AnalysisTree" ,"isRealData", m_bcc.isRealData);

    if(m_readCommonInfo){
        ConnectVariable( "AnalysisTree", "triggerResults" , m_bcc.triggerResults);
        ConnectVariable( "AnalysisTree", "triggerNames" , m_bcc.triggerNames);
        ConnectVariable( "AnalysisTree" ,"beamspot_x0", m_bcc.beamspot_x0);
        ConnectVariable( "AnalysisTree" ,"beamspot_y0", m_bcc.beamspot_y0);
        ConnectVariable( "AnalysisTree" ,"beamspot_z0", m_bcc.beamspot_z0);
    }
}

namespace {
    
template<typename T>                                                                                                                                                                                                      
void setup_branch(TTree * tree, const char * branchname, T * t){
    if(!t) return;
    tree->Branch(branchname, t);
}

}

void AnalysisModuleRunner::setup_output(){
    TTree * outtree = 0;
    try{
      outtree = GetOutputTree("AnalysisTree");
    }
    catch(const SError &){
         // if output tree is not found, we don't write anything.
        return;
    }
    assert(outtree!=0);

    // We want to write everything we have.
    // a) For all pointers in the base-cycle container, this means we want to write everyhing that is not 0. Note 
    //    that this rule covers the case of writing recoHyps which are not present int the input but were created by an
    //    AnalysisModule 'on-the-fly'
    // b) For all "plain" data, use the configuration to see whether it was read and only write what has been read.
    //
    // Note that we might write more than we read in some circumstances (especially for recoHyps): if
    // an AnalysisModule produced some data not there in the input, we still want to write it. So for all
    // pointers, create a branch fro each non-0 pointer:

    
    // a.:
    setup_branch(outtree, m_ElectronCollection.c_str(), m_bcc.electrons);
    setup_branch(outtree, m_MuonCollection.c_str(), m_bcc.muons);
    setup_branch(outtree, m_TauCollection.c_str(), m_bcc.taus);
    setup_branch(outtree, m_JetCollection.c_str(), m_bcc.jets);
    setup_branch(outtree, m_GenJetCollection.c_str(), m_bcc.genjets);
    setup_branch(outtree, m_PhotonCollection.c_str(), m_bcc.photons);
    setup_branch(outtree, m_METName.c_str(), m_bcc.met);
    setup_branch(outtree, m_PrimaryVertexCollection.c_str() , m_bcc.pvs);
    setup_branch(outtree, m_TopJetCollection.c_str(), m_bcc.topjets);
    setup_branch(outtree, m_TopJetCollectionGen.c_str() , m_bcc.topjetsgen);
    setup_branch(outtree, m_PrunedJetCollection.c_str() , m_bcc.prunedjets);
    setup_branch(outtree, m_GenParticleCollection.c_str() , m_bcc.genparticles);
    setup_branch(outtree, m_PFParticleCollection.c_str() , m_bcc.pfparticles);
    setup_branch(outtree, "genInfo" , m_bcc.genInfo);
    setup_branch(outtree, "recoHyps", m_bcc.recoHyps);
    setup_branch(outtree, "triggerResults", m_bcc.triggerResults);
    setup_branch(outtree, "triggerNames", m_bcc.triggerNames);

    // b.:
    // these are always read:
    setup_branch(outtree, "run", &m_bcc.run);
    setup_branch(outtree, "luminosityBlock", &m_bcc.luminosityBlock);
    setup_branch(outtree, "event", &m_bcc.event);
    setup_branch(outtree, "rho", &m_bcc.rho);
    setup_branch(outtree, "isRealData", &m_bcc.isRealData);
    
    if(m_readCommonInfo){
        setup_branch(outtree, "beamspot_x0", &m_bcc.beamspot_x0);
        setup_branch(outtree, "beamspot_y0", &m_bcc.beamspot_y0);
        setup_branch(outtree, "beamspot_z0", &m_bcc.beamspot_z0);
    }
}



void AnalysisModuleRunner::ExecuteEvent( const SInputData&, Double_t ) throw( SError ){
    EventCalc * ec = EventCalc::Instance();
    ec->Reset();
    assert(context.get()!=0);
    analysis->process(*ec, *context);
    
    // note that it is important to do the output setup *after* the event has been processed
    // as only then we know which information in the event is available (calculated), and
    // we want to write it all ...
    if(first_event_inputdata){
        first_event_inputdata = false;
        setup_output();
    }
    
    // prevent writing events not selected:
    if(!ec->selection_passed(selid_output)){
        throw SError( SError::SkipEvent );
    }
}

void AnalysisModuleRunner::EndMasterInputData(const SInputData & d) throw (SError){
    TList * l = GetHistOutput();
    TIter next(l);
    TObject * obj;
    
    // populate cutflows:
    map<string, pair<TH1D*, TH1D*> > cutflows; // first: weighted, second: unweighted
    while((obj = next())){
        string name(obj->GetName());
        if(name.find("cf_") != 0) continue;
        SCycleOutput * out = dynamic_cast<SCycleOutput*>(obj);
        if(!out) continue;
        obj = out->GetObject();
        TH1D * cutflow = dynamic_cast<TH1D*>(obj);
        if(!cutflow)continue;
        // get the selection name by cutting off "cf_" at the start:
        string sel_name = name.substr(3);
        // cut off "_raw in the end:"
        bool raw = false;
        if(name.size() > 4 and name.find("_raw")==name.size()-4){
            sel_name = sel_name.substr(0, sel_name.size()-4);
            raw = true;
        }
        if(raw){
            cutflows[sel_name].second = cutflow;
        }
        else{
            cutflows[sel_name].first = cutflow;
        }
    }
    
    // double to string; long uint to string:
    auto d2s = [](double d)-> string{
        char s[20];
        snprintf(s, 20, "%.6g", d);
        return s;
    };
    auto ul2s = [](unsigned long i)-> string{
        char s[20];
        snprintf(s, 20, "%lu", i);
        return s;
    };
    
    // print them:
    for(auto & it : cutflows){
        TH1D * cf = it.second.first;
        TH1D * cf_raw = it.second.second;
        if(cf==0 or cf_raw == 0 or cf->GetNbinsX() != cf_raw->GetNbinsX()){
            m_logger << WARNING << " did not find all cutflows (or inconsistent cutflows)" << endl;
            continue;
        }
        cout << endl << "Cutflow for selection '" << it.first << "':" << endl;
        TableOutput out({"Selection", "N_raw", "N_weighted"});
        TAxis * xax = cf->GetXaxis();
        for(int ibin=1; ibin<=cf->GetNbinsX(); ++ibin){
            out.add_row({xax->GetBinLabel(ibin), ul2s(cf_raw->GetBinContent(ibin)), d2s(cf->GetBinContent(ibin))});
        }
        out.print(cout);
    }
}

ClassImp(AnalysisModuleRunner);

