#ifndef ANALYSISMODULE_H
#define ANALYSISMODULE_H

#include "SFrameTools/include/EventCalc.h"
#include "SFrameTools/include/Selection.h"
#include "SFrameTools/include/boost_includes.h" // for static assert and type traits, ptr_vector, and lexical_cast

#include "core/include/SCycleBaseHist.h"
#include "core/include/SCycleBase.h"

#include "registry.h"
#include "identifier.h"

#include <memory>

class Context;

/**  \brief Abstract base class for all analysis modules, independent of SFrame
 * 
 * Using it comprises two parts: implementing the analysis logic and calling it within the event loop to do the actual work.
 * 
 * For the first part, derive from this class and implement the 'begin_dataset' and 'process' methods.
 * 
 * Calling 'begin_dataset' or 'process' can be done either
 *   - "manually" from another AnalysisModule or from a SFrame cycle
 *   - "automatically" by using the AnalysisModuleRunner SFrame cycle which does nothing but call these
 *     methods in the correct order, after setting up the input.
 */
class AnalysisModule {
public:
    
    virtual ~AnalysisModule();
    
    /** \brief Method called at the beginning of a dataset
     * 
     * Use the Context to 
     *  - access the configuration
     *  - create output
     */
    virtual void begin_dataset(Context & ctx) = 0;
    
    /** \brief Method called for each event in the dataset
     * 
     * Do the filling of histograms, calculating new variables, etc. here.
     */
    virtual void process(EventCalc & event, Context & ctx) = 0;
};


// for modules which should be directly runnable with the AnalysisModuleRunner, call this macro in a C++ file:
#define REGISTER_ANALYSIS_MODULE(T) REGISTER(T, AnalysisModule, T)


/** \brief Abstract utility class for configuration and I/O; passed to AnalysisModule::begin_dataset and AnalysisModule::process.
 * 
 * This class provides methods to put histograms or tree data in the output root file.
 * 
 * There are two types of trees: event trees and user-defined trees. Event trees are managed by the framework
 * in the sense that each entry in the input event tree is processed once (by calling AnalysisModule::process)
 * and each selected event is written to the output event tree.
 * User-defined trees, on the other hand, only exist as output trees and can contain arbitrary, non-event data. An example
 * would be to fill data on a per-object basis. To use user-defined trees, call the 'declare_output' method, and
 * for each entry, call the 'write_output' method.
 */
class Context{
public:
    
    /** \brief Get a setting from the configuration
     * 
     * Throws a runtime_error if the setting with this key is not available, i.e. iff has_setting(key)==false.
     */
    std::string get_setting(const std::string & key) const{
        std::map<std::string, std::string>::const_iterator it = settings.find(key);
        if(it==settings.end()) fail(key);
        return it->second;
    }
    
    /** \brief Get a setting from the configuration, returning a default value if it does not exist
     */
    std::string get_setting(const std::string & key, const std::string & def) const{
        std::map<std::string, std::string>::const_iterator it = settings.find(key);
        if(it==settings.end()) return def;
        else return it->second;
    }
    
    /** \brief Get all settings from the configuration
     */
    std::map<std::string, std::string> get_all_settings() const{
        return settings;
    }
    
    /** \brief Test whether a setting is available
     */
    bool has_setting(const std::string & key) const{
        return settings.find(key) != settings.end();
    }
    
    /** \brief Set a key to the given value
     *
     * Set a configuration value, which can be retrieved later with get_setting; can also be used to overwrite
     * an existing setting.
     */
    void set_setting(const std::string & key, const std::string & value){
        settings[key] = value;
    }
        
    /** \brief Put a histogram in the output root file at the specified path
     * 
     * By calling this routine, you pass memory management of t to the framework.
     * 
     * To create the object in a subdirectory of the output, use an id with the corresponding
     * full name, i.e. "dir/subdirB/name"; the required directories will be created automatically.
     */
    virtual void put(const identifier & id, TH1 * t) = 0;
    
    /** \brief Declare an input variable in the event tree
     * 
     * T must not be a pointer type. Make sure that the lifetime of t exceeds the reading of the input event tree.
     */
    template<typename T>
    void declare_event_input(const char * name, T & t){
        BOOST_STATIC_ASSERT_MSG(!boost::is_pointer<T>::value, "T must not be of pointer type");
        do_declare_event_input(name, static_cast<void*>(&t), typeid(T));
    }
    
    /** \brief Declare an output variable in the event tree
     * 
     * T must not be a pointer type. Make sure that the lifetime of t exceeds the writing of the last event to the output event tree.
     */
    template<typename T>
    void declare_event_output(const char * name, const T & t){
        BOOST_STATIC_ASSERT_MSG(!boost::is_pointer<T>::value, "T must not be of pointer type");
        do_declare_event_output(name, static_cast<const void*>(&t), typeid(T));
    }
    
    /** \brief Declare a non-event variable which should be written to the tree identified by tree_id instead.
     * 
     * T must not be of pointer type. Make sure that the lifetime of t exceeds the last call of write_output of the corresponding tree
     */
    template<typename T>
    void declare_output(const identifier & tree_id, const char * name, T & t){
        BOOST_STATIC_ASSERT_MSG(!boost::is_pointer<T>::value, "T must not be of pointer type");
        do_declare_output(tree_id, name, static_cast<const void*>(&t), typeid(T));
    }
    
    /// Write a single entry to the output tree identified by tree_id with the current contents of the declared output addresses.
    virtual void write_output(const identifier & tree_id) = 0;
    
    /// declare the destructor virtual, as this is a purely virtual base class.
    virtual ~Context();
    
private:
    
    virtual void do_declare_event_input(const char * name, void * addr, const type_info & ti) = 0;
    virtual void do_declare_event_output(const char * name, const void * addr, const type_info & ti) = 0;
    virtual void do_declare_output(const identifier & tree_id, const char * name, const void * t, const type_info & ti) = 0;
    
    void fail(const std::string & key) const{
        throw std::runtime_error("did not find setting '" + key + "'");
    }
    
    std::map<std::string, std::string> settings;
};


// alternative implementation of BaseHists, NOT inheriting from SCycleBase
//
// typical use: derive from this class and implement the constructor with the same signature as the one here.
// In the constructor, book the histograms to produce via the protected book template method.
//
// In your analysis module in AnalysisModule::begin_dataset, instantiate a number
// of 'Hists' (derived) classes, using the name of the selection as directory name.
// In your AnalysisModule::process routine, call the 'Hists::fill' method for each event passing the selection.
class Hists{
public:
    
    Hists(Context & ctx_, const string & dirname_): ctx(ctx_), dirname(dirname_){
        // check that the dirname does not start or end with '/':
        if(dirname.find('/')==0) throw std::runtime_error("SimpleBaseHists: provided directory starts with '/'; this is not allowed.");
        size_t p = dirname.rfind('/');
        if(p != string::npos and p + 1 == dirname.size()){
            throw std::runtime_error("SimpleBaseHists: provided directory ends with '/'; this is not allowed.");
        }
    }
    
    virtual ~Hists();
    
    // - fill the histograms based on the event content; this is called for each (selected) event.
    //   Use the get method with the same name as used in book to retrieve the histogram
    virtual void fill(EventCalc & ev) = 0;
    
protected:
    
    // in C++11, do this:
    //template<typename T, typename... targs>
    // book(...)
    // T * t = new T(name.c_str(), args...);
    
    void book(const identifier & id, TH1 * h){
        std::string name = id.name();
        if(name.find('/')!=string::npos){
            throw std::runtime_error(" name '" + name + "' illegal: '/' in histogram names not allowed in Hists (use Context::put directly for putting histograms in subdirectories)");
        }
        h->SetName(name.c_str());
        h->Sumw2();
        ctx.put(dirname + "/" + name, h);
        histos[id] = h;
    }
    
    TH1* get_hist(const identifier & id){
        std::map<identifier, TH1*>::iterator it = histos.find(id);
        if(it==histos.end()){
            throw std::runtime_error("Hists::get_hist: did not find histogram '" + id.name() + "'");
        }
        return it->second;
    }
    
private:
    Context & ctx;
    std::string dirname;
    std::map<identifier, TH1*> histos;
};


/** \brief Apply an logical AND of a given list of selection modules, make a cut flow, and store the result in the EventCalc container
 * 
 * This is a replacement for Selection, featuring:
 *  - simpler interface
 *  - clear memory ownership semantics
 *  - a cutflow which is correctly merged in PROOF mode
 * 
 * The basic usage is:
 *  - in AnalysisModule::begin_dataset of your analysis class, create an instance of AndSelection with some id, and call the 'begin_dataset' method of this class.
 *  - in the 'process' method of your analysis class, call the 'AndSelection::process' method
 *  - to get the result of the selection, ask the EventCalc container for passed(id), using the same id as used in the construction
 */
class AndSelection: public AnalysisModule{
public:
    explicit AndSelection(const identifier & selection_id, bool create_cutflow_histo = true);
    
    /// add the selection module, transferring memory ownership. If alternative_description is given, it is used in the histogram and output instead of module->description()
    void add(std::auto_ptr<SelectionModule> & module, const string & alternative_description = "");

    /// initialize histograms based on the current list of modules; call this exactly once per dataset.
    virtual void begin_dataset(Context & ctx);
    
    /// calculates the result of the selection and writes it to the event, using the selection_id passed in the constructor.
    virtual void process(EventCalc & event, Context & ctx);
    
private:
    
    // constructor variables:
    identifier selid;
    bool create_cutflow;
    
    // module info:
    boost::ptr_vector<SelectionModule> modules;
    std::vector<std::string> descriptions;
    
    TH1D * cutflow_raw, * cutflow_weighted; // owned by Context
};


/** some utility methods to convert string to various other objects; mainly for interpreting the configuration
 * file. 
 *
 * Will throw a runtime_error exception if a conversion is not possible.
 */
/** \brief Interpret the string as a boolean
 * 
 * Accepts "1"/"0"  "on"/"off"  "true"/"false"  "yes"/"no", ignoring case.
 */
bool string2bool(const string & s);

#endif

