#if !defined(REGISTRY_HPP) && !defined(__CINT__) // completely hide this from CINT, which does not like boost constructs
#define REGISTRY_HPP

#include <stdint.h>
#include <memory>
#include <stdexcept>

#include <boost/ptr_container/ptr_map.hpp>
#include <boost/static_assert.hpp>
#include <boost/type_traits.hpp>

template<typename base_type_>
class Registry{
public:
    static Registry & instance(){
        static Registry i;
        return i;
    }
    
    typedef base_type_ base_type;
        
    std::auto_ptr<base_type> build(const std::string & key) const;
    
    template<typename T>
    int register_(const std::string & key);
    
private:
    Registry(){}
    
    struct Factory{
        virtual std::auto_ptr<base_type> operator()() const = 0;
        virtual ~Factory(){}
    };
    
    template<typename T>
    struct FactoryDefault: public Factory{
        virtual std::auto_ptr<base_type> operator()() const{
            return std::auto_ptr<base_type>(new T());
        }
        virtual ~FactoryDefault(){}
    };
    
    boost::ptr_map<std::string, Factory> id_to_factory;
};


template<typename base_type>
template<typename T>
int Registry<base_type>::register_(const std::string & key){
    const typename boost::is_base_of<base_type, T>::value_type is_base = boost::is_base_of<base_type, T>::value;
    BOOST_STATIC_ASSERT(is_base);
    typename boost::ptr_map<std::string, Factory>::const_iterator it = id_to_factory.find(key);
    if(it != id_to_factory.end()){
        throw std::runtime_error("Registry: tried to register type with same name twice (name: '" + key + "')");
    }
    std::string keytmp = key;
    id_to_factory.insert(keytmp, new FactoryDefault<T>());
    return 0;
}

template<typename base_type>
std::auto_ptr<base_type> Registry<base_type>::build(const std::string & key) const {
    typename boost::ptr_map<std::string, Factory>::const_iterator it = id_to_factory.find(key);
    if(it==id_to_factory.end()) throw std::runtime_error("Registry: did not find registered type of name '" + key + "'");
    try{
        return it->operator()();
    }
    catch(std::exception & ex){
        throw std::runtime_error("Registry: exception while trying to build type '" + key + "': " + ex.what());
    }
}

// register type T with base tye BT and given name
#define REGISTER(T, BT, name) namespace { int dummy##T = ::Registry<BT>::instance().register_<T>(#name); }

#endif
