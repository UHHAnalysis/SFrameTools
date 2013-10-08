#include "SFrameTools/include/identifier.h"
#include <stdexcept>
#include <map>

namespace {

struct identifier_registry{
    
    identifier_registry(): next_id(0){}
    
    static identifier_registry & instance(){
        static identifier_registry i;
        return i;
    }
    
    size_t register_or_get(const std::string & name){
        std::map<std::string, size_t>::const_iterator it = name_to_id.find(name);
        if(it==name_to_id.end()){
            size_t id = next_id++;
            name_to_id[name] = id;
            id_to_name[id] = name;
            return id;
        }
        else{
            return it->second;
        }
    }
    
    std::map<std::string, size_t> name_to_id;
    std::map<size_t, std::string> id_to_name;
    size_t next_id;
};

}

identifier::identifier(const std::string & s): id_(identifier_registry::instance().register_or_get(s)){
}

identifier::identifier(const char * c): id_(identifier_registry::instance().register_or_get(c)){
}

std::string identifier::name() const{
    std::map<size_t, std::string> & id_to_name = identifier_registry::instance().id_to_name;
    std::map<size_t, std::string>::const_iterator it = id_to_name.find(id_);
    if(it==id_to_name.end()) return "<<invalid id>>";
    return it->second;
}

identifier identifier::empty("");
