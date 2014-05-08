#ifndef IDENTIFIER_HPP
#define IDENTIFIER_HPP

#include <string>
#include <inttypes.h>

/** \brief Program-unique identifiers
 *
 * Identifiers can be used as faster replacement for strings used to identify all sorts of things
 * such as histograms or event data in the input / output, etc. As with strings, the user is who
 * gives meaning to the identifier strings by convention.
 * 
 * identifier objects internally only consist of an id and are can be used as key in a map. They are thus
 * usually *much* faster than performing the same operation with strings.
 * 
 * Note that identifiers are (currently) *not* thread-safe.
 */
class identifier{
public:
    /// constant identifier constructed with an empty string
    static identifier empty;
    
    identifier(): id_(-1){}
    identifier(const char * c);
    identifier(const std::string & s);
    
    bool operator==(const identifier & other) const{
        return id_ == other.id_;
    }
    
    bool operator!=(const identifier & other) const{
        return id_ != other.id_;
    }
    
    bool operator<(const identifier & other) const{
        return id_ < other.id_;
    }

    std::string name() const;
    
    bool valid() const{
        return id_ >= 0;
    }
        
private:
   
   int64_t id_;
   
};


// short cut to define a variable of same variable name as identifier name.
// Use this either at global scope or in a function
#define ID(id) static identifier id(#id)


#endif
