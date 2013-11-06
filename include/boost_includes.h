// This file includes the boost dependencies of SFrame that do not play well with CINT./
// If run with CINT, some dummy definitions are provided.
#ifndef SFRAMETOOLS_BOOST_INCLUDES_H
#define SFRAMETOOLS_BOOST_INCLUDES_H


#ifndef __CINT__

#include <boost/shared_array.hpp>
#include <boost/static_assert.hpp>
#include <boost/type_traits.hpp>
#include <boost/ptr_container/ptr_vector.hpp>
#include <boost/lexical_cast.hpp>

#else

// some summy definitions to make cint happy:
namespace boost{
    template<typename T> struct shared_array{};
    template<typename T> struct ptr_vector{};
}

#define BOOST_STATIC_ASSERT_MSG(a,b) {}

#endif

#endif
