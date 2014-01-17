#include "RVersion.h"

#if ROOT_VERSION_CODE == 336391
// see https://sft.its.cern.ch/jira/browse/ROOT-5959
#warning "applying fix for ReadCollectionBool"

#include "TBuffer.h"
#include "TStreamerInfoActions.h"

using namespace TStreamerInfoActions;


namespace TStreamerInfoActions {
    
// see TStreamerInfoActions.cxx
class TConfigSTL : public TConfiguration {
      // Configuration object for the kSTL case
   private:
      void Init();
   
   public:
      TClass          *fOldClass;   // Class of the content on file
      TClass          *fNewClass;   // Class of the content in memory.
      TMemberStreamer *fStreamer;
      const char      *fTypeName;   // Type name of the member as typed by ther user.
      Bool_t          fIsSTLBase;  // aElement->IsBase() && aElement->IsA()!=TStreamerBase::Class()

      TVirtualCollectionProxy::CreateIterators_t    fCreateIterators;
      TVirtualCollectionProxy::CopyIterator_t       fCopyIterator;
      TVirtualCollectionProxy::DeleteIterator_t     fDeleteIterator;
      TVirtualCollectionProxy::DeleteTwoIterators_t fDeleteTwoIterators;

      TConfigSTL(TVirtualStreamerInfo *info, UInt_t id, Int_t offset, UInt_t length, TClass *oldClass, const char *type_name, Bool_t isbase);

      TConfigSTL(TVirtualStreamerInfo *info, UInt_t id, Int_t offset, UInt_t length, TClass *oldClass, TClass *newClass, const char *type_name, Bool_t isbase);

      TConfigSTL(TVirtualStreamerInfo *info, UInt_t id, Int_t offset, UInt_t length, TClass *oldClass, TMemberStreamer* streamer, const char *type_name, Bool_t isbase);

      TConfigSTL(TVirtualStreamerInfo *info, UInt_t id, Int_t offset, UInt_t length, TClass *oldClass, TClass *newClass, TMemberStreamer* streamer, const char *type_name, Bool_t isbase);

      virtual TConfiguration *Copy();
   };
}


extern "C" int _ZN20TStreamerInfoActions12VectorLooper18ReadCollectionBoolER7TBufferPvPKNS_14TConfigurationE(TBuffer &buf, void *addr, const TConfiguration *conf)
      {
         // Collection of numbers.  Memberwise or not, it is all the same.

         TConfigSTL *config = (TConfigSTL*)conf;
         UInt_t start, count;
         /* Version_t vers = */ buf.ReadVersion(&start, &count, config->fOldClass);

         std::vector<bool> *const vec = (std::vector<bool>*)(((char*)addr)+config->fOffset);
         Int_t nvalues;
         buf.ReadInt(nvalues);
         vec->resize(nvalues);

         bool *items = new bool[nvalues];
         buf.ReadFastArray(items, nvalues);
         for(Int_t i = 0 ; i < nvalues; ++i) {
            (*vec)[i] = items[i];
         }
         delete[] items;

         // We could avoid the call to ReadFastArray, and we could
         // the following, however this breaks TBufferXML ...
         // for(Int_t i = 0 ; i < nvalues; ++i) {
         //    bool tmp; buf >> tmp;
         //    (*vec)[i] = tmp;
         // }

         buf.CheckByteCount(start,count,config->fTypeName);
         return 0;
}
#else
#warning "NOT applying fix for ReadCollectionBool: not sure whether it's Ok for this version of ROOT (please investigate!)"
#endif
