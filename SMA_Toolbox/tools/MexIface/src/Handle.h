/** @file Handle.h
 * @author Mark J. Olah (mjo\@cs.unm.edu)
 * @date 09-2014
 * @brief Helper class and templated functions to represent and manipulate Handles to C++ objects as Matlab mxArrays
 *
 * @copyright Mark J. Olah and The Regents of the University of New Mexico (2014).
 *            This code is free for non-commercial use and modification, provided
 *            this copyright notice remains unmodified and attached to the code
 */

#ifndef _HANDLE_H
#define _HANDLE_H
#include "mex.h"
#include <cstdint>
#include <string>
#include <typeinfo>


/** @brief A class to represent and manipulate handles to C++ classes that can be wrapped as
 * Matlab arrays, allowing C++ objects to persist between Mex calls.
 *
 * This allows Matlab to hold a handle to a C++ object allocated during one Mex call, but used
 * during subsequent calls.
 */
template<class T> class Handle
{
public:
    /**
     * @brief The C++ datatype of corresponding to the array type a Handle pointer will be stored
     *  in for the Matlab side of things
     */
    typedef uint64_t HandlePtrT;

    Handle(T *obj);
    ~Handle();
    
    bool is_valid() const;
    T* object() const;

    static mxArray* makeHandle(T *obj);
    static Handle<T>* getHandle(const mxArray *arr);
    static inline T* getObject(const mxArray *in);
    static void destroyObject(const mxArray *in);

private:
    static const uint32_t class_handle_signature=0xFF00F0A5; /**< Matlab signature for class handle type */
    uint32_t signature; /**< The stored signature (should match class_handle_signature */
    std::string name; /**< The stored name of the class, should match with templated type T */
    T *obj; /**< The actual pointer to the object stored */
};

/* Templated Member Functions */


/**
 * @brief Make a new handle object to hold a pointer to given object of class T
 * @param obj The object to wrap in a handle and make persistent
 *
 * Note: we assume ownership of obj and will free it with a call to delete when the
 * handle itself is deleted.  This assumes obj was created with a call to new.
 *
 */
template<class T>
Handle<T>::Handle(T *obj)
    : signature(class_handle_signature),
      name(typeid(T).name()),
      obj(obj)
{
}


/**
 * @brief Delete the object which was assumed to have been created with new.
 */
template<class T>
Handle<T>::~Handle()
{
    signature=0;
    delete obj;
}


/**
 * @brief Check that this is a valid handle to a valid C++ object
 * @returns True if valid
 */
template<class T>
bool Handle<T>::is_valid() const
{
    bool sig_ok= (signature == class_handle_signature);
    bool name_ok= (name == typeid(T).name());
    return sig_ok && name_ok;
}


/**
 * @brief Retrieve a pointer to the object stored by this handle
 * @returns The pointer to the object
 */
template<class T>
inline
T* Handle<T>::object() const
{
    return obj;
}

/* Templated Static Member Functions */

/**
 * @brief Given a pointer to a C++ object, make a new Handle object and save that as
 *    a uint64_t in a Matlab mxArray object.
 * @param obj The object to wrap.
 * @returns A mxArray that contains the handle as a numeric scalar uint64_t
 */
template<class T>
mxArray* Handle<T>::makeHandle(T *obj)
{
    mexLock(); /* Increment the lock count */
    auto arr = mxCreateNumericMatrix(1, 1, mxUINT64_CLASS, mxREAL); //Make a new numeric array to hold the handle address
    auto handle_data = static_cast<HandlePtrT*>(mxGetData(arr)); //Get pointer to arrays internal storage
    *handle_data = reinterpret_cast<HandlePtrT>(new Handle<T>(obj)); //Save the pointer as a uint64_t
    return arr;
}

/**
 * @brief Given a Matlab mxArray object pointer to data that represents a handle, return a Handle object pointer.
 * @param arr A Matlab mxArray where the handle is stored as a uint64_t scalar.
 * @returns A pointer to the Handle object
 */
template<class T>
Handle<T>* Handle<T>::getHandle(const mxArray *arr)
{
    if (mxGetNumberOfElements(arr) != 1) mexErrMsgIdAndTxt("ClassHandle:covertMat2HandlePtr","Non scalar.");
    if (mxGetClassID(arr) != mxUINT64_CLASS) mexErrMsgIdAndTxt("ClassHandle:covertMat2HandlePtr","not UINT64");
    if (mxIsComplex(arr)) mexErrMsgIdAndTxt("ClassHandle:covertMat2HandlePtr","is complex");
    auto handle_data=static_cast<HandlePtrT*>(mxGetData(arr));
    auto handle = reinterpret_cast<Handle<T>*>( *handle_data );
    if (!handle->is_valid()) mexErrMsgIdAndTxt("ClassHandle:covertMat2HandlePtr","Handle not valid.");
    return handle;
}

/**
 * @brief Given a matlab mxArray object pointer to data that represents a handle, retrieve the object pointer for the C++ object.
 * @param arr A Matlab mxArray where the handle is stored as a uint64_t scalar.
 * @returns A pointer to the actual object that was wrapped in the Handle object that was itself stored in the array numerically
 */
template<class T>
T* Handle<T>::getObject(const mxArray *arr)
{
    return getHandle(arr)->object();
}

/**
 * @brief Given a matlab mxArray object pointer to data that represents a handle, delete the handle which also implies deleting
 *    object itself.
 * @param arr The Matlab mxArray that contains the numerically encoded pointer to the handle we wish to destroy
 *
 * The Handle object as well as the object that wrapped it are assumed to have been created with a call to new, and we are assuming
 * that the Handle object now "owns" the memory of the wrapped object and thus is responsible for freeing it.
 *
 * This also decrements the mexLock count, as we have freed one of the persistent object we previously created.
 */
template<class T>
void Handle<T>::destroyObject(const mxArray *arr)
{
    delete getHandle(arr); //This destroys the handle object, and the Handle's destructor will destroy the actual wrapped object.
    mexUnlock(); /* Decrement the lock count */
}

#endif // _HANDLE_H
