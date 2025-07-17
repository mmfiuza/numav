// Copyright (c) 2025 Matheus Machado Fiuza <matheusmachadofiuza@gmail.com>

#pragma once

#ifdef SAFE_PTR_DEBUG
    #define SAFE_PTR_DEBUG_BOOL 1
#else
    #define SAFE_PTR_DEBUG_BOOL 0
#endif

#ifdef SAFE_PTR_TEST
    #define SAFE_PTR_TEST_BOOL 1
#else
    #define SAFE_PTR_TEST_BOOL 0
#endif

#include <iostream>
#include <algorithm>
#include <stdexcept>
#if SAFE_PTR_DEBUG_BOOL
    #include <unordered_map>
    #include <mutex>
#endif

namespace fz {

#if SAFE_PTR_TEST_BOOL
    struct _SafePtrWarning {};
#endif

template<typename T>
class SafePtr
{
public:
    // constructor
    SafePtr() {
        #if SAFE_PTR_DEBUG_BOOL
            std::lock_guard<std::mutex> lock(_mtx);
            _memory_id = _get_new_memory_id();
            _ref_count[_memory_id] = 1;
            _is_deleted[_memory_id] = true;
        #endif
    }

    // constructor
    SafePtr(const size_t& size) {
        #if SAFE_PTR_DEBUG_BOOL
            std::lock_guard<std::mutex> lock(_mtx);
            _memory_id = _get_new_memory_id();
            _ref_count[_memory_id] = 1;
            _is_deleted[_memory_id] = false;
        #endif
        _begin = new T[size];
        _end = _begin + size;
    }

    // constructor
    SafePtr(const std::initializer_list<T>& il) {
        #if SAFE_PTR_DEBUG_BOOL
            std::lock_guard<std::mutex> lock(_mtx);
            _memory_id = _get_new_memory_id();
            _ref_count[_memory_id] = 1;
            _is_deleted[_memory_id] = false;
        #endif
        _begin = new T[il.size()];
        _end = _begin + il.size();
        std::copy(il.begin(), il.end(), this->_begin);
    }

    // destructor
    ~SafePtr() noexcept(!SAFE_PTR_DEBUG_BOOL) {
        #if SAFE_PTR_DEBUG_BOOL
            std::lock_guard<std::mutex> lock(_mtx);
            --_get_ref_count_nothrow();
            if (_get_ref_count_nothrow() == 0) {
                if(!_get_is_deleted_nothrow()) {
                    _warning("Memory was leaked.");
                }
                _ref_count.erase(_memory_id);
                _is_deleted.erase(_memory_id);
            }
        #endif
    }

    // copy constructor
    SafePtr(const SafePtr& other) {
        #if SAFE_PTR_DEBUG_BOOL
            std::lock_guard<std::mutex> lock(_mtx);
            other._check_for_usage_after_free();
            this->_memory_id = _get_new_memory_id();
            _ref_count[this->_memory_id] = 1;
            _is_deleted[this->_memory_id] = false;
        #endif
        this->_begin = new T[other.size()];
        this->_end = this->_begin + other.size();
        std::copy(other.begin(), other.end(), this->_begin);
    }
    
    // move constructor
    SafePtr(SafePtr&& other)
    noexcept(!SAFE_PTR_DEBUG_BOOL && !SAFE_PTR_TEST_BOOL) {
        #if SAFE_PTR_DEBUG_BOOL
            std::lock_guard<std::mutex> lock(_mtx);
            other._check_for_usage_after_free();
            this->_memory_id = other._memory_id;
            ++_get_ref_count_nothrow();
        #endif
        this->_begin = other._begin;
        this->_end = other._end;
    }

    // copy assignment operator
    SafePtr& operator=(const SafePtr& other) {
        #ifndef SAFE_PTR_DISABLE_SELF_ASSIGNING_CHECKING
            if (this != &other) {
        #endif
        #if SAFE_PTR_DEBUG_BOOL
            std::lock_guard<std::mutex> lock(_mtx);
            other._check_for_usage_after_free();
            --_get_ref_count_nothrow();
            if (_get_ref_count_nothrow()==0 && !_get_is_deleted_nothrow()) {
                _warning("Memory was leaked.");
            }
            this->_memory_id = _get_new_memory_id();
            _ref_count[this->_memory_id] = 1;
            _is_deleted[this->_memory_id] = false;
        #endif
        this->_begin = new T[other.size()];
        this->_end = this->_begin + other.size();
        std::copy(other.begin(), other.end(), this->_begin);
        #ifndef SAFE_PTR_DISABLE_SELF_ASSIGNING_CHECKING
            }
        #endif
        return *this;
    }

    // move assignment operator
    SafePtr& operator=(SafePtr&& other)
    noexcept(!SAFE_PTR_DEBUG_BOOL && !SAFE_PTR_TEST_BOOL) {
        #ifndef SAFE_PTR_DISABLE_SELF_ASSIGNING_CHECKING
            if (this != &other) {
        #endif
        #if SAFE_PTR_DEBUG_BOOL
            std::lock_guard<std::mutex> lock(_mtx);
            other._check_for_usage_after_free();
            --_get_ref_count_nothrow();
            if (_get_ref_count_nothrow()==0 && !_get_is_deleted_nothrow()) {
                _warning("Memory was leaked.");
            }
            this->_memory_id = other._memory_id;
            ++_get_ref_count_nothrow();
        #endif
        this->_begin = other._begin;
        this->_end = other._end;
        #ifndef SAFE_PTR_DISABLE_SELF_ASSIGNING_CHECKING
            }
        #endif
        return *this;
    }

    void free() const noexcept(!SAFE_PTR_DEBUG_BOOL && !SAFE_PTR_TEST_BOOL) {
        #if SAFE_PTR_DEBUG_BOOL
            std::lock_guard<std::mutex> lock(_mtx);
            if(_get_is_deleted_nothrow() == true) {
                throw std::logic_error(
                    "it was tried to free the same memory pointer twice"
                );
            }
            _get_is_deleted_nothrow() = true;
        #endif
        delete[] _begin;
    }

    size_t size() const noexcept(!SAFE_PTR_TEST_BOOL) {
        #if SAFE_PTR_DEBUG_BOOL
            _check_for_usage_after_free();
        #endif
        return _end - _begin;
    }

    const T* begin() const noexcept(!SAFE_PTR_TEST_BOOL) {
        #if SAFE_PTR_DEBUG_BOOL
            _check_for_usage_after_free();
        #endif
        return _begin;
    }

    T* begin() noexcept(!SAFE_PTR_TEST_BOOL) {
        return const_cast<T*>(
            const_cast<const SafePtr<T>&>(*this).begin()
        );
    }

    const T* cbegin() const noexcept(!SAFE_PTR_TEST_BOOL) {
        return begin();
    }

    const T* cbegin() noexcept(!SAFE_PTR_TEST_BOOL) {
        return const_cast<const SafePtr<T>&>(*this).cbegin();
    }

    const T* end() const noexcept(!SAFE_PTR_TEST_BOOL) {
        #if SAFE_PTR_DEBUG_BOOL
            _check_for_usage_after_free();
        #endif
        return _end;
    }

    T* end() noexcept(!SAFE_PTR_TEST_BOOL) {
        return const_cast<T*>(
            const_cast<const SafePtr<T>&>(*this).end()
        );
    }

    const T* cend() const noexcept(!SAFE_PTR_TEST_BOOL) {
        return end();
    }

    const T* cend() noexcept(!SAFE_PTR_TEST_BOOL) {
        return const_cast<const SafePtr<T>&>(*this).cend();
        
    }

    const T& operator[](const size_t& index)
    const noexcept(!SAFE_PTR_TEST_BOOL) {
        #if SAFE_PTR_DEBUG_BOOL
            _check_for_usage_after_free();
        #endif
        return *(_begin + index);
    }

    T& operator[](const size_t& index) noexcept(!SAFE_PTR_TEST_BOOL) {
        return const_cast<T&>(
            const_cast<const SafePtr<T>&>(*this)[index]
        );
    }

    const T& at(const size_t& index) const {
        #if SAFE_PTR_DEBUG_BOOL
            _check_for_usage_after_free();
        #endif
        if (index >= this->size()) {
            throw std::out_of_range(
                "tried to access SafePtr element out of range"
            );
        }
        return *(_begin + index);
    }

    T& at(const size_t& index) {
        return const_cast<T&>(
            const_cast<const SafePtr<T>&>(*this).at(index)
        );
    }

    bool empty() const noexcept(!SAFE_PTR_TEST_BOOL) {
        #if SAFE_PTR_DEBUG_BOOL
            _check_for_usage_after_free();
        #endif
        return this->size() == 0;
    }

    const T* data() const noexcept(!SAFE_PTR_TEST_BOOL) {
        #if SAFE_PTR_DEBUG_BOOL
            _check_for_usage_after_free();
        #endif
        return _begin;
    }

    T* data() noexcept(!SAFE_PTR_TEST_BOOL) {
        return const_cast<T*>(
            const_cast<const SafePtr<T>&>(*this).data()
        );
    }

    const T& front() const noexcept(!SAFE_PTR_TEST_BOOL) {
        #if SAFE_PTR_DEBUG_BOOL
            _check_for_usage_after_free();
        #endif
        return *(this->_begin);
    }

    T& front() noexcept(!SAFE_PTR_TEST_BOOL) {
        return const_cast<T&> (
            const_cast<const SafePtr<T>&>(*this).front()
        );
    }

    const T& back() const noexcept(!SAFE_PTR_TEST_BOOL) {
        #if SAFE_PTR_DEBUG_BOOL
            _check_for_usage_after_free();
        #endif
        return *(this->_end-1);
    }

    T& back() noexcept(!SAFE_PTR_TEST_BOOL) {
        return const_cast<T&>(
            const_cast<const SafePtr<T>&>(*this).back()
        );
    }

    void fill(const T& value) {
        #if SAFE_PTR_DEBUG_BOOL
            _check_for_usage_after_free();
        #endif
        for (auto& p : *this) {
            p = value;
        }
    }

    void print_all(const char* const variable_name = "SafePtr::print_all") 
    const {
        #if SAFE_PTR_DEBUG_BOOL
            _check_for_usage_after_free();
        #endif
        std::cout << variable_name << ": {\n";
        if (!empty()) {
            size_t idx = 0;
            for (const T* it=cbegin(); it!=cend()-1; ++it) {
                std::cout << "    " << idx << ": " << *it << ",\n";
                ++idx;
            }
            std::cout << "    " << idx << ": " << back() << "\n";
        }
        std::cout << "}\n";
    }

    void print(const char* const variable_name = "SafePtr::print") const {
        #if SAFE_PTR_DEBUG_BOOL
            _check_for_usage_after_free();
        #endif

        constexpr size_t MAX_ELEMENTS_TO_PRINT = 26;
        if (size() <= MAX_ELEMENTS_TO_PRINT) {
            print_all(variable_name);
            return;
        }
        
        constexpr size_t ELEMENTS_BEFORE_DOTS = 13;
        constexpr size_t ELEMENTS_AFTER_DOTS = 12;
        std::cout << variable_name << ": {\n";
        size_t idx = 0;
        for (const T* it=cbegin(); it!=cbegin()+ELEMENTS_BEFORE_DOTS; ++it) {
            std::cout << "    " << idx << ": " << *it << ",\n";
            ++idx;
        }
        idx = size() - ELEMENTS_AFTER_DOTS;
        std::cout << "    ...\n";
        for (const T* it=cend()-ELEMENTS_AFTER_DOTS; it!=cend()-1; ++it) {
            std::cout << "    " << idx << ": " << *it << ",\n";
            ++idx;
        }
        std::cout << "    " << idx << ": " << back() << "\n";
        std::cout << "}\n";
    }

private:
    T* _begin; // points to the first element
    T* _end;   // points to the byte after the last byte of the element

    #if SAFE_PTR_DEBUG_BOOL
        size_t _memory_id;
        static size_t _next_available_memory_id;
        static std::unordered_map<size_t,size_t> _ref_count;
        static std::unordered_map<size_t,bool> _is_deleted;
        static bool _id_overflow_occurred;
        static std::mutex _mtx;

        static size_t& _get_new_memory_id() noexcept {
            ++_next_available_memory_id;
            if (_next_available_memory_id == 0) {
                _id_overflow_occurred = true;
            }
            if (!_id_overflow_occurred) {
                return _next_available_memory_id;
            }
            while (true) { // handle overflow
                if (_ref_count.count(_next_available_memory_id) == 0) {
                    return _next_available_memory_id;
                }
                ++_next_available_memory_id;
            }
        }

        void _check_for_usage_after_free()
        const noexcept(!SAFE_PTR_TEST_BOOL) {
            if (_get_is_deleted_nothrow() == true) {
                _warning("Tried to access data after free() was called.");
            }
        }

        size_t& _get_ref_count_nothrow() const noexcept {
            typename std::unordered_map<size_t,size_t>::iterator it = 
                _ref_count.find(_memory_id);
            return it->second;
        }

        bool& _get_is_deleted_nothrow() const noexcept {
            typename std::unordered_map<size_t,bool>::iterator it =
                _is_deleted.find(_memory_id);
            return it->second;
        }

        void _warning(const char* const msg) const {
            #if SAFE_PTR_TEST_BOOL
                throw _SafePtrWarning();
            #endif
            std::cerr << 
                "\033[33m" << "SafePtr warning: " << "\033[0m" << msg << "\n";
        }
    #endif
};

#if SAFE_PTR_DEBUG_BOOL
    template<typename T> size_t SafePtr<T>::_next_available_memory_id = 0;

    template<typename T>
    std::unordered_map<size_t,size_t> SafePtr<T>::_ref_count;
    template<typename T>
    std::unordered_map<size_t,bool> SafePtr<T>::_is_deleted;
    template<typename T> bool SafePtr<T>::_id_overflow_occurred = false;
    
    template<typename T> std::mutex SafePtr<T>::_mtx;
#endif

} // namespace fz
