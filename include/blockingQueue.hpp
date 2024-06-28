//
//  blockingQueue.hpp
//  MatrixCalculation
//
//  Created by Moreau Enzo on 10/04/2024.
//

#ifndef blockingQueue_hpp
#define blockingQueue_hpp
#include <boost/thread.hpp>
#include <deque>
#include <fstream>

#include <iostream>


template <typename T> class BlockingQueue {
  public:
    explicit BlockingQueue(size_t capacity) : _buffer(), _capacity(capacity) {
        assert(capacity>0);
    }

    void push(const T &elem) {
   
        _buffer.push_back(elem);
     ;
    }

    T pop() {
     

        if(!_buffer.empty())
        {
            T elem = _buffer.front();
            std::cout<<"Popped: ("<<std::get<0>(elem)<<" , "<<std::get<1>(elem)<<")\n";
            _buffer.pop_front();
            return elem;
        }
        else
        {
            return (T) std::tuple<int,int>(-1,1);
        }
     
       
        
    }

  private:
  

    std::deque<T> _buffer;
    size_t _capacity;
};


#endif
