#ifndef PREDEFINES_HPP_
#define PREDEFINES_HPP_

#include <memory>

class Grid ;
class Pictures ;
class Observer;
class Photon ;
class Source ;
class Sources ;
class Directions ;  // directions grid

class Dust;
using DustCPtr = std::shared_ptr<Dust const>;
using DustCRef = DustCPtr const &;

#endif
