#ifndef PREDEFINES_HPP_
#define PREDEFINES_HPP_

#include <memory>

class Pictures ;
class Observer;
class Photon ;
class Source ;
class Sources ;
class Directions ;  // directions grid

class FlaredDisk;
using FlaredDiskCPtr = std::shared_ptr<FlaredDisk const>;

class Grid;
using GridCPtr = std::shared_ptr<Grid const>;
using GridCRef = GridCPtr const &;

class Dust;
using DustCPtr = std::shared_ptr<Dust const>;
using DustCRef = DustCPtr const &;

#endif
