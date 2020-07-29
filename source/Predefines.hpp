#ifndef PREDEFINES_HPP_
#define PREDEFINES_HPP_

#include <memory>

class Pictures ;
class Observer;
class Photon ;
class Source ;
class Directions ;  // directions grid
class Random;

class Sources;
using SourcesPtr = std::shared_ptr<Sources>;
using SourcesRef = SourcesPtr const &;

class IMatter;
using IMatterCPtr = std::shared_ptr<IMatter const>;

class MatterTranslation;
using MatterTranslationCPtr = std::shared_ptr<MatterTranslation const>;

class Grid;
using GridCPtr = std::shared_ptr<Grid const>;
using GridCRef = GridCPtr const &;

class Dust;
using DustCPtr = std::shared_ptr<Dust const>;
using DustCRef = DustCPtr const &;

#endif
