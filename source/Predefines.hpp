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

class IDiskHump;
using IDiskHumpCPtr = std::shared_ptr<IDiskHump const>;

class MatterTranslation;
using MatterTranslationCPtr = std::shared_ptr<MatterTranslation const>;

class CartesianGrid;
using CartesianGridCPtr = std::shared_ptr<CartesianGrid const>;
using CartesianGridCRef = CartesianGridCPtr const &;

class Dust;
using DustCPtr = std::shared_ptr<Dust const>;
using DustCRef = DustCPtr const &;

#endif
