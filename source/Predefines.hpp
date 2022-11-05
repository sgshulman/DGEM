#ifndef PREDEFINES_HPP_
#define PREDEFINES_HPP_

#include <memory>

class Vector3d;

class Pictures ;
class Observer;
class Photon ;
class Source ;
class Directions ;  // directions grid

class Sources;
using SourcesPtr = std::shared_ptr<Sources>;
using SourcesRef = SourcesPtr const &;
using SourcesCPtr = std::shared_ptr<const Sources>;

class IMatter;
using IMatterCPtr = std::shared_ptr<IMatter const>;

class IDiskHump;
using IDiskHumpCPtr = std::shared_ptr<IDiskHump const>;

class MatterTranslation;
using MatterTranslationCPtr = std::shared_ptr<MatterTranslation const>;

class IGrid;
using IGridPtr = std::shared_ptr<IGrid>;
using IGridCPtr = std::shared_ptr<IGrid const>;
using IGridCRef = IGridCPtr const &;

class IDust;
using IDustCPtr = std::shared_ptr<IDust const>;
using IDustCRef = IDustCPtr const &;

class IRandomGenerator;

#endif
