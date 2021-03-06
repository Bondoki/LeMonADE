/*--------------------------------------------------------------------------------
    ooo      L   attice-based  |
  o\.|./o    e   xtensible     | LeMonADE: An Open Source Implementation of the
 o\.\|/./o   Mon te-Carlo      |           Bond-Fluctuation-Model for Polymers
oo---0---oo  A   lgorithm and  |
 o/./|\.\o   D   evelopment    | Copyright (C) 2013-2015 by
  o/.|.\o    E   nvironment    | LeMonADE Principal Developers (see AUTHORS)
    ooo                        |
----------------------------------------------------------------------------------

This file is part of LeMonADE.

LeMonADE is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

LeMonADE is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with LeMonADE.  If not, see <http://www.gnu.org/licenses/>.

--------------------------------------------------------------------------------*/

#ifndef LEMONADE_UPDATER_MOVES_MOVELABELBASE_H
#define LEMONADE_UPDATER_MOVES_MOVELABELBASE_H

#include <LeMonADE/updater/moves/MoveBase.h>
#include <LeMonADE/utility/Vector3D.h>
#include <LeMonADE/utility/RandomNumberGenerators.h>

/*****************************************************************************/
/**
 * @file
 *
 * @class MoveLabelBase
 *
 * @brief Base class for all simple local bfm-moves. Specialized versions is MoveLabelAlongChain.
 *
 * @details This class provides a base for local bfm-moves. The implementation
 * details (i.e. MoveLabelAlongChain,  optimized versions,etc.) are given by the template parameter.
 * This way, using the "Curiously Recurring Template Pattern (CRTP) to avoid virtual functions
 * but still providing their functionality. A specialized local
 * move type derived from this class must then be constructed in the following
 * way: "class MySpecialMove:public MoveLabelBase<MySpecialMove>".
 * It must implement the functions init,check and apply. Calls to the corresponding
 * functions in this base class are then redirected to the specialized implementation.
 * See for an example the class MoveLabelAlongChain.
 * Since this class simply serves as a common base, the functions apply(), check(), and init() don't do
 * anything particular.
 *
 * @tparam <SpecializedMove> name of the specialized move.
 *
 **/
/*****************************************************************************/
template <class SpecializedMove>
class MoveLabelBase:public MoveBase
{
 public:
	//! Returns the index of the Vertex (monomer) in the graph which should be moved
	uint32_t getIndex() const {return index;}
	
	//! Returns the direction of the Vertex (monomer) either right or left 
	const int32_t& getDir() const { return dir;}
	
	//!
	uint32_t getConnectedLabel() const {return connectedLabel;}

	//here come the functions that are implemented by the specialization
	template <class IngredientsType> void init(const IngredientsType& ingredients);
	template <class IngredientsType> void init(const IngredientsType& ing, uint32_t index);
	template <class IngredientsType> void init(const IngredientsType& ing, uint32_t index, int32_t dir);
	template <class IngredientsType> void check(const IngredientsType& ingredients);
	template <class IngredientsType> void apply(IngredientsType& ingredients);

 protected:
	/**
	 * @brief Set the index of the Vertex (monomer) in the graph which should be moved
	 * @param i the index in the graph
	 */
	void setIndex(uint32_t i) {index=i;}

	/**
	 * @brief Set the move direction of the Vertex (monomer) which should be moved
	 * @param dir The displacement in direction on the Cartesian-space.
	 */
	void setDir(const int32_t& dir_) {dir=dir_;}
	
	/**
	 * 
	 */
	void setConnectedLabel(uint32_t connectedLabel_) { connectedLabel=connectedLabel_;}


	//! Random Number Generator (RNG)
	RandomNumberGenerators randomNumbers;


 private:

	//! Index of the Vertex (monomer) in the graph which should be moved
	uint32_t index;

	//! Direction for the move
	int32_t dir;
	
	//!starting from 1 
	uint32_t connectedLabel;
};



////////////////////////////////////////////////////////////////////////////////
// implementation of the members
////////////////////////////////////////////////////////////////////////////////
/*****************************************************************************/
/**
 * @brief Initialize the move. Done by the SpecializedMove.
 *
 * @details Here, this is only redirected to the implementation
 * given in the template parameter
 *
 * @tparam <SpecializedMove> name of the specialized move.
 **/
/*****************************************************************************/
template <class SpecializedMove>
template <class IngredientsType>
void MoveLabelBase<SpecializedMove>::init(const IngredientsType& ingredients)
{
  static_cast<SpecializedMove*>(this)->init(ingredients);
}

/*****************************************************************************/
/**
 * @brief Initialize the move with a given monomer index. Done by the SpecializedMove.
 *
 * @details Here, this is only redirected to the implementation
 * given in the template parameter
 *
 * @tparam <SpecializedMove> name of the specialized move.
 **/
/*****************************************************************************/
template <class SpecializedMove>
template <class IngredientsType>
void MoveLabelBase<SpecializedMove>::init(const IngredientsType& ingredients, uint32_t index)
{
  static_cast<SpecializedMove*>(this)->init(ingredients, index);
}

/*****************************************************************************/
/**
 * @brief Initialize the move with a given monomer index and direction. Done by the SpecializedMove.
 *
 * @details Here, this is only redirected to the implementation
 * given in the template parameter
 *
 * @tparam <SpecializedMove> name of the specialized move.
 **/
/*****************************************************************************/
template <class SpecializedMove>
template <class IngredientsType>
void MoveLabelBase<SpecializedMove>::init(const IngredientsType& ingredients, uint32_t index, int32_t dir)
{
  static_cast<SpecializedMove*>(this)->init(ingredients, index, dir);
}
/*****************************************************************************/
/**
 * @brief Check if the move is accepted by the system. Done by the SpecializedMove.
 *
 * @details Here, this is only redirected to the implementation
 * given in the template parameter.
 *
 * @tparam <SpecializedMove> name of the specialized move.
 **/
/*****************************************************************************/
template <class SpecializedMove>
template <class IngredientsType>
void MoveLabelBase<SpecializedMove>::check(const IngredientsType& ingredients)
{
  static_cast<SpecializedMove*>(this)->check(ingredients);
}

/*****************************************************************************/
/**
 * @brief Apply the move to the system.
 *
 * @details Here, this is only redirected to the implementation
 * given in the template parameter
 *
 * @tparam <SpecializedMove> name of the specialized move.
 * */
/*****************************************************************************/
template <class SpecializedMove>
template <class IngredientsType>
void MoveLabelBase<SpecializedMove>::apply(IngredientsType& ingredients)
{
  static_cast<SpecializedMove*>(this)->apply(ingredients);
}

#endif
