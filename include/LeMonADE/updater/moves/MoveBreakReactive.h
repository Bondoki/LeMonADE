/*--------------------------------------------------------------------------------
    ooo      L   attice-based  |
  o\.|./o    e   xtensible     | LeMonADE: An Open Source Implementation of the
 o\.\|/./o   Mon te-Carlo      |           Bond-Fluctuation-Model for Polymers
oo---0---oo  A   lgorithm and  |
 o/./|\.\o   D   evelopment    | Copyright (C) 2013-2015,2021 by
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

#pragma once
#include <limits>
#include <vector>
#include <LeMonADE/updater/moves/MoveBreakBase.h>
#include <LeMonADE/utility/DistanceCalculation.h>

/*****************************************************************************/
/**
 * @file
 *
 * @class MoveBreakReactive
 *
 * @brief Standard local bfm-move on simple cubic lattice for the scBFM.
 *
 * @details The class is a specialization of MoveLocalBase using the (CRTP) to avoid virtual functions.
 **/
/*****************************************************************************/

class MoveBreakReactive:public MoveBreakBase<MoveBreakReactive>
{
public:
  MoveBreakReactive(){};

  // overload initialise function to be able to set the moves index and direction if neccessary
  template <class IngredientsType> void init(const IngredientsType& ing);
  template <class IngredientsType> void init(const IngredientsType& ing, uint32_t index);
  template <class IngredientsType> void init(const IngredientsType& ing, uint32_t index, uint32_t partner);

  template <class IngredientsType> bool check(IngredientsType& ing);
  template< class IngredientsType> void apply(IngredientsType& ing);

};
/////////////////////////////////////////////////////////////////////////////
/////////// implementation of the members ///////////////////////////////////

/*****************************************************************************/
/**
 * @brief Initialize the move.
 *
 * @details Resets the move probability to unity. Get a random bond vector from 
 * the feature FeatureReactiveBonds and set the corresponding IDs. 
 * @param ing A reference to the IngredientsType - mainly the system
 **/
template <class IngredientsType>
void MoveBreakReactive::init(const IngredientsType& ing)
{
  this->resetProbability();

  //draw index
  auto nReactiveBonds(ing.getNReactedBonds()); 
  if (nReactiveBonds==0) {
    this->setIndex(0);
    this->setPartner(std::numeric_limits<uint32_t>::max());    
  }else {
    auto index(this->randomNumbers.r250_rand32() % nReactiveBonds);
    auto it(ing.getBondedMonomers().cbegin());    
    std::advance( it, index);
    auto BondPair(it->first);
    this->setIndex(BondPair.first);
    this->setPartner(BondPair.second);
  }
}

/*****************************************************************************/
/**
 * @brief Initialize the move with a given monomer index.
 *
 * @details Resets the move probability to unity. Dice a new random direction.
 *
 * @param ing A reference to the IngredientsType - mainly the system
 * @param index index of the monomer to be connected
 **/
template <class IngredientsType>
void MoveBreakReactive::init(const IngredientsType& ing, uint32_t index)
{
  this->resetProbability();

  //set index
  if( (index >= 0) && (index < ing.getMolecules().size() ) )
    this->setIndex( index );
  else
    throw std::runtime_error("MoveBreakReactive::init(ing, index): index out of range !");


  std::vector<uint32_t> neighborIdx;
  for (auto i =0 ; i < ing.getMolecules().getNumLinks(this->getIndex() );i++ ) 
  { 
    auto neighbor(ing.getMolecules().getNeighborIdx(this->getIndex(), i));
    if (ing.getMolecules()[neighbor].isReactive()) 
      neighborIdx.push_back(neighbor);
  }
  auto nNeighbor(neighborIdx.size());
  if ( nNeighbor == 0  )
  {
    this->setPartner(std::numeric_limits<uint32_t>::max());
//     throw std::runtime_error("MoveBreakReactive::init(ing, index): index does not have a reactive bond partner!");
  }else 
    this->setPartner(neighborIdx[ this->randomNumbers.r250_rand32() % nNeighbor ] );
  if (!  ing.getMolecules()[index].isReactive())
    this->setPartner(std::numeric_limits<uint32_t>::max());
}

/*****************************************************************************/
/**
 * @brief Initialize the move with a given monomer index.
 *
 * @details Resets the move probability to unity. Dice a new random direction.
 *
 * @param ing A reference to the IngredientsType - mainly the system
 * @param index index of the monomer to be connected
 * @param partner index of the monomer to which index is connected
 **/
template <class IngredientsType>
void MoveBreakReactive::init(const IngredientsType& ing, uint32_t index, uint32_t partner)
{
  this->resetProbability();
  if ( ing.checkReactiveBondExists(index,partner) ) 
  {
    this->setIndex( index );
    this->setPartner(partner);
  }
  else 
  {
    this->setIndex( index );
    this->setPartner(std::numeric_limits<uint32_t>::max()) ;
//     throw std::runtime_error("MoveBreakReactive::init(ing, index, partner): index/partner out of range, not reactive or have no connection!");
  }

}
/*****************************************************************************/
/**
 * @brief Check if the move is accepted by the system.
 *
 * @details This function delegates the checking to the Feature.
 *
 * @param ing A reference to the IngredientsType - mainly the system
 * @return True if move is valid. False, otherwise.
 **/
template <class IngredientsType>
bool MoveBreakReactive::check(IngredientsType& ing)
{
    // check for valid id
  if (std::numeric_limits<uint32_t>::max() == this->getPartner() ) return false ; 
  
  // check if bond vector is out of P(2,2,1); P(3,0,0); P(3,1,0)
  //auto bv = ing.getMolecules()[this->getIndex()].getVector3D()-ing.getMolecules()[this->getPartner()].getVector3D();
  
  // fold the vector to if neccessary
  VectorInt3 bv = LemonadeDistCalcs::MinImageVector(ing.getMolecules()[this->getIndex()].getVector3D(), ing.getMolecules()[this->getPartner()].getVector3D(), ing);
  
  // reject short bond vectors P(2,0,0); P(2,1,0); P(2,1,1)
  if(bv*bv < 9)
      return false;
  
  //send the move to the Features to be checked
  return ing.checkMove(ing,*this);
}
/*****************************************************************************/
/**
 * @brief Apply the move to the system , e.g. add the displacement to Vertex (monomer) position.
 *
 * @details As first step: all Feature should apply the move using applyMove().\n
 * Second: Modify the positions etc. of the Vertex etc.
 *
 * @param ing A reference to the IngredientsType - mainly the system
 **/
template< class IngredientsType>
void MoveBreakReactive::apply(IngredientsType& ing)
{
	///@todo Think about the applying of move. Esp. make this independent of the order to avoid confusion!!
	///@todo check if it makes any difference in this case?! 
	//FIRST disconnect
	ing.modifyMolecules().disconnect(this->getIndex(),this->getPartner());	
	//THEN all features are applied 
  	ing.applyMove(ing,*this);
}

