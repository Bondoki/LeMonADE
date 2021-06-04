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

/*****************************************************************************/
/**
 * @file
 * @brief Tests for the class MoveConnectSc
 * @author Toni 
 * */
/*****************************************************************************/

#include <limits>

#include "gtest/gtest.h"

#include <LeMonADE/core/Molecules.h>
#include <LeMonADE/core/Ingredients.h>
#include <LeMonADE/feature/FeatureMoleculesIO.h>
#include <LeMonADE/feature/FeatureExcludedVolumeSc.h>
#include <LeMonADE/feature/FeatureConnectionSc.h>
#include <LeMonADE/feature/FeatureReactiveBonds.h>

#include <LeMonADE/updater/moves/MoveBreakReactive.h>
#include <LeMonADE/feature/FeatureReactiveBonds.h>


class TestMoveReactiveBreak: public ::testing::Test{
public:
  typedef LOKI_TYPELIST_3( FeatureMoleculesIO,FeatureReactiveBonds, FeatureExcludedVolumeSc<>) Features;
  typedef ConfigureSystem<VectorInt3,Features,3> Config;
  typedef Ingredients<Config> IngredientsType;

  IngredientsType ingredients;

  //redirect cout output
  virtual void SetUp(){
    originalBuffer=std::cout.rdbuf();
    std::cout.rdbuf(tempStream.rdbuf());
  };

  //restore original output
  virtual void TearDown(){
    std::cout.rdbuf(originalBuffer);
  };

private:
  std::streambuf* originalBuffer;
  std::ostringstream tempStream;
};


TEST_F(TestMoveReactiveBreak, checkAll)
{
  ingredients.setBoxX(16);
  ingredients.setBoxY(16);
  ingredients.setBoxZ(16);
  ingredients.setPeriodicX(true);
  ingredients.setPeriodicY(true);
  ingredients.setPeriodicZ(true);
  ingredients.modifyBondset().addBFMclassicBondset();
  ingredients.modifyMolecules().addMonomer(8,8,8);
  ingredients.modifyMolecules().addMonomer(10,8,8);
  ingredients.modifyMolecules().addMonomer(10,10,8);
  ingredients.modifyMolecules()[0].setReactive(true);
  ingredients.modifyMolecules()[1].setReactive(true);
  ingredients.modifyMolecules()[2].setReactive(false);
  ingredients.modifyMolecules()[0].setNumMaxLinks(1);
  ingredients.modifyMolecules()[1].setNumMaxLinks(2);
  ingredients.modifyMolecules()[2].setNumMaxLinks(2);
  ingredients.modifyMolecules().connect(0,1);
  ingredients.modifyMolecules().connect(1,2);
  EXPECT_NO_THROW(ingredients.synchronize());
  EXPECT_EQ(3,ingredients.getMolecules().size());
  EXPECT_TRUE(ingredients.getMolecules().areConnected(0,1));
  EXPECT_TRUE(ingredients.getMolecules().areConnected(1,2));
  
  MoveBreakReactive move;
  
  move.init(ingredients);
  EXPECT_EQ(0, move.getIndex());
  EXPECT_EQ(1, move.getPartner());
  EXPECT_FALSE( move.check(ingredients));
  
  move.init(ingredients,2);
  EXPECT_EQ(2, move.getIndex());
  EXPECT_EQ(std::numeric_limits<uint32_t>::max(), move.getPartner());
  EXPECT_FALSE( move.check(ingredients));
  
  move.init(ingredients,1);
  EXPECT_EQ(1, move.getIndex());
  EXPECT_EQ(0, move.getPartner());
  EXPECT_FALSE( move.check(ingredients));
  
  move.init(ingredients,2,1);
  EXPECT_EQ(2, move.getIndex());
  EXPECT_EQ(std::numeric_limits<uint32_t>::max(), move.getPartner());
  EXPECT_FALSE( move.check(ingredients));
  
  move.init(ingredients,0,1);
  EXPECT_EQ(0, move.getIndex());
  EXPECT_EQ(1, move.getPartner());
  EXPECT_FALSE( move.check(ingredients));
  move.apply(ingredients);
  
  EXPECT_FALSE(ingredients.getMolecules().areConnected(0,1));

  // check for connection across periodic boundary conditions
  ingredients.modifyMolecules().addMonomer(-3     ,1,1); // id 3 -> 3,0,0
  ingredients.modifyMolecules().addMonomer(0      ,1,1); // id 4 -> 2,1,1
  ingredients.modifyMolecules().addMonomer(2      ,3,2); // id 5 -> 3,1,0
  ingredients.modifyMolecules().addMonomer(5      ,4,2); // id 6 -> 2,0,0
  ingredients.modifyMolecules().addMonomer(7      ,4,2); // id 7 -> 2,1,0
  ingredients.modifyMolecules().addMonomer(9      ,5,2); // id 8 -> 2,1,1
  ingredients.modifyMolecules().addMonomer(11     ,6,3); // id 9
  
  
  ingredients.modifyMolecules()[3].setReactive(true);
  ingredients.modifyMolecules()[4].setReactive(true);
  ingredients.modifyMolecules()[5].setReactive(true);
  ingredients.modifyMolecules()[6].setReactive(true);
  ingredients.modifyMolecules()[7].setReactive(true);
  ingredients.modifyMolecules()[8].setReactive(true);
  ingredients.modifyMolecules()[9].setReactive(true);
  
  ingredients.modifyMolecules()[3].setNumMaxLinks(1);
  ingredients.modifyMolecules()[4].setNumMaxLinks(2);
  ingredients.modifyMolecules()[5].setNumMaxLinks(2);
  ingredients.modifyMolecules()[6].setNumMaxLinks(2);
  ingredients.modifyMolecules()[7].setNumMaxLinks(2);
  ingredients.modifyMolecules()[8].setNumMaxLinks(2);
  ingredients.modifyMolecules()[9].setNumMaxLinks(1);
  
  ingredients.modifyMolecules().connect(3,4);
  ingredients.modifyMolecules().connect(4,5);
  ingredients.modifyMolecules().connect(5,6);
  ingredients.modifyMolecules().connect(6,7);
  ingredients.modifyMolecules().connect(7,8);
  ingredients.modifyMolecules().connect(8,9);
  
  EXPECT_NO_THROW(ingredients.synchronize());
  EXPECT_EQ(10,ingredients.getMolecules().size());
  EXPECT_TRUE(ingredients.getMolecules().areConnected(3,4));
  EXPECT_TRUE(ingredients.getMolecules().areConnected(4,5));
  EXPECT_TRUE(ingredients.getMolecules().areConnected(5,6));
  EXPECT_TRUE(ingredients.getMolecules().areConnected(6,7));
  EXPECT_TRUE(ingredients.getMolecules().areConnected(7,8));
  EXPECT_TRUE(ingredients.getMolecules().areConnected(8,9));
 
  move.init(ingredients);
  EXPECT_TRUE((move.getIndex() >= 3) && (move.getIndex() <= 9)); // 0 and 1 are gone
  EXPECT_TRUE((move.getPartner() >= 3) && (move.getPartner() <= 9)); // 0 and 1 are gone
  
  // bonds between (3,4), (4,5) are allowed to break
  // bonds between (5,6), (6,7), (7,8), (8,9) are not allowed to break
  // bond with getIndex()==5 is indifferent
  EXPECT_TRUE((move.getIndex() <= 4) ? move.check(ingredients) : (move.getIndex() >= 6) ? !move.check(ingredients) : true );
  
}
