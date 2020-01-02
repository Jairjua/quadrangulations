/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   test_help_functions.cpp
 * Author: eu
 *
 * Created on December 5, 2019, 2:32 PM
 */

#include <stdlib.h>
#include <iostream>
#include "help_functions_for_tests.h"
#include "vcg/complex/algorithms/bitquad_creation.h"
/*
 * Simple C++ Test Suite
 */ //pq cout não funciona?
using std::string;

void testPairToBitQuad() {
    MeshType m;
    CreateTestMesh(m,1,true,true,false);
    FaceType* fi = &*m.face.begin();
    
    int whichEdge = 0 ; //borda
    PairToBitQuad<false>(fi, whichEdge);
    if (fi->IsF(whichEdge)) {
        std::cout << "%TEST_FAILED% time=0 testname=testPairToBitQuad (test_bq_creation) message=emparelhou borda=" << whichEdge<< std::endl;
    }
    
    whichEdge = 1 ; //cluster
    PairToBitQuad<false>(fi, whichEdge);
    if (!fi->IsF(whichEdge)) {
        std::cout << "%TEST_FAILED% time=0 testname=testPairToBitQuad (test_bq_creation) message= não emparelhou =" << whichEdge<< std::endl;
    }
    
    fi->ClearF(whichEdge);
    PairToBitQuad<false>(fi, whichEdge);
    if (fi->IsF(whichEdge)) {
        std::cout << "%TEST_FAILED% time=0 testname=testPairToBitQuad (test_bq_creation) message= Não detectou faux edge =" << whichEdge<< std::endl;
    }
}

void testPair(){
  MeshType m;
  CreateTestMesh(m,1,false,true,false);
  if (Pair(&*m.face.begin(),0)) {
        std::cout << "%TEST_FAILED% time=0 testname=testPair (test_help_functions) message=não detectou a borda"<< std::endl;
    }
  
  Pair(&*m.face.begin(),1);
  if (!(m.face.begin()->IsF(1) && m.face.begin()->FFp(1)->IsF(m.face.begin()->FFi(1)))) {
        std::cout << "%TEST_FAILED% time=0 testname=testPair (test_help_functions) message=não fez faux"<< std::endl;
    }
}

void testCountTriAndQuad(){
  MeshType m;
  CreateTestMesh(m,0,false,true,false);
  int nquads=0;
  int ntris=0;
  int n = CountTriAndQuad(m,nquads,ntris);
  if (!(n==0 && nquads==0 && ntris==0)) {
        std::cout << "%TEST_FAILED% time=0 testname=testCountTriAndQuad (test_help_functions) message=m0 . nquads=" <<nquads<<" ntris="<<ntris<<" n="<<n << std::endl;
    }
  
  CreateTestMesh(m,1,false,true,false);
  nquads=0;
  ntris=0;
  n = CountTriAndQuad(m,nquads,ntris);
  if (!(n==0 && nquads==0 && ntris==2)) {
        std::cout << "%TEST_FAILED% time=0 testname=testCountTriAndQuad (test_help_functions) message=m1 . nquads=" <<nquads<<" ntris="<<ntris<<" n="<<n << std::endl;
    }
  
  CreateTestMesh(m,3,false,true,false);
  nquads=0;
  ntris=0;
  n = CountTriAndQuad(m,nquads,ntris);
  if (!(n==0 && nquads==0 && ntris==6)) {
        std::cout << "%TEST_FAILED% time=0 testname=testCountTriAndQuad (test_help_functions) message=m3.1 . nquads=" <<nquads<<" ntris="<<ntris<<" n="<<n << std::endl;
    }
  
  vcg::tri::BitQuadCreation<MeshType>::MakeDominant(m,2);
  n = CountTriAndQuad(m,nquads,ntris);
  if (!(n==2 && nquads==2 && ntris==2)) {
        std::cout << "%TEST_FAILED% time=0 testname=testCountTriAndQuad (test_help_functions) message=m3.2 . nquads=" <<nquads<<" ntris="<<ntris<<" n="<<n << std::endl;
    }
  
  m.face.begin()->ClearAllF();
  nquads=0;
  ntris=0;
  n = CountTriAndQuad(m,nquads,ntris);
  if (!(n==-1 && nquads==-1 && ntris==0)) {
        std::cout << "%TEST_FAILED% time=0 testname=testCountTriAndQuad (test_help_functions) message=m3.3 . nquads=" <<nquads<<" ntris="<<ntris<<" n="<<n << std::endl;
    }
  
  m.face.begin()->SetF(1);
  m.face.begin()->FFp(1)->SetF( m.face.begin()->FFi(1) );
  m.face.begin()->SetF(2);
  m.face.begin()->FFp(2)->SetF( m.face.begin()->FFi(2) );
  nquads=0;
  ntris=0;
  n = CountTriAndQuad(m,nquads,ntris);
  if (!(n==-2 && nquads==0 && ntris==0)) {
        std::cout << "%TEST_FAILED% time=0 testname=testCountTriAndQuad (test_help_functions) message=m3.4 . nquads=" <<nquads<<" ntris="<<ntris<<" n="<<n << std::endl;
    }
  
  
}
void testCreateTestMesh() {
    MeshType m;

    int i=0;
    bool showCoord=false;
    bool updateFFAdj=false;
    bool updateNormalFace=false;
    if (string("No mesh created") != CreateTestMesh(m, i, showCoord, updateFFAdj, updateNormalFace)) {
        std::cout << "%TEST_FAILED% time=0 testname=testCreateTestMesh (test_help_functions) message=texto da malha 0 errado" << std::endl;
    }
    
    i=1;
    showCoord=true;
    CreateTestMesh(m, i, showCoord, updateFFAdj, updateNormalFace);
    if ( m.FN()!= 2) {
        std::cout << "%TEST_FAILED% time=0 testname=testCreateTestMesh (test_help_functions) message= nº faces diferente de 2. FN=" << m.FN() << std::endl;
    }
    if ( m.face.begin()->P(1)[0] != 2) {
        std::cout << "%TEST_FAILED% time=0 testname=testCreateTestMesh (test_help_functions) message= coordenada errada =" << m.face.begin()->P(1)[0] << std::endl;
    }
    if ( m.face.begin()->FFp(0) != 0) {
        std::cout << "%TEST_FAILED% time=0 testname=testCreateTestMesh (test_help_functions) message= ponteiro errado" << m.face.begin()->FFp(0) << std::endl;
    }
    
    i=1;
    showCoord=false;
    updateFFAdj=true;
    CreateTestMesh(m, i, showCoord, updateFFAdj, updateNormalFace);
    if ( m.FN()!= 2) {
        std::cout << "%TEST_FAILED% time=0 testname=testCreateTestMesh (test_help_functions) message= nº faces diferente de 2. FN=" << m.FN() << std::endl;
    }
    if ( !(m.face.begin()->P(1)[0] == 2)) {
        std::cout << "%TEST_FAILED% time=0 testname=testCreateTestMesh (test_help_functions) message= coordenada errada =" << m.face.begin()->P(1)[0] << std::endl;
    }
    if ( !(m.face.begin()->FFp(0) == &*m.face.begin())) {
        std::cout << "%TEST_FAILED% time=0 testname=testCreateTestMesh (test_help_functions) message= ponteiro errado" << m.face.begin()->FFp(0) << std::endl;
    }
    if ( !(m.face.begin()->FFp(1) != &*m.face.begin())) {
        std::cout << "%TEST_FAILED% time=0 testname=testCreateTestMesh (test_help_functions) message= ponteiro errado" << m.face.begin()->FFp(0) << std::endl;
    }
     //    if ( m.face.begin()->FFp(0) != &*m.face.begin()) {
     
}

int main(int argc, char** argv) {
    std::cout << "%SUITE_STARTING% test_help_functions" << std::endl;
    std::cout << "%SUITE_STARTED%" << std::endl;

    std::cout << "%TEST_STARTED% testPairToBitQuad (test_help_functions)" << std::endl;
    testPairToBitQuad(); //ok
    std::cout << "%TEST_FINISHED% time=0 testPairToBitQuad (test_help_functions)" << std::endl;

    std::cout << "%TEST_STARTED% testPair (test_help_functions)" << std::endl;
    testPair(); //ok
    std::cout << "%TEST_FINISHED% time=0 testPair (test_help_functions)" << std::endl;

    std::cout << "%TEST_STARTED% testCreateTestMesh (test_help_functions)" << std::endl;
    testCountTriAndQuad(); //ok
    std::cout << "%TEST_FINISHED% time=0 testCreateTestMesh (test_help_functions)" << std::endl;

    std::cout << "%TEST_STARTED% testCountTriAndQuad (test_help_functions)" << std::endl;
    testCreateTestMesh(); //ok
    std::cout << "%TEST_FINISHED% time=0 testCountTriAndQuad (test_help_functions)" << std::endl;

    std::cout << "%SUITE_FINISHED% time=0" << std::endl;

    return (EXIT_SUCCESS);
}

