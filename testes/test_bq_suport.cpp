/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   test_bitquad_suport.cpp
 * Author: eu
 *
 * Created on December 5, 2019, 4:22 PM
 */

#include <stdlib.h>
#include <iostream>
#include "help_functions_for_tests.h"
#include "vcg/complex/algorithms/clean.h" //meshclass
#include <vcg/complex/algorithms/bitquad_support.h>
//#include <vcg/complex/algorithms/bitquad_creation.h> //PairToBitQuad
/*
 * Simple C++ Test Suite
 */

void testUpdateValencyInQualityWithBorders()
{
  MeshType m,m2;
  CreateTestMesh(m,1,true,true,false);
  //vcg::tri::UpdateFlags<MeshType>::FaceBorderFromFF(m);
  
  vcg::tri::BitQuad<MeshType>::UpdateValencyInQualityWithBorders(m);
  
  VertexIterator vi = m.vert.begin();
  if (vi->Q() != 2) {
      std::cout << "%TEST_FAILED% time=0 testname=testUpdateValencyInQualityWithBorders (test_bq_suport) message=m1 v0 valency 2!=" << vi->Q()<<std::endl;
    }
  vi++;
  if (vi->Q() != 3) {
      std::cout << "%TEST_FAILED% time=0 testname=testUpdateValencyInQualityWithBorders (test_bq_suport) message=m1 v1 valency 3!=" << vi->Q()<<std::endl;
    }
  
  PairToBitQuad<false>(&*(m.face.begin()),1);
  vcg::tri::BitQuad<MeshType>::UpdateValencyInQualityWithBorders(m);

  for (VertexIterator vi = m.vert.begin();  vi!=m.vert.end(); vi++){
    if (vi->Q() != 2) {
        std::cout << "%TEST_FAILED% time=0 testname=testUpdateValencyInQualityWithBorders (test_bq_suport) message=m1 v1 valency 2!=" << vi->Q()<<std::endl;
      }
  }
  
  CreateTestMesh(m,4,false,true,false);
  CreateTestMesh(m2,4,false,true,false);

  vcg::tri::BitQuadCreation<MeshType>::MakeDominant(m,2);
  vcg::tri::BitQuadCreation<MeshType>::MakeDominant(m2,2);
  //cout <<vcg::tri::BitQuad<MeshType>::CountFauxs(m)<<endl;
  
  vcg::tri::BitQuad<MeshType>::UpdateValencyInQualityWithBorders(m);
  vcg::tri::BitQuad<MeshType>::UpdateValencyInQualityWithBorders(m2);
  
  VertexIterator vi2 = m2.vert.begin();
  for (VertexIterator vi = m.vert.begin() ;  vi!=m.vert.end(); vi++,vi2++){
    //cout<<vi->Q()<<" "<<vi2->Q()<<" ";
    if (vi->Q() != vi2->Q()) {
        std::cout << "%TEST_FAILED% time=0 testname=testUpdateValencyInQualityWithBorders (test_bq_suport) message=m4 valencys of closed surface disferent " << vi->Q()<<std::endl;
      }
  }

  cout <<endl;
  
  //para superfícies fechadas comparar com UpdateValencyInQuality
}
void testMeasureQuality() {
    MeshType m;
    int quality=1;
    
    CreateTestMesh(m,1,false,true,false);
    vcg::tri::BitQuadCreation<MeshType>::MakeDominant(m,1);
    cout<<m.face.begin()->IsF(1)<<endl;
    
    ScalarType result = vcg::tri::BitQuad<MeshType>::MeasureQuality(m,quality);
    cout << result<<" "<<m.face.begin()->Q()<<" " <<m.face.begin()->FFp(1)->Q()<<endl;
    if (false) {
      std::cout << "%TEST_FAILED% time=0 testname=testMeasureQuality (test_bq_suport) message=error message sample" << std::endl;
    }
}

void testQuadQualitys() {
    MeshType m;
    CreateTestMesh(m,3,false, true, false);
    FaceType* f = &*m.face.begin();
    int version=0;
    int edgeInd=2;
    ScalarType result = vcg::tri::BitQuad<MeshType>::QuadQualitys(f, edgeInd, version);
    if (result != 1) {
        std::cout << "%TEST_FAILED% time=0 testname=testQuadQualitys (test_bitquad_suport) message=falha version0 no quadrado=" << result << std::endl;
    }
    
    version=2;
    edgeInd=2;
    result = vcg::tri::BitQuad<MeshType>::QuadQualitys(f, edgeInd, version);
    if (result != 1) {
        std::cout << "%TEST_FAILED% time=0 testname=testQuadQualitys (test_bitquad_suport) message=falha version2 no quadrado=" << result << std::endl;
    }
    
    CreateTestMesh(m,1,false, true, false);
    f = &*m.face.begin();
    version=0;
    edgeInd=1;
    result = vcg::tri::BitQuad<MeshType>::QuadQualitys(f, edgeInd, version);
    if ( result - 0.3964466 > TOL &&  0.3964467 - result > 0 ) {
        std::cout << "%TEST_FAILED% time=0 testname=testQuadQualitys (test_bitquad_suport) message=falha version0 no quad 180º plano=" << result << std::endl;
    }
    
    version=2;
    edgeInd=1;
    result = vcg::tri::BitQuad<MeshType>::QuadQualitys(f, edgeInd, version);
    if ( result - 0.3964466 > TOL &&  0.3964467 - result > 0 ) {
        std::cout << "%TEST_FAILED% time=0 testname=testQuadQualitys (test_bitquad_suport) message=falha version1 no quad 180º plano=" << result << std::endl;
    }
    
    version=1;
    edgeInd=0;
    result = vcg::tri::BitQuad<MeshType>::QuadQualitys(f, edgeInd, version);
    if ( result != 2) {
        std::cout << "%TEST_FAILED% time=0 testname=testQuadQualitys (test_bitquad_suport) message=falha version1 compr 2!=" << result << std::endl;
    }
    
    version=1;
    edgeInd=1;
    result = vcg::tri::BitQuad<MeshType>::QuadQualitys(f, edgeInd, version);
    if ( result - 1.41421356 > TOL && 1.41421357 - result > 0) {
        std::cout << "%TEST_FAILED% time=0 testname=testQuadQualitys (test_bitquad_suport) message=falha version1 compr sqrt(2) !=" << result << std::endl;
    }
    
    CreateTestMesh(m,4,false, true, false);
    f = &*m.face.begin();
    version=0;
    edgeInd=0;
    result = vcg::tri::BitQuad<MeshType>::QuadQualitys(f, edgeInd, version);
    if ( result - 0.6464466 > TOL &&  0.6464467 - result > 0 ) {
        std::cout << "%TEST_FAILED% time=0 testname=testQuadQualitys (test_bitquad_suport) message=falha version0 no quad 45º não plano=" << result << std::endl;
    }
    
    version=2;
    edgeInd=0;
    result = vcg::tri::BitQuad<MeshType>::QuadQualitys(f, edgeInd, version);
    if ( result - 0.7071067 > TOL &&  0.7071068 - result > 0 ) {
        std::cout << "%TEST_FAILED% time=0 testname=testQuadQualitys (test_bitquad_suport) message=falha version2 no quad 45º não plano=" << result << std::endl;
    }
}


void testCountFauxs() {
    MeshType m;
    //criar m
    int result = 0;
    result = vcg::tri::BitQuad<MeshType>::CountFauxs(m);
    if (result != 0) {
        std::cout << "%TEST_FAILED% time=0 testname=testCountFauxs (test_bitquad_suport) message=contando como faux em malha vazia = " << result << std::endl;
    }
    
    CreateTestMesh(m,3,false, true, false); //mude o primeiro false para ver as coordenadas
    result = vcg::tri::BitQuad<MeshType>::CountFauxs(m);
    if (result != 0) {
        std::cout << "%TEST_FAILED% time=0 testname=testCountFauxs (test_bitquad_suport) message=malha sem faux identificada com=" << result << std::endl;
    }
    
    PairToBitQuad<false>(&*m.face.begin(), 2);
    result = vcg::tri::BitQuad<MeshType>::CountFauxs(m);
    if (result != 1) {
        std::cout << "%TEST_FAILED% time=0 testname=testCountFauxs (test_bitquad_suport) message=contagem errada 1 != " << result << std::endl;
    }
    
    for (FaceIterator fi = m.face.begin();  fi!=m.face.end(); fi++) {
      fi->Q()=0.0;
      vcg::tri::BitQuadCreation<MeshType>::selectBestDiag<false>(&*fi); //precisa da componente  vcg::face::Qualityf
    }
    result = vcg::tri::BitQuad<MeshType>::CountFauxs(m);
    if (result != 3) {
        std::cout << "%TEST_FAILED% time=0 testname=testCountFauxs (test_bitquad_suport) message=contagem errada 3 != " << result << std::endl;
    }
    
    vcg::tri::BitQuadCreation<MeshType>::MakeBitTriOnly(m);
    m.face.begin()->FFp(1)->SetF(m.face.begin()->FFi(1));
    result = vcg::tri::BitQuad<MeshType>::CountFauxs(m);
    if (result != -1) {
        std::cout << "%TEST_FAILED% time=0 testname=testCountFauxs (test_bitquad_suport) message=falha ao identificar inconsistência= " << result << std::endl;
    }
}


void testTestQuadConvex() {
    MeshType m;
    CreateTestMesh(m,3,false, true, false); //mude o primeiro false para ver as coordenadas
    FaceType* f = &*m.face.begin();
    int ed=2;
    bool result = vcg::tri::BitQuad<MeshType>::TestQuadConvex(f, ed);
    if (!result) {
        std::cout << "%TEST_FAILED% time=0 testname=testTestQuadConvex (test_bitquad_suport) message=falha testando quadrado" << std::endl;
    }
    
    f++;f++;f++;f++; //penúltima face com última. supondo vetor de faces igual sequência de definição das faces.
    result = vcg::tri::BitQuad<MeshType>::TestQuadConvex(f, 0);
    if (result) {
        std::cout << "%TEST_FAILED% time=0 testname=testTestQuadConvex (test_bitquad_suport) message=falha em não convexo" << std::endl;
    }
    
    CreateTestMesh(m,1,false, true, false); //mude o primeiro false para ver as coordenadas
    f = &*m.face.begin();
    ed=1;
    result = vcg::tri::BitQuad<MeshType>::TestQuadConvex(f, ed);
    if (!result) {
        std::cout << "%TEST_FAILED% time=0 testname=testTestQuadConvex (test_bitquad_suport) message=falha em pi" << std::endl;
    }
}

int main(int argc, char** argv) {
    std::cout << "%SUITE_STARTING% test_bitquad_suport" << std::endl;
    std::cout << "%SUITE_STARTED%" << std::endl;
    
    std::cout << "%TEST_STARTED% testQuadQualitys (test_bitquad_suport)" << std::endl;
    testUpdateValencyInQualityWithBorders(); //OK
    std::cout << "%TEST_FINISHED% time=0 testQuadQualitys (test_bitquad_suport)" << std::endl;
    
    std::cout << "%TEST_STARTED% testQuadQualitys (test_bitquad_suport)" << std::endl;
    testMeasureQuality(); //OK
    std::cout << "%TEST_FINISHED% time=0 testQuadQualitys (test_bitquad_suport)" << std::endl;
    
    std::cout << "%TEST_STARTED% testQuadQualitys (test_bitquad_suport)" << std::endl;
    testQuadQualitys(); //OK
    std::cout << "%TEST_FINISHED% time=0 testQuadQualitys (test_bitquad_suport)" << std::endl;

    std::cout << "%TEST_STARTED% testCountFauxs (test_bitquad_suport)" << std::endl;
    testCountFauxs(); //ok
    std::cout << "%TEST_FINISHED% time=0 testCountFauxs (test_bitquad_suport)" << std::endl;

    std::cout << "%TEST_STARTED% testTestQuadConvex (test_bitquad_suport)" << std::endl;
    testTestQuadConvex(); //ok
    std::cout << "%TEST_FINISHED% time=0 testTestQuadConvex (test_bitquad_suport)" << std::endl;

    std::cout << "%SUITE_FINISHED% time=0" << std::endl;

    return (EXIT_SUCCESS);
}

