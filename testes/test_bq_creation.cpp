/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   test_bq_creation.cpp
 * Author: eu
 *
 * Created on December 5, 2019, 5:52 PM
 */

#include <stdlib.h>
#include <iostream>
#include "help_functions_for_tests.h" //meshclass
#include <vcg/complex/algorithms/bitquad_creation.h>
/*
 * Simple C++ Test Suite
 */



void testMakeDominantByTriOrder() {
    MeshType m;
    int quality=1;
    ScalarType limit=0;
    bool descedent = false;
    bool testConvex = false;
    
    CreateTestMesh(m,2,true,true,false);
    int result = vcg::tri::BitQuadCreation<MeshType>::MakeDominantByTriOrder(m, quality,limit, descedent,testConvex);
    if (result != vcg::tri::BitQuad<MeshType>::CountFauxs(m)) {
        std::cout << "%TEST_FAILED% time=0 testname=testMakeDominantByTriOrder (test_bq_creation) message=falha na contageem m2 n=" << result<<std::endl;
    }
}

void testMakeDominantByOrder() {
    MeshType m;
    int quality;
    ScalarType limit;
    bool descedent = true;
    bool testConvex = false;
    
    for (int malha=0;malha<7;malha++){
      CreateTestMesh(m,malha,false, true,false);
      for (quality = 0; quality < 3; quality++){
        vcg::tri::BitQuadCreation<MeshType>::MakeBitTriOnly(m);
        int result = vcg::tri::BitQuadCreation<MeshType>::MakeDominantByOrder(m, quality, limit, descedent, testConvex);
        if (result != vcg::tri::BitQuad<MeshType>::CountFauxs(m)) {
            std::cout << "%TEST_FAILED% time=0 testname=testMakeDominantByOrder (test_bq_creation) message=falha nacontagem q"<<quality<<" m"<<malha<<" qtd = " << result << std::endl;
        }
      }
    }
    
}

int main(int argc, char** argv) {
    std::cout << "%SUITE_STARTING% test_bq_creation" << std::endl;
    std::cout << "%SUITE_STARTED%" << std::endl;

    
    std::cout << "%TEST_STARTED% testMakeDominantByTriOrder (test_bq_creation)" << std::endl;
    testMakeDominantByTriOrder();
    std::cout << "%TEST_FINISHED% time=0 testMakeDominantByTriOrder (test_bq_creation)" << std::endl;

    std::cout << "%TEST_STARTED% testMakeDominantByOrder (test_bq_creation)" << std::endl;
    testMakeDominantByOrder();
    std::cout << "%TEST_FINISHED% time=0 testMakeDominantByOrder (test_bq_creation)" << std::endl;

    std::cout << "%SUITE_FINISHED% time=0" << std::endl;

    return (EXIT_SUCCESS);
}

