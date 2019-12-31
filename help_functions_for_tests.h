/*
 * 
 */

#ifndef TEST_HELP_FUNCTIONS_H
#define TEST_HELP_FUNCTIONS_H

#include "meshclass.h" //#include <vector> #include <vcg/complex/complex.h>
#include <string>
#include <chrono>       //high_resolution_clock::now()
#include "vcg/complex/algorithms/bitquad_creation.h" // #include "meshclass.h"   #include <iostream> istriquadonly
#include <vcg/complex/algorithms/create/platonic.h> //BuildMeshFromCoordVectorIndexVector
#include <vcg/complex/algorithms/bitquad_support.h> //countfaux
#include <wrap/io_trimesh/export.h>
#include <wrap/io_trimesh/export_vmi.h>
#include <wrap/io_trimesh/import_vmi.h>
#include <wrap/io_trimesh/import_off.h>
#include <vcg/complex/algorithms/stat.h> //histogram
#include <vcg/complex/algorithms/clean.h> //IsBitTriOnly, CountNonManifoldVertexFF


#define TOL 0.0000001 //tolerância para testar igualdades

/**
 aplica algoritmos de emparelhamento
 0 - nenhum
 1 - MakeDominantByOrder(m,0); //maior qualidade orto planar
 2 - MakeDominantByOrder(m,1); //Luis Velho
 3 - MakeDominantByOrder(m,2); //Maior Qualidade Orto Planar
 4 - MakeDominantByTriOrder(m,1); //Menor Tri Aspect Ratio Maior Aresta
 5 - MakeDominantByTriOrder(m,1,0,false,true); //Menor Tri Aspect Ratio Maior Aresta com TestConvex
 6 - MakeDominant(m,0); //linear level 0 - 1 iteração
 7 - MakeDominant(m,2); //linear level 2 - 4 iterações
 8 - MakeDominant(m,1); //linear level 1 - 3 iterações
 9 - MakeDominant(m,3); //linear level 3 - 2 iterações
 */
int ApplyPairing(MeshType& m, int nalg);

/**
 copia de vcg::tri::Stat<MeshType>::ComputePerFaceQualityHistogram, mas alterando a inicialização  SetRangeSetRange( minmax.first-0.01, minmax.second+0.1, HistSize );
 */
void ComputePerFaceQualityHistogram2( MeshType & m, vcg::Histogram<ScalarType> &h, bool selectionOnly=false,int HistSize=10);

/**
   cópia de vcg::tri::Stat<MeshType>::ComputePerVertexQualityHistogram, , mas alterando a inicialização  SetRange (SetRange( 0,9, 8);
 */
void ComputePerVertexQualityHistogram2( MeshType & m, vcg::Histogram<ScalarType> &h, bool selectionOnly = false, int HistSize=10 );

/**
 apply refinament and calculate face and vertex qualitys. Save in file with identifier 
 */
void GenerateHistogramsData (MeshType& m, string identifier, int quality, vector<vcg::Histogram<ScalarType>>& histFaces, vector<vcg::Histogram<ScalarType>>& histVerts);

/**
 save 2 files with qualitys of all faces and vertex
 */
void SaveQualitys(MeshType& m, string nameFileFaces="face qualitys", string nameFileVert="vertex qualitys");

/**
 Make the no border edge false for f and FFp.
 If edge is a border, return false.
 */
bool Pair(FaceType* f, int edge);

/**
 Save in parameters the number of quads and tris 
 tests if IsTriQuadOnly, if false return -2.
 if mesh has inconsistency return -1.
 Otherwise return number of quads
 */
 int CountTriAndQuad(MeshType& m, int& nQuads, int& nTris);

/**
 mostra as coordenas
 add param para mostrar lista de vértices de cada face tbm
 */
void PrintVertexCoord(const MeshType& m);

 /** Create a simple triangle mesh using just a vector of coords and a vector of indexes
 *  retorna a descrição da malha
 *  "Mesh tipo 1(Lista de vértices não convexo:4v2f)"
 *  "Mesh tipo 2(Lista de vértices não convexo e impar:5v3f)"
 *  "Mesh tipo 3(Lista de vértices:7v6f)"
 *  "Mesh tipo 4(Lista de vértices sem borda:6v8f)"
 *  "Mesh tipo 5(Tetrahedron)"
 *  "Mesh tipo 6(Octahedron)"  
 */
std::string CreateMesh(MeshType& m, int i);

/** Create a simple triangle mesh using just a vector of coords and a vector of indexes
 *  aplica alguns updates opcionalmente
 *  retorna a descrição da malha
 *  "Mesh tipo 1(Lista de vértices não convexo:4v2f)"
 *  "Mesh tipo 2(Lista de vértices não convexo e impar:5v3f)"
 *  "Mesh tipo 3(Lista de vértices:7v6f)"
 *  "Mesh tipo 4(Lista de vértices sem borda:6v8f)"
 *  "Mesh tipo 5(Tetrahedron)"
 *  "Mesh tipo 6(Octahedron)"  
 */
std::string CreateTestMesh(MeshType& m, int nMesh, bool showCoord, bool updateFFAdj, bool updateNormalFace);


/**
 * save the mesh m on a file with name "name[.format]"
 * @param m
 * @param name don't use extension
 * @param type =1. Define format.
 * \n1 .off com mask polygonal
 * \n2 .off 
 * \n3 .vmi (dump format)
 * \n4 .off with polygonal mask and it saves in order
 */
void SaveMesh(MeshType& m, std::string name, int type=4);

/** helper function:
// Create a bitQuad. given a triangle, merge it with its neightboord to form a quad through whichEdge
// return true if do it, false if edge is a border or ffp IsAnyF
// override to clear fauxs of both faces
*/
template <bool override>
bool PairToBitQuad(FaceType *fi, int whichEdge){
  if (vcg::face::IsBorder(*fi, whichEdge)) return false; // 
  if (override) {
    // clear any faux edge of the other face
    for (int k=0; k<3; k++){
      if (fi->FFp(whichEdge)->IsF(k)) {
        fi->FFp(whichEdge)->ClearF(k);
        fi->FFp(whichEdge)->FFp(k)->ClearF( fi->FFp(whichEdge)->FFi(k) ); // other face's ex-buddy is now single and sad :(
      }
    }
    // clear all faux edges of this face...
    for (int k=0; k<3; k++){
      if (fi->IsF(k)) {
        fi->ClearF(k);
        fi->FFp(k)->ClearF( fi->FFi(k) );
      }
    }
  }
  if (fi->FFp(whichEdge)->IsAnyF()) return false;
  fi->SetF(whichEdge);
  fi->FFp(whichEdge)->SetF( fi->FFi(whichEdge) );
  return true;
}
#endif /* TEST_HELP_FUNCTIONS_H */

